"""
Tests for the Python interface (pygenarris_mpi) of cgenarris.

Covers both entry points, crystal and asymmetric-unit generation, and checks
that their output files share one block format readable by the same parser.

Run on one rank:   pytest tests/test_python_api.py
Run on two ranks:  mpirun -n 2 python -m pytest tests/test_python_api.py
Requires the in-place build: cd python && python setup.py build_ext --inplace
"""

from __future__ import annotations

import os
import shutil
import sys

import ase.io
import numpy as np
import pytest
from ase import Atoms
from ase.data import atomic_numbers, vdw_radii
from mpi4py import MPI

_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, "data")
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..", "python")))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..", "examples", "asu")))

import pygenarris_mpi as pg_mpi  # noqa: E402
from generate_asu_example import molecules_to_arrays, read_geometry_out  # noqa: E402

SR_MIN, SR_MAX = 0.75, 1.30
N_ASU = 12


@pytest.fixture(scope="module")
def comm():
    return MPI.COMM_WORLD


@pytest.fixture(scope="module")
def molecules():
    return [
        ase.io.read(os.path.join(_DATA, "asu", "geometry_0.in"), format="aims"),
        ase.io.read(os.path.join(_DATA, "asu", "geometry_1.in"), format="aims"),
    ]


@pytest.fixture
def run_dir(tmp_path_factory, comm):
    """
    One scratch directory shared by all ranks (rank 0 creates it).
    """
    path = str(tmp_path_factory.mktemp("run")) if comm.rank == 0 else None
    return comm.bcast(path, root=0)


# ---------------------------------------------------------------------------
# Crystal generation


def test_crystal_round_trip(run_dir, comm, monkeypatch):
    geometry_in = os.path.join(_DATA, "regression1", "geometry.in")
    if comm.rank == 0:
        shutil.copy(geometry_in, run_dir)
    comm.Barrier()
    monkeypatch.chdir(run_dir)  # the crystal generator reads/writes in cwd
    mol: Atoms = ase.io.read(geometry_in, format="aims")
    z, volume_mean, volume_std, sr = 2, 600.0, 40.0, 0.85

    # Pairwise cutoff sr * (r_i + r_j) for z copies of the molecule.
    symbols = mol.get_chemical_symbols()
    radii = np.tile([vdw_radii[atomic_numbers[s]] for s in symbols], z)
    cutoff = sr * (radii[:, None] + radii[None, :])
    cutoff = np.ascontiguousarray(cutoff, dtype=np.float32)

    pg_mpi.mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(
        cutoff,
        1,            # structures per space group
        z,
        volume_mean,
        volume_std,
        0.1,          # tolerance
        2000,         # max attempts per space group
        "uniform",
        10000,        # volume attempts
        223344,       # random seed
        0.4,          # lattice norm deviation
        8.0,          # lattice angle std
        comm,
    )
    comm.Barrier()

    crystals = read_geometry_out(os.path.join(run_dir, "geometry.out"))
    assert len(crystals) > 0
    for i, xtal in enumerate(crystals, start=1):
        assert xtal.pbc.all()
        assert len(xtal) == z * len(mol)
        assert xtal.get_chemical_symbols() == mol.get_chemical_symbols() * z
        assert xtal.info["structure_number"] == i
        assert xtal.info["Z"] == z
        assert xtal.info["number_of_atoms_in_molecule"] == len(mol)
        assert abs(xtal.get_volume() - volume_mean) < 5 * volume_std


# ---------------------------------------------------------------------------
# Asymmetric-unit generation


def _generate_asus(molecules, stoichiometry, output_file, comm, **kwargs):
    positions, species, n_atoms_per_mol = molecules_to_arrays(molecules)
    args = {
        "num_structures": N_ASU,
        "sr_min": SR_MIN,
        "sr_max": SR_MAX,
        "max_attempts": 100000,
        "random_seed": 7,
    }
    args.update(kwargs)
    return pg_mpi.generate_asymmetric_units(
        positions,
        species,
        n_atoms_per_mol,
        np.asarray(stoichiometry, dtype=np.int32),
        args["num_structures"],
        args["sr_min"],
        args["sr_max"],
        args["max_attempts"],
        args["random_seed"],
        output_file,
        comm,
    )


def test_asu_round_trip_two_components(molecules, run_dir, comm):
    output_file = os.path.join(run_dir, "asu.out")
    assert _generate_asus(molecules, [1, 1], output_file, comm) == N_ASU
    comm.Barrier()

    asus = read_geometry_out(output_file)
    assert len(asus) == N_ASU
    n0, n1 = len(molecules[0]), len(molecules[1])
    for i, unit in enumerate(asus, start=1):
        assert len(unit) == n0 + n1
        assert not unit.pbc.any()
        # Numbering is consecutive across ranks in the merged file.
        assert unit.info["structure_number"] == i
        assert unit.info["number_of_atoms"] == n0 + n1
        assert unit.info["number_of_molecules"] == 2
        assert unit.info["number_of_molecule_types"] == 2
        assert SR_MIN < unit.info["sr"] < SR_MAX
        np.testing.assert_array_equal(unit.info["stoichiometry"], [1, 1])
        np.testing.assert_array_equal(unit.info["molecule_types"], [0, 1])
        np.testing.assert_array_equal(unit.info["molecule_index"], [0, n0])
        np.testing.assert_array_equal(
            unit.info["number_of_atoms_in_molecule_type"], [n0, n1]
        )
        assert unit.get_chemical_symbols() == (
            molecules[0].get_chemical_symbols() + molecules[1].get_chemical_symbols()
        )
        # Rigid bodies: intramolecular distances are preserved.
        for mol, start in zip(molecules, unit.info["molecule_index"]):
            block = unit[start : start + len(mol)]
            np.testing.assert_allclose(
                block.get_all_distances(), mol.get_all_distances(), atol=2e-3
            )
        # Recentred at the origin.
        np.testing.assert_allclose(unit.get_positions().mean(axis=0), 0, atol=1e-4)


def test_asu_stoichiometry_two_to_one(molecules, run_dir, comm):
    output_file = os.path.join(run_dir, "asu.out")
    assert _generate_asus(molecules, [2, 1], output_file, comm, num_structures=4) == 4
    comm.Barrier()
    unit = read_geometry_out(output_file)[0]
    assert len(unit) == 2 * len(molecules[0]) + len(molecules[1])
    np.testing.assert_array_equal(unit.info["molecule_types"], [0, 0, 1])
    np.testing.assert_array_equal(unit.info["stoichiometry"], [2, 1])


def test_asu_seed_reproducibility(molecules, run_dir, comm):
    output_file = os.path.join(run_dir, "asu.out")
    runs = []
    for _ in range(2):
        _generate_asus(
            molecules, [1, 1], output_file, comm, num_structures=6, random_seed=3
        )
        comm.Barrier()
        runs.append(read_geometry_out(output_file))
        comm.Barrier()
    for a, b in zip(*runs):
        np.testing.assert_array_equal(a.get_positions(), b.get_positions())


@pytest.mark.skipif(not os.path.exists("/dev/full"), reason="requires /dev/full")
@pytest.mark.parametrize("fail_merge", [False, True])
def test_asu_write_failure(
    molecules: list[Atoms], run_dir: str, comm: MPI.Comm, fail_merge: bool
) -> None:
    """
    Report shard and merge write failures on every rank; retain merge inputs.
    """
    output_file = os.path.join(run_dir, "asu.out")
    if comm.rank == 0:
        target = output_file if fail_merge else output_file + ".rank0"
        os.symlink("/dev/full", target)
    comm.Barrier()
    result = _generate_asus(
        molecules, [1, 1], output_file, comm, num_structures=comm.size
    )
    assert comm.allgather(result) == [-1] * comm.size
    if fail_merge:
        shard = f"{output_file}.rank{comm.rank}"
        assert len(read_geometry_out(shard)) == 1


def test_asu_invalid_input_returns_minus_one(molecules, run_dir, comm):
    output_file = os.path.join(run_dir, "asu.out")
    assert (
        _generate_asus(molecules, [1, 1], output_file, comm, sr_min=1.3, sr_max=0.75)
        == -1
    )
    assert _generate_asus(molecules, [1, 1], output_file, comm, max_attempts=0) == -1
    assert _generate_asus(molecules, [1], output_file, comm) == -1  # too short
    assert _generate_asus(molecules, [0, 1], output_file, comm) == -1
    # A single molecule has no pair to score.
    assert _generate_asus(molecules[:1], [1], output_file, comm) == -1


def test_asu_bad_array_shapes_raise(molecules, run_dir, comm):
    output_file = os.path.join(run_dir, "asu.out")
    positions, species, n_atoms_per_mol = molecules_to_arrays(molecules)
    stoic = np.array([1, 1], dtype=np.int32)
    with pytest.raises(TypeError):
        pg_mpi.generate_asymmetric_units(
            positions.ravel(), species, n_atoms_per_mol, stoic,
            N_ASU, SR_MIN, SR_MAX, 1000, 1, output_file, comm,
        )
    with pytest.raises(TypeError):
        pg_mpi.generate_asymmetric_units(
            positions, species, n_atoms_per_mol.astype(np.float64), stoic,
            N_ASU, SR_MIN, SR_MAX, 1000, 1, output_file, comm,
        )
    # Inconsistent species length or atom counts are caught in C rather than
    # read out of bounds.
    assert pg_mpi.generate_asymmetric_units(
        positions, species[:-2], n_atoms_per_mol, stoic,
        N_ASU, SR_MIN, SR_MAX, 1000, 1, output_file, comm,
    ) == -1
    assert pg_mpi.generate_asymmetric_units(
        positions, species, n_atoms_per_mol + 1, stoic,
        N_ASU, SR_MIN, SR_MAX, 1000, 1, output_file, comm,
    ) == -1
    # A negative count must not pass the sum check by offsetting another entry.
    offset = np.array([-1, 1], dtype=np.int32) * (n_atoms_per_mol[0] + 1)
    assert pg_mpi.generate_asymmetric_units(
        positions, species, n_atoms_per_mol + offset, stoic,
        N_ASU, SR_MIN, SR_MAX, 1000, 1, output_file, comm,
    ) == -1
