"""
Example of asymmetric-unit (ASU) generation from Python.

Run:  mpirun -n 2 python generate_asu_example.py
"""

from __future__ import annotations

import os
import sys

import ase.io
import numpy as np
from ase import Atoms
from mpi4py import MPI

# Make the in-place-built extension importable from this directory.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..", "..", "src")))

import pygenarris_mpi as pg_mpi  # noqa: E402


def molecules_to_arrays(molecules: list[Atoms]) -> tuple[np.ndarray, str, np.ndarray]:
    """
    Marshal a list of ASE Atoms into the flat arrays the C function expects.

    Returns:
        Concatenated (n_atoms, 3) float64 positions, a two-chars-per-atom
        species string ("C ", "Cl", ...), and the per-molecule-type atom counts.
    """
    positions = np.concatenate(
        [np.asarray(m.get_positions(), dtype=np.float64) for m in molecules], axis=0
    )
    species = "".join(s.ljust(2) for m in molecules for s in m.get_chemical_symbols())
    n_atoms_per_mol = np.array([len(m) for m in molecules], dtype=np.int32)
    return np.ascontiguousarray(positions), species, n_atoms_per_mol


def read_geometry_out(path: str) -> list[Atoms]:
    """
    Read a cgenarris output file (crystals or asymmetric units) into ASE Atoms.

    Each "####### BEGIN/END STRUCTURE #######" block becomes one Atoms object.
    "#key = value" lines land in Atoms.info (ints, floats, or int lists);
    blocks without lattice_vector lines are non-periodic.
    """
    structures: list[Atoms] = []
    with open(path) as f:
        blocks = f.read().split("#######  END  STRUCTURE #######")[:-1]
    for block in blocks:
        symbols, positions, cell, info = [], [], [], {}
        for line in block.splitlines():
            parts = line.split()
            if not parts or parts[0].startswith("#######"):
                continue
            if parts[0] == "atom":
                positions.append([float(v) for v in parts[1:4]])
                symbols.append(parts[4])
            elif parts[0] == "lattice_vector":
                cell.append([float(v) for v in parts[1:4]])
            elif parts[0].startswith("#") and len(parts) > 2 and parts[1] == "=":
                values = parts[2:]
                try:
                    parsed = [int(v) for v in values]
                    info[parts[0][1:]] = parsed[0] if len(parsed) == 1 else np.array(parsed)
                except ValueError:
                    try:
                        info[parts[0][1:]] = float(values[0])
                    except ValueError:
                        info[parts[0][1:]] = " ".join(values)
        atoms = Atoms(symbols=symbols, positions=positions)
        if cell:
            atoms.set_cell(cell)
            atoms.set_pbc(True)
        atoms.info.update(info)
        structures.append(atoms)
    return structures


def main() -> None:
    """
    Generate asymmetric units.
    """
    comm = MPI.COMM_WORLD

    molecules = [ase.io.read(os.path.join(_HERE, "molecule.xyz"))]  # one type ...
    stoichiometry = np.array([2], dtype=np.int32)                   # ... taken twice -> dimer

    run_dir = os.path.join(_HERE, "run")
    if comm.rank == 0:
        os.makedirs(run_dir, exist_ok=True)
    comm.Barrier()
    output_file = os.path.join(run_dir, "asu.out")

    positions, species, n_atoms_per_mol = molecules_to_arrays(molecules)

    num_structures = 20
    n_generated = pg_mpi.generate_asymmetric_units(
        positions,
        species,
        n_atoms_per_mol,
        stoichiometry,
        num_structures,
        0.75,        # sr_min
        1.30,        # sr_max
        100000,      # max_attempts per asymmetric unit
        42,          # random_seed
        output_file,
        comm,
    )

    if comm.rank == 0:
        asus = read_geometry_out(output_file)
        print(f"\nrequested {num_structures}, generated {n_generated}, "
              f"read back {len(asus)} asymmetric units from {output_file}")
        first = asus[0]
        print("first asymmetric unit:")
        print(f"  n_atoms             = {len(first)}")
        print(f"  pbc                 = {tuple(bool(p) for p in first.pbc)}")
        print(f"  formula             = {first.get_chemical_formula()}")
        print(f"  info.sr             = {first.info['sr']}")
        print(f"  info.stoichiometry  = {first.info['stoichiometry']}")
        print(f"  info.molecule_types = {first.info['molecule_types']}")
        print(f"  info.molecule_index = {first.info['molecule_index']}")


if __name__ == "__main__":
    main()
