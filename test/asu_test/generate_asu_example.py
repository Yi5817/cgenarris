"""
Example of asymmetric-unit (ASU) generation.
"""

from __future__ import annotations

import os
import sys

import numpy as np
from mpi4py import MPI

import ase.io

# Make the in-place-built extension importable from this directory.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..", "..", "src")))

import pygenarris_mpi as pg_mpi


def molecules_to_arrays(molecules):
    """
    Marshal a list of ASE Atoms into the flat arrays the C function expects:
    concatenated (n_atoms x 3) positions, a two-char-per-atom species string
    ("C ", "Cl", ...), and the per-molecule-type atom counts.
    """
    positions = np.concatenate(
        [np.asarray(m.get_positions(), dtype=np.float64) for m in molecules], axis=0
    )
    species = "".join(s.ljust(2) for m in molecules for s in m.get_chemical_symbols())
    n_atoms_per_mol = np.array([len(m) for m in molecules], dtype=np.int32)
    return np.ascontiguousarray(positions), species, n_atoms_per_mol


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
    output_file = os.path.join(run_dir, "asu.extxyz")

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
        100000,      # max_attempts
        42,          # random_seed
        output_file,
        comm,
    )

    if comm.rank == 0:
        structures = ase.io.read(output_file, index=":", format="extxyz")
        print(f"\nrequested {num_structures}, generated {n_generated}, "
              f"read back {len(structures)} structures from {output_file}")
        first = structures[0]
        print("first structure:")
        print(f"  n_atoms      = {len(first)}")
        print(f"  pbc          = {tuple(bool(p) for p in first.pbc)}")
        print(f"  formula      = {first.get_chemical_formula()}")
        print(f"  info.sr      = {first.info.get('sr')}")
        print(f"  info.n_mols  = {first.info.get('n_mols')}")
        print(f"  info.stoic   = {first.info.get('stoic')}")


if __name__ == "__main__":
    main()
