import numpy as np

from . import rigid_press as rp
from gnrs.cgenarris import pygenarris_mpi as pg


def optimize_structure(struct, Z, cutoff_mat, spg=0, max_iter=400):
    """
    Optimizes the geometry of structure using rigid press
    """

    # print(struct.get_positions())
    # Create a pygenarris crystal
    xtal = pg.crystal()
    pos = struct.get_positions()
    xc, yc, zc = pos[:, 0], pos[:, 1], pos[:, 2]
    species = struct.get_chemical_symbols()
    lattice = struct.cell.array
    species_str = ""
    for i, at in enumerate(species):
        if len(at) == 1:
            species_str += at + " "
        else:
            species_str += at
    # spg = get_spacegroup(struct).no
    n_atoms = len(xc)
    pg.create_crystal_from_array(
        xtal, lattice, xc, yc, zc, species_str, n_atoms, Z, spg
    )

    # Create settings structure
    opt_set = rp.Opt_settings()
    opt_set.spg = spg
    opt_set.max_iteration = max_iter
    cutoff_mat = cutoff_mat.flatten()
    status = rp.optimize_crystal(xtal, cutoff_mat, opt_set)
    # status = 0

    # Update ASE structure
    lattice = np.ascontiguousarray(lattice)
    xc, yc, zc = (
        np.ascontiguousarray(xc),
        np.ascontiguousarray(yc),
        np.ascontiguousarray(zc),
    )

    pg.create_array_from_crystal(
        xtal, lattice, xc, yc, zc, species_str, n_atoms, Z, spg
    )

    struct.set_cell(lattice)
    pos[:, 0], pos[:, 1], pos[:, 2] = xc, yc, zc
    struct.set_positions(pos)

    if status == 0:
        return True
    else:
        return False


def test_optimize():
    from ase.io import read

    cutoff_mat = np.loadtxt(
       "../sample_structures/Example4/cutoff_matrix.txt", dtype="float32"
    )
    struct = read("../sample_structures/Example4/geometry.in")
    print(struct.positions, struct.cell)
    optimize_structure(struct, 2, cutoff_mat, spg=2)
    print(struct.positions, struct.cell)


if __name__ == "__main__":
    test_optimize()
