<h1 align="center">cgenarris</h1>

<p align="center">
Random molecular crystal &amp; layer structure generator written in C with Python bindings.
</p>

<p align="center">
<a href="LICENSE"><img src="https://img.shields.io/badge/license-BSD--3--Clause-blue.svg" alt="License"></a>
</p>

## Overview

**cgenarris** is an MPI-parallel program that generates random molecular crystal structures based on space group symmetry. It supports all 230 space groups for 3D crystals and layer groups for 2D structures.

**Key capabilities:**

- Random crystal structure generation with van der Waals distance validation
- Layer crystal generation with substrate lattice constraints
- Rigid-press optimization for improved packing acceptance rates
- Python API via SWIG with NumPy and mpi4py integration

## Contributors

| Name | Contribution | GitHub |
|------|-------------|--------|
| Rithwik Tom | Original author | [@ritwit](https://github.com/ritwit) |
| Yi Yang | Developer | [@Yi5817](https://github.com/Yi5817) |
| Haoran Ni | Layer group generator | [@haoran-ni](https://github.com/haoran-ni) |

## License

cgenarris is available under the [BSD-3-Clause License](LICENSE).
