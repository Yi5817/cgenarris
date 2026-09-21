"""
setup.py file
"""

import os

from setuptools import setup, Extension
from distutils import sysconfig

mpicompiler = "mpicc"

version_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "VERSION")
with open(version_file) as f:
    version = f.read().strip()

# These flags may conflict with other compilers
ccvars = sysconfig.get_config_vars()
key_list1 = [
    "BASECFLAGS",
    "CFLAGS",
    "OPT",
    "PY_CFLAGS",
    "CCSHARED",
    "CFLAGSFORSHARED",
    "LINKFORSHARED",
    "LIBS",
    "SHLIBS",
]
for key in key_list1:
    if key in ccvars:
        ccvars[key] = " "

key_list2 = ["CC", "LDSHARED"]
for key in key_list2:
    if key in ccvars:
        value = ccvars[key].split()
        value[0] = mpicompiler
        ccvars[key] = " ".join(value)

packages = ["numpy", "mpi4py"]
for package in packages:
    try:
        __import__(package)
    except ImportError:
        print("Please install", package)
        exit()

import mpi4py
import numpy

sources_spglib = [
    "arithmetic.c",
    "cell.c",
    "delaunay.c",
    "determination.c",
    "hall_symbol.c",
    "kgrid.c",
    "kpoint.c",
    "mathfunc.c",
    "niggli.c",
    "overlap.c",
    "pointgroup.c",
    "primitive.c",
    "refinement.c",
    "sitesym_database.c",
    "site_symmetry.c",
    "spacegroup.c",
    "spin.c",
    "spg_database.c",
    "spglib.c",
    "symmetry.c",
]

source_dir = "spglib_src"
include_dirs = [
    source_dir,
]
for i, s in enumerate(sources_spglib):
    sources_spglib[i] = "%s/%s" % (source_dir, s)

pygenarris_mpi = Extension(
    "_pygenarris_mpi",
    include_dirs=[numpy.get_include(), mpi4py.get_include(), "./"],
    sources=[
        "pygenarris_mpi.i",
        "pygenarris_mpi.c",
        "combinatorics.c",
        "molecule_placement.c",
        "algebra.c",
        "molecule_utils.c",
        "spg_generation.c",
        "lattice_generator.c",
        "crystal_utils.c",
        "check_structure.c",
        "read_input.c",
        "randomgen.c",
        "lattice_generator_layer.c",
        "pygenarris_mpi_utils.c",
        "asu_generation.c",
        "asu_utils.c",
    ]
    + sources_spglib,
    extra_compile_args=["-std=gnu99", "-fPIC", "-O3", "-Wno-error=int-conversion"],
    swig_opts=[
        "-I./",
        f"-I{mpi4py.get_include()}",
    ],
    define_macros=[("CGENARRIS_VERSION", '"%s"' % version)],
)

setup(
    name="pygenarris_mpi",
    version=version,
    author="Rithwik Tom, Yi Yang",
    maintainer="Yi Yang, Haoran Ni",
    maintainer_email="yiy5@andrew.cmu.edu",
    ext_modules=[pygenarris_mpi],
    py_modules=["pygenarris_mpi"],
)
