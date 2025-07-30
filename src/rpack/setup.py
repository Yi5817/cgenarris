"""
setup.py file
"""

from setuptools import setup, Extension, find_packages
from distutils import sysconfig
import os

packages = ["numpy"]
for package in packages:
    try:
        __import__(package)
    except ImportError:
        print("Please install", package)
        exit()

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

sources_rpress = [
    "rigid_press.i",
    "rigid_press.c",
    "symmetrization.c",
    "d_algebra.c",
]

spg_source_dir = "../spglib_src"
rpress_source_dir = "rigid_press"

for i, src in enumerate(sources_spglib):
    sources_spglib[i] = f"{spg_source_dir}/{src}"
for i, src in enumerate(sources_rpress):
    sources_rpress[i] = f"{rpress_source_dir}/{src}"


rigid_press = Extension(
    "rigid_press._rigid_press",
    include_dirs=["./", numpy.get_include()],
    sources=sources_rpress + sources_spglib,
    extra_compile_args=["-std=gnu99", "-fPIC", "-O3", "-DROPT_DEBUG"],
    libraries=["lapack", "blas"],
)

setup(
    name="rigid_press",
    version="0.1",
    author="Jonathan Moussa, Rithwik Tom",
    description="""email:rtom@andrew.cmu.edu""",
    ext_modules=[rigid_press],
    py_modules=["rigid_press"],
)
