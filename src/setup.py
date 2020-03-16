
"""
setup.py file 
"""

from setuptools import setup, Extension
from distutils import sysconfig
import os

mpicompiler = 'mpicc'

# These flags may conflict with other compilers
ccvars = sysconfig.get_config_vars()
key_list1 = ['BASECFLAGS', 'CFLAGS', 'OPT', 'PY_CFLAGS',
            'CCSHARED', 'CFLAGSFORSHARED', 'LINKFORSHARED',
            'LIBS', 'SHLIBS']
for key in key_list1:
    if key in ccvars:
        ccvars[key] = ' '

key_list2 = ['CC', 'LDSHARED']
for key in key_list2:
    if key in ccvars:
        value = ccvars[key].split()
        value[0] = mpicompiler
        ccvars[key] = ' '.join(value)

packages = ['numpy', 'mpi4py']
for package in packages:
    try:
        __import__(package)
    except ImportError:
        print("Please install", package)
        exit()

import mpi4py
import numpy

sources_spglib = ['arithmetic.c',
           'cell.c',
           'delaunay.c',
           'determination.c',
           'hall_symbol.c',
           'kgrid.c',
           'kpoint.c',
           'mathfunc.c',
           'niggli.c',
           'overlap.c',
           'pointgroup.c',
           'primitive.c',
           'refinement.c',
           'sitesym_database.c',
           'site_symmetry.c',
           'spacegroup.c',
           'spin.c',
           'spg_database.c',
           'spglib.c',
           'symmetry.c']

source_dir = "spglib_src"
include_dirs = [source_dir, ]
for i, s in enumerate(sources_spglib):
    sources_spglib[i] = "%s/%s" % (source_dir, s)

pygenarris_mpi = Extension('_pygenarris_mpi',
                           include_dirs= ['./', numpy.get_include(), mpi4py.get_include()],
                           sources=['pygenarris_mpi.i', 'pygenarris_mpi.c', 'combinatorics.c', 'molecule_placement.c',
                           'algebra.c', 'molecule_utils.c','spg_generation.c', 'lattice_generator.c', 'crystal_utils.c',
                           'check_structure.c', 'read_input.c', 'randomgen.c']+sources_spglib,
                           extra_compile_args=["-std=gnu99", "-fPIC", "-O3"])

setup (name = 'pygenarris_mpi',
       version = '0.1',
       author      = "Rithwik Tom",
       description = """email:rtom@andrew.cmu.edu""",
       ext_modules = [pygenarris_mpi],
       py_modules = ["pygenarris_mpi"],
       )
