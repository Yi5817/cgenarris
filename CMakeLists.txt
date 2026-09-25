cmake_minimum_required(VERSION  3.10.2)
# The VERSION file is the single source of truth; release tags mirror it.
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION CGENARRIS_VERSION)
project("cgenarris" VERSION ${CGENARRIS_VERSION} LANGUAGES C)
add_definitions(-DCGENARRIS_VERSION="${PROJECT_VERSION}")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_C_FLAGS "-Wall -Wpedantic -O3 -std=gnu99")

# Layout
#   include/cgenarris/   public headers  (#include "cgenarris/crystal.h")
#   src/                 library sources and private headers
#   third_party/spglib/  vendored spglib
#   apps/                executables
#   tests/               ctest sources; fixtures in tests/data/
#   python/              SWIG bindings (built by python/setup.py, not CMake)
include_directories("${PROJECT_SOURCE_DIR}/include"
                    "${PROJECT_SOURCE_DIR}/src"
                    "${PROJECT_SOURCE_DIR}/third_party/spglib")
set(SOURCES
src/algebra.c
src/check_structure.c
src/combinatorics.c
src/crystal_utils.c
src/lattice_generator.c
src/molecule_placement.c
src/molecule_utils.c
src/pygenarris_mpi.c
src/pygenarris_mpi_utils.c
src/randomgen.c
src/read_input.c
src/spg_generation.c
src/asu_generation.c
src/asu_utils.c

# Layer groups
src/lattice_generator_layer.c

# spglib sources
third_party/spglib/arithmetic.c
third_party/spglib/cell.c
third_party/spglib/debug.c
third_party/spglib/delaunay.c
third_party/spglib/determination.c
third_party/spglib/hall_symbol.c
third_party/spglib/kgrid.c
third_party/spglib/kpoint.c
third_party/spglib/mathfunc.c
third_party/spglib/niggli.c
third_party/spglib/overlap.c
third_party/spglib/pointgroup.c
third_party/spglib/primitive.c
third_party/spglib/refinement.c
third_party/spglib/sitesym_database.c
third_party/spglib/site_symmetry.c
third_party/spglib/spacegroup.c
third_party/spglib/spg_database.c
third_party/spglib/spglib.c
third_party/spglib/spglib_f.c
third_party/spglib/spin.c
third_party/spglib/symmetry.c
third_party/spglib/tetrahedron_method.c
)

# One static library; every executable links it instead of recompiling SOURCES.
add_library(cgenarris_core STATIC ${SOURCES})
target_link_libraries(cgenarris_core m ${MPI_LIBRARIES})

add_executable(${PROJECT_NAME} apps/cgenarris_mpi.c)
target_link_libraries(${PROJECT_NAME} cgenarris_core)
add_executable(cgenarris_test tests/test_main.c)
target_link_libraries(cgenarris_test cgenarris_core)
add_executable(test_check_structure tests/test_check_structure.c)
target_link_libraries(test_check_structure cgenarris_core)
add_executable(test_lattice_generator src/lattice_generator.c src/randomgen.c
               tests/test_lattice_generator.c)
target_link_libraries(test_lattice_generator m)
add_executable(test_asu tests/test_asu.c)
target_link_libraries(test_asu cgenarris_core)

enable_testing()

set(TEST_DATA ${PROJECT_SOURCE_DIR}/tests/data)

set(REGRESSION_DIR ${CMAKE_BINARY_DIR}/regression_test1)
file(MAKE_DIRECTORY ${REGRESSION_DIR})
configure_file(${TEST_DATA}/regression1/geometry.in
               ${REGRESSION_DIR}/geometry.in COPYONLY)
add_test(NAME cgenarris_test
         COMMAND cgenarris_test
         WORKING_DIRECTORY ${REGRESSION_DIR})

# Bad-structure rejection test: a Z=4 overlapping structure produced when
# the user requested Z=2. The vdw matrix is sized for Z=2, so the check
# must reject the mismatched structure. Expected exit code: 1 (FAIL).
add_test(NAME check_structure_rejects_bad
         COMMAND test_check_structure
            ${TEST_DATA}/check_structure/geometry.in
            ${TEST_DATA}/check_structure/cutoff_matrix.txt
            ${TEST_DATA}/check_structure/geometry.out)
set_tests_properties(check_structure_rejects_bad PROPERTIES WILL_FAIL TRUE)

# Lattice generator: volume must equal the target and cell parameters must
# obey the constraints of each crystal system.
add_test(NAME lattice_generator COMMAND test_lattice_generator)

# Asymmetric units: unit checks plus an end-to-end run on the two JESMOK
# components; also exercises the control.in-driven executable.
set(ASU_TEST_DIR ${CMAKE_BINARY_DIR}/asu_test)
file(MAKE_DIRECTORY ${ASU_TEST_DIR})
foreach(f geometry_0.in geometry_1.in control.in)
    configure_file(${TEST_DATA}/asu/${f}
                   ${ASU_TEST_DIR}/${f} COPYONLY)
endforeach()
add_test(NAME asu_generation
         COMMAND test_asu
         WORKING_DIRECTORY ${ASU_TEST_DIR})
add_test(NAME cgenarris_asu_executable
         COMMAND ${PROJECT_NAME}
         WORKING_DIRECTORY ${ASU_TEST_DIR})
set_tests_properties(cgenarris_asu_executable PROPERTIES
                     PASS_REGULAR_EXPRESSION "Generated 20/20 asymmetric units")

add_executable(test_geometry_bugs tests/test_geometry_bugs.c)
target_link_libraries(test_geometry_bugs cgenarris_core)
add_test(NAME geometry_bugs COMMAND test_geometry_bugs)
add_test(NAME lattice_small_volumes COMMAND test_geometry_bugs small-volumes)
set_tests_properties(lattice_small_volumes PROPERTIES TIMEOUT 5)
