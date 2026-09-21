cmake_minimum_required(VERSION  3.10.2)
# The VERSION file is the single source of truth; release tags mirror it.
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION CGENARRIS_VERSION)
project("cgenarris" VERSION ${CGENARRIS_VERSION} LANGUAGES C)
add_definitions(-DCGENARRIS_VERSION="${PROJECT_VERSION}")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_C_FLAGS "-Wall -Wpedantic -O3 -std=gnu99")

include_directories("${PROJECT_SOURCE_DIR}/src")
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

# Layer groups
src/lattice_generator_layer.c

# spglib sources
src/spglib_src/arithmetic.c
src/spglib_src/cell.c
src/spglib_src/debug.c
src/spglib_src/delaunay.c
src/spglib_src/determination.c
src/spglib_src/hall_symbol.c
src/spglib_src/kgrid.c
src/spglib_src/kpoint.c
src/spglib_src/mathfunc.c
src/spglib_src/niggli.c
src/spglib_src/overlap.c
src/spglib_src/pointgroup.c
src/spglib_src/primitive.c
src/spglib_src/refinement.c
src/spglib_src/sitesym_database.c
src/spglib_src/site_symmetry.c
src/spglib_src/spacegroup.c
src/spglib_src/spg_database.c
src/spglib_src/spglib.c
src/spglib_src/spglib_f.c
src/spglib_src/spin.c
src/spglib_src/symmetry.c
src/spglib_src/tetrahedron_method.c
)

add_executable(${PROJECT_NAME} ${SOURCES} src/cgenarris_mpi.c)
add_executable(cgenarris_test ${SOURCES} src/test_main.c)
add_executable(test_check_structure ${SOURCES} src/test_check_structure.c)
target_link_libraries(${PROJECT_NAME} m)
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})
target_link_libraries(cgenarris_test m)
target_link_libraries(cgenarris_test ${MPI_LIBRARIES})
target_link_libraries(test_check_structure m)
target_link_libraries(test_check_structure ${MPI_LIBRARIES})
add_executable(test_lattice_generator src/lattice_generator.c src/randomgen.c
               src/test_lattice_generator.c)
target_link_libraries(test_lattice_generator m)

enable_testing()

set(REGRESSION_DIR ${CMAKE_BINARY_DIR}/regression_test1)
file(MAKE_DIRECTORY ${REGRESSION_DIR})
configure_file(${PROJECT_SOURCE_DIR}/test/regression_test1/geometry.in
               ${REGRESSION_DIR}/geometry.in COPYONLY)
add_test(NAME cgenarris_test
         COMMAND cgenarris_test
         WORKING_DIRECTORY ${REGRESSION_DIR})

# Bad-structure rejection test: a Z=4 overlapping structure produced when
# the user requested Z=2. The vdw matrix is sized for Z=2, so the check
# must reject the mismatched structure. Expected exit code: 1 (FAIL).
add_test(NAME check_structure_rejects_bad
         COMMAND test_check_structure
            ${PROJECT_SOURCE_DIR}/test/check_structure_test/geometry.in
            ${PROJECT_SOURCE_DIR}/test/check_structure_test/cutoff_matrix.txt
            ${PROJECT_SOURCE_DIR}/test/check_structure_test/geometry.out)
set_tests_properties(check_structure_rejects_bad PROPERTIES WILL_FAIL TRUE)

# Lattice generator: volume must equal the target and cell parameters must
# obey the constraints of each crystal system.
add_test(NAME lattice_generator COMMAND test_lattice_generator)
