cmake_minimum_required(VERSION 3.16)
project(PrimHyb)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-march=native -O2 -std=c++2a)
file(GLOB_RECURSE SRCS src/*.cpp)

include_directories(src)

# add_executable(nc_test exec/near_collision_test.cpp)
# target_sources(nc_test PRIVATE ${SRCS})

# add_executable(fcons_test exec/fixed_constraint_test.cpp)
# target_sources(fcons_test PRIVATE ${SRCS})

# add_executable(adm_test exec/admissible_test.cpp)
# target_sources(adm_test PRIVATE ${SRCS})

# add_executable(cons_test exec/constraint_test.cpp)
# target_sources(cons_test PRIVATE ${SRCS})

add_executable(unif_test exec/uniform_test.cpp)
target_sources(unif_test PRIVATE ${SRCS})

add_executable(cons_simulation exec/cons_simulation.cpp)
target_sources(cons_simulation PRIVATE ${SRCS})

add_executable(lsh_simulation exec/lsh_simulation.cpp)
target_sources(lsh_simulation PRIVATE ${SRCS})

# add_executable(var_test exec/variance_test.cpp)
# target_sources(var_test PRIVATE ${SRCS})

