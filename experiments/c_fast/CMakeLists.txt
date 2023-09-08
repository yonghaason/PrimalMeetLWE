cmake_minimum_required(VERSION 3.16)
project(PrimHyb)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-march=native -O2 -std=c++2a)
file(GLOB_RECURSE SRCS src/*.cpp)

include_directories(src)

add_executable(jensen_gap exec/jensen_gap.cpp)
target_sources(jensen_gap PRIVATE ${SRCS})

add_executable(ncf_uniform_test exec/ncf_uniform_test.cpp)
target_sources(ncf_uniform_test PRIVATE ${SRCS})

add_executable(ncf_near_col_test exec/ncf_near_col_test.cpp)
target_sources(ncf_near_col_test PRIVATE ${SRCS})

add_executable(split_uniform_test exec/split_uniform_test.cpp)
target_sources(split_uniform_test PRIVATE ${SRCS})

add_executable(primal_meet_lwe exec/primal_meet_lwe.cpp)
target_sources(primal_meet_lwe PRIVATE ${SRCS})

