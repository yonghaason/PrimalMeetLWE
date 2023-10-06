cmake_minimum_required(VERSION 3.16)
project(PrimHyb)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-march=native -O2 -std=c++2a)
file(GLOB_RECURSE SRCS src/*.cpp)

include_directories(src)

add_executable(expectation_exchange exec/expectation_exchange.cpp)
target_sources(expectation_exchange PRIVATE ${SRCS})

add_executable(ncf_uniform_test exec/ncf_uniform_test.cpp)
target_sources(ncf_uniform_test PRIVATE ${SRCS})

add_executable(split_uniform_test exec/split_uniform_test.cpp)
target_sources(split_uniform_test PRIVATE ${SRCS})

add_executable(meet_lwe_matrix exec/meet_lwe_matrix.cpp)
target_sources(meet_lwe_matrix PRIVATE ${SRCS})