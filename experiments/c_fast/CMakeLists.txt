cmake_minimum_required(VERSION 3.16)
project(PrimHyb)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-march=native -O2 -std=c++2a)
file(GLOB_RECURSE SRCS src/*.cpp)

include_directories(src)

add_executable(jensen_gap exec/jensen_gap.cpp)
target_sources(jensen_gap PRIVATE ${SRCS})

add_executable(ncf_model_test exec/ncf_model_test.cpp)
target_sources(ncf_model_test PRIVATE ${SRCS})

add_executable(pair_model_test exec/pair_model_test.cpp)
target_sources(pair_model_test PRIVATE ${SRCS})


