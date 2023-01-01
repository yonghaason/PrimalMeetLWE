cmake_minimum_required(VERSION 3.19)
project(PrimHyb)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-march=native -O2 -std=c++2a)
add_executable(main main.cpp)
file(GLOB_RECURSE SRCS src/*.cpp)
target_sources(main PRIVATE ${SRCS})
