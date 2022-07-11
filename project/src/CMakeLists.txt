cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -ffast-math -std=c++11 -fopenmp")

find_package(OpenMP REQUIRED)
add_executable(groupsplit groupsplit.cpp)
