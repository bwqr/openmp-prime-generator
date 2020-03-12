cmake_minimum_required(VERSION 3.15)
project(prime-generator)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11")

add_executable(prime-generator main.cpp)