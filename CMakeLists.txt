cmake_minimum_required(VERSION 2.8.9)
project(evo_q)
set (CMAKE_CXX_STANDARD 11)

#include pybind11 for python
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG        v2.2.3 )

include_directories(
    include )

FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

set(CMAKE_BUILD_TYPE Debug)

pybind11_add_module(evo_p src/util.cpp src/phenotype.cpp src/parse.cpp src/gene.cpp src/organism.cpp src/population.cpp src/convergence.cpp src/python.cpp)

add_library( evo_q SHARED
             src/util.cpp src/phenotype.cpp src/parse.cpp src/gene.cpp src/organism.cpp src/population.cpp src/convergence.cpp)

set(APP_EXE main_test)
add_executable( ${APP_EXE}
                src/main_test.cpp )
target_link_libraries( ${APP_EXE}
                       evo_q )

enable_testing()
add_test(NAME evo_q_test COMMAND ${APP_EXE})
