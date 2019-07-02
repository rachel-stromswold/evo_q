cmake_minimum_required(VERSION 2.8.9)
project(evo_q)
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

# check if Doxygen is installed
if(CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
      # set input and output files
      set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile)
      set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

      # request to configure the file
      configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
      message("Doxygen build started")

      # note the option ALL which allows to build the docs together with the application
      add_custom_target( doc_doxygen ALL
	  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
	  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	  COMMENT "Generating API documentation with Doxygen"
	  VERBATIM )
  else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
  endif (DOXYGEN_FOUND)
endif()

#set build preprocessor macros
add_compile_definitions(USE_EXCEPTIONS=1)

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