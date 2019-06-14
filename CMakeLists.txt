cmake_minimum_required(VERSION 2.8.9)
project(evo_q)
set (CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Debug)

include_directories(
    include
  )
add_library( evo_q SHARED
             src/util.cpp src/phenotype.cpp src/parse.cpp src/gene.cpp src/organism.cpp src/population.cpp )

set(APP_EXE main_test)
add_executable( ${APP_EXE}
                src/main_test.cpp )
target_link_libraries( ${APP_EXE}
                       evo_q )

enable_testing()
add_test(NAME evo_q_test COMMAND ${APP_EXE})
