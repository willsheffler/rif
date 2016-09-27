cmake_minimum_required(VERSION 3.3)
project(riflib)

message("CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -DCXX11 -DCXX14")


#add_subdirectory(external/gmock)
#
#include_directories("external/gmock/gtest/src")
#include_directories("external/gmock/gtest/include")
#include_directories("external/gmock/include")
#list( APPEND EXTRA_LIBS gmock )

add_subdirectory(external/googletest)
include_directories("external/googletest/include")
list( APPEND EXTRA_LIBS gtest )

include_directories("external/include")

include_directories("src")


add_subdirectory(external/pybind11)



FILE(GLOB_RECURSE PYTHON_CPP "./src/[0-9a-zA-Z_]*.pybind.cpp")
pybind11_add_module( riflib ${PYTHON_CPP} )
# pybind11_add_module( riflib src/test.py.cpp )




########################### tests #########################################################
# shitty "data" files to compile faster
# TODO remove these
list(APPEND LIB_CPP src/riflib/rosetta/score/EtableParams_init.lib.cpp
        src/riflib/rosetta/atype_map.lib.cpp src/riflib/chemical/ligand_factory.lib.cpp)
SET_SOURCE_FILES_PROPERTIES(${LIB_CPP} PROPERTIES COMPILE_FLAGS -O0)
add_library(scheme ${LIB_CPP})

list(APPEND EXTRA_LIBS boost_system)

include_directories(".")

FILE(GLOB_RECURSE TEST_CPP "src/[0-9a-zA-Z_]*.gtest.cpp")

add_executable(gtest_all src/test/gtest_all.cpp ${TEST_CPP})
target_link_libraries(gtest_all scheme ${EXTRA_LIBS})

add_executable(gtest_wip src/test/gtest_wip.cpp)
target_link_libraries(gtest_wip scheme ${EXTRA_LIBS})

