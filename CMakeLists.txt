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



FILE(GLOB L1 "./[0-9a-zA-Z_]*.py.cpp")
FILE(GLOB L2 "./*/[0-9a-zA-Z_]*.py.cpp")
FILE(GLOB L3 "./*/*/[0-9a-zA-Z_]*.py.cpp")
FILE(GLOB L4 "./*/*/*/[0-9a-zA-Z_]*.py.cpp")
FILE(GLOB L5 "./*/*/*/*/[0-9a-zA-Z_]*.py.cpp")

pybind11_add_module( riflib ${L1} ${L2} ${L3} ${L4} ${L5} )
# pybind11_add_module( riflib src/test.py.cpp )




########################### tests #########################################################
# shitty "data" files to compile faster
# TODO remove these
list(APPEND LIB_CC src/scheme/rosetta/score/EtableParams_init.lib.cc
        src/scheme/rosetta/atype_map.lib.cc src/scheme/chemical/ligand_factory.lib.cc)
SET_SOURCE_FILES_PROPERTIES(${LIB_CC} PROPERTIES COMPILE_FLAGS -O0)
add_library(scheme ${LIB_CC})

list(APPEND EXTRA_LIBS boost_system)

include_directories(".")

FILE(GLOB L1 "src/[0-9a-zA-Z_]*.gtest.cc")
FILE(GLOB L2 "src/*/[0-9a-zA-Z_]*.gtest.cc")
FILE(GLOB L3 "src/*/*/[0-9a-zA-Z_]*.gtest.cc")
FILE(GLOB L4 "src/*/*/*/[0-9a-zA-Z_]*.gtest.cc")
FILE(GLOB L5 "src/*/*/*/*/[0-9a-zA-Z_]*.gtest.cc")

add_executable(gtest_all src/test/gtest_all.cc ${L1} ${L2} ${L3} ${L4} ${L5})
target_link_libraries(gtest_all scheme ${EXTRA_LIBS})

add_executable(gtest_wip src/test/gtest_wip.cc)
target_link_libraries(gtest_wip scheme ${EXTRA_LIBS})

