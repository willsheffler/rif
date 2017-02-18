cmake_minimum_required(VERSION 3.3)
project(riflib)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -DCXX11 -DCXX14")
if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")    
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math") # eigen xform perf worse!!
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

include_directories("external/include")
include_directories("src")
include_directories("src/riflib")
include_directories("/usr/local/include")
include_directories("/opt/local/include")
link_directories("/usr/local/lib")
#link_directories("/opt/local/lib")

function(PRINT_CONFIG)
    message("CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
    message("CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
    message("EXTRA_LIBS: " ${EXTRA_LIBS})
endfunction()


function(build_python_bindings)
    execute_process(COMMAND python src/sourcegen.py WORKING_DIRECTORY ../..)
    FILE(GLOB_RECURSE PY_CPP "./src/[0-9a-zA-Z_]*.pybind.cpp")
    pybind11_add_module( riflib src/riflib.gen.cpp ${PY_CPP} )
endfunction()


function(BUILD_GTESTS)
    add_subdirectory(external/googletest)
    include_directories("external/googletest/include")
    list(APPEND EXTRA_LIBS gtest )
    list(APPEND EXTRA_LIBS boost_system)
    list(APPEND EXTRA_LIBS boost_iostreams)

    # shitty "data" files to compile faster TODO: remove these
    list(APPEND LIB_CPP src/riflib/rosetta/score/EtableParams_init.lib.cpp
                        src/riflib/rosetta/atype_map.lib.cpp
                        src/riflib/chemical/ligand_factory.lib.cpp)
    set_source_files_properties(${LIB_CPP} PROPERTIES COMPILE_FLAGS -O0)
    add_library(riflib_data ${LIB_CPP})


    file(GLOB_RECURSE TEST_CPP "src/[0-9a-zA-Z_]*.gtest.cpp")
    add_executable(gtest_all src/test/gtest_all.cpp ${TEST_CPP})
    target_link_libraries(gtest_all riflib_data ${EXTRA_LIBS})

    add_executable(gtest_wip src/test/gtest_wip.cpp)
    target_link_libraries(gtest_wip riflib_data ${EXTRA_LIBS})
endfunction()


PRINT_CONFIG()

add_subdirectory(external/pybind11)
build_python_bindings("")

BUILD_GTESTS()

