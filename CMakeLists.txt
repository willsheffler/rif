cmake_minimum_required(VERSION 3.0)
project(rif)

if(NOT DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    message("CMAKE_LIBRARY_OUTPUT_DIRECTORY not defined, setting to .")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY .)
endif()

if(DEFINED BOOST_ROOT)
    message("cmake: using BOOST_ROOT " ${BOOST_ROOT})
    include_directories(${BOOST_ROOT}/include)
    link_directories(${BOOST_ROOT}/lib)
else()
    message("cmake: NOT using BOOST_ROOT")
endif()

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    message("cmake: using ccache ${CCACHE_PROGRAM}")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
else()
    message("cmake: not using ccache")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} $ENV{CXXFLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCXX11 -DCXX14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} $ENV{LDFLAGS}")

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math") # eigen xform perf worse!!
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

include_directories("external/include")
include_directories("src/rif")
include_directories("/usr/local/include")
include_directories("/opt/local/include")
link_directories("/usr/local/lib")
#link_directories("/opt/local/lib")
list(APPEND extra_libs boost_system)
list(APPEND extra_libs boost_iostreams)

function(print_config)
    message("cmake: CC:  " $ENV{CC})
    message("cmake: CXX: " $ENV{CXX})
    message("cmake: CXXFLAGS: " $ENV{CXXFLAGS})
    message("cmake: CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER})
    message("cmake: CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
    message("cmake: CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
    message("cmake: extra_libs: " ${extra_libs})
    message("cmake: CMAKE_OPTIONS: " $ENV{CMAKE_OPTIONS})
    message("cmake: BOOST_ROOT: " $ENV{BOOST_ROOT})
endfunction()

function(make_any_missing__init__py_and_copy_py_to_liboutdir)
    execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/tools/gen__init__.py src/rif
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} RESULT_VARIABLE STATUS)
    if(NOT "${STATUS}" STREQUAL "0")
        message(FATAL_ERROR "gen__init__.py failed: ${STATUS}")
    endif()
    execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/tools/cp_filtered.py .py
        ${CMAKE_SOURCE_DIR}/src/rif ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/rif
        RESULT_VARIABLE STATUS)
    if(NOT "${STATUS}" STREQUAL "0")
        message(FATAL_ERROR "cp_filtered.py failed: ${STATUS}")
    endif()
endfunction()

function(build_python_bindings)
    make_any_missing__init__py_and_copy_py_to_liboutdir()
    set(SRCGEN ${CMAKE_SOURCE_DIR}/tools/pybind_source_gen.py)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/rif.gen.cpp COMMAND python ${SRCGEN} ${CMAKE_SOURCE_DIR}/src/rif ${CMAKE_CURRENT_BINARY_DIR} DEPENDS src/rif/rif.gen.cpp.jinja)
    FILE(GLOB_RECURSE PY_CPP "./src/rif/[0-9a-zA-Z_]*.pybind.cpp")
    pybind11_add_module(rif_cpp ${CMAKE_CURRENT_BINARY_DIR}/rif.gen.cpp ${PY_CPP})
    target_link_libraries( rif_cpp PRIVATE ${extra_libs} )
endfunction()


function(build_gtests)
    add_subdirectory(external/googletest)
    include_directories("external/googletest/include")
    list(APPEND extra_libs gtest)

    # shitty "data" files to compile faster TODO: remove these
    list(APPEND LIB_CPP src/rif/rosetta/score/EtableParams_init.lib.cpp
                        src/rif/rosetta/atype_map.lib.cpp
                        src/rif/chem/ligand_factory.lib.cpp)
    set_source_files_properties(${LIB_CPP} PROPERTIES COMPILE_FLAGS -O0)
    add_library(rif_data ${LIB_CPP})

    file(GLOB_RECURSE TEST_CPP "src/rif/[0-9a-zA-Z_]*.gtest.cpp")
    add_executable(gtest_all src/rif/test/gtest_all.cpp ${TEST_CPP})
    target_link_libraries(gtest_all rif_data ${extra_libs})

    add_executable(gtest_wip src/rif/test/gtest_wip.cpp)
    target_link_libraries(gtest_wip rif_data ${extra_libs})
endfunction()

function(SUMMARY)
    if(CCACHE_PROGRAM)
        execute_process(COMMAND ccache -s
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} RESULT_VARIABLE STATUS)
    endif()
endfunction()









print_config()

add_subdirectory(external/pybind11)
build_python_bindings("")

build_gtests()

SUMMARY()


