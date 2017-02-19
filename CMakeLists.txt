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
include_directories("/usr/local/include")
include_directories("/opt/local/include")
link_directories("/usr/local/lib")
#link_directories("/opt/local/lib")
list(APPEND EXTRA_LIBS boost_system)
list(APPEND EXTRA_LIBS boost_iostreams)

function(PRINT_CONFIG)
    message("CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
    message("CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
    message("EXTRA_LIBS: " ${EXTRA_LIBS})
endfunction()


function(build_python_bindings)
    execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/devtools/sourcegen.py ${CMAKE_SOURCE_DIR}/src)
    FILE(GLOB_RECURSE PY_CPP "./src/[0-9a-zA-Z_]*.pybind.cpp")
    pybind11_add_module( riflib_cpp src/riflib.gen.cpp ${PY_CPP} )
    target_link_libraries( riflib_cpp PRIVATE ${EXTRA_LIBS} )
    execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/devtools/cp_filtered.py .py
        ${CMAKE_SOURCE_DIR}/src ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/riflib)
endfunction()


function(BUILD_GTESTS)
    add_subdirectory(external/googletest)
    include_directories("external/googletest/include")
    list(APPEND EXTRA_LIBS gtest )

    # shitty "data" files to compile faster TODO: remove these
    list(APPEND LIB_CPP src/rosetta/score/EtableParams_init.lib.cpp
                        src/rosetta/atype_map.lib.cpp
                        src/chemical/ligand_factory.lib.cpp)
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

