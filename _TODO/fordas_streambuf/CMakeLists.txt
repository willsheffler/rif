cmake_minimum_required(VERSION 3.1)
project(1705_ostream_example)

#find_package(FindPythonInterp REQUIRED)

#find_package(PythonLibs REQUIRED)
#include_directories(${PYTHON_INCLUDE_DIR})
execute_process(
  COMMAND python -c "from distutils import sysconfig; print(sysconfig.get_python_inc())"
  OUTPUT_VARIABLE PYTHON_INCLUDE_DIR
)
include_directories(${PYTHON_INCLUDE_DIR})

execute_process(
  COMMAND python -c "import pybind11; print(pybind11.get_include())"
  OUTPUT_VARIABLE PYBIND11_INCLUDE_DIR
)
include_directories(${PYBIND11_INCLUDE_DIR})

include_directories(.)

set (CMAKE_CXX_STANDARD 11)

add_library(ostream_example SHARED
  ostream_example.cpp
)

set_target_properties(ostream_example PROPERTIES PREFIX "")
