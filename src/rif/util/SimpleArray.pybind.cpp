// this was mostly copied from pybind11/eigen.h
//
#include "SimpleArray.hpp"
#include "util/pybind_simplearray.hpp"

template <class T>
T echo(T const &t) {
  return t;
}

using namespace rif::util;
namespace py = pybind11;

void RIFLIB_PYBIND_util_SimpleArray(py::module &m) {
  m.def("echo_SimpleArray_1_int", &echo<SimpleArray<1, int>>);
  m.def("echo_SimpleArray_3_int", &echo<SimpleArray<3, int>>);
}
