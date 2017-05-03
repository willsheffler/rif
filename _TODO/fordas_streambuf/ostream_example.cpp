#include <pybind11/pybind11.h>

#include <iostream>
#include <ostream>

#include <pystreambuf.h>

void testprint(std::ostream & o) {
  o << "testprint" << std::endl;
}

void testprint_noflush(std::ostream & o) {
  o << "testprint_noflush";
}

int testparse(std::istream & i) {
  int result;
  i >> result;
  return result;
}

namespace py = pybind11;

PYBIND11_PLUGIN(ostream_example) {
    py::module m("ostream_example");

    m.def("testprint", &testprint);
    m.def("testprint_noflush", &testprint_noflush);
    m.def("testparse", &testparse);

    return m.ptr();
}
