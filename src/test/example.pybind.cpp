#include <pybind11/pybind11.h>

int add(int i, int j) { return i + j; }

namespace py = pybind11;

void init_test_example(py::module & riflib){
    py::module m = riflib.def_submodule("example", "pybind11 example plugin");

    m.def("add", &add, "A function which adds two numbers");

}
