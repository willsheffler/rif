#include <pybind11/pybind11.h>

int add(int i, int j) { return i + j; }
int sub(int i, int j) { return i - j; }
int mul(int i, int j) { return i * j; }

namespace py = pybind11;

void RIFLIB_PYBIND_test_example_add(py::module &m) { m.def("add", &add); }
void RIFLIB_PYBIND_test_example_sub(py::module &m) { m.def("sub", &sub); }
void RIFLIB_PYBIND_test_example_mul(py::module &m) { m.def("mul", &mul); }
void RIFLIB_PYBIND_test_example_dummy(py::module &m) {}
