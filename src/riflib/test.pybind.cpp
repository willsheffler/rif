#include <pybind11/pybind11.h>

int add(int i, int j) {
    return i + j;
}

int subtract(int i, int j) {
    return i - j;
}

int mult(int i, int j) {
    return i * j;
}

namespace py = pybind11;


PYBIND11_PLUGIN(riflib) {
    py::module m("riflib", R"pbdoc(
        riflib docs
        -----------------------

        .. currentmodule:: riflib

        .. autosummary::
           :toctree: _generate

           add
           subtract
           mult
    )pbdoc");

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", &subtract, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("mult", &mult, R"pbdoc(
        Multiply two numbers

        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = py::str(VERSION_INFO);
#else
    m.attr("__version__") = py::str("dev");
#endif

    return m.ptr();

}
