#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_riflib_numeric_pigen(py::module & riflib);
void init_test_example(py::module & riflib);


PYBIND11_PLUGIN(riflib) {
    py::module riflib("riflib", R"pbdoc(
        riflib docs
        -----------------------

        .. currentmodule:: riflib

        .. autosummary::
           :toctree: _generate

    )pbdoc");


    init_riflib_numeric_pigen(riflib);
    init_test_example(riflib);


#ifdef VERSION_INFO
    riflib.attr("__version__") = py::str(VERSION_INFO);
#else
    riflib.attr("__version__") = py::str("dev");
#endif

    return riflib.ptr();

}
