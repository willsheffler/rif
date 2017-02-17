#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_riflib_numeric_pigen(py::module & riflib);
void init_test_example(py::module & riflib);
void init_test_numpy(py::module & riflib_test);


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

    {
        auto m = riflib.def_submodule("test");
        init_test_numpy(m);
    }


#ifdef VERSION_INFO
    riflib.attr("__version__") = py::str(VERSION_INFO);
#else
    riflib.attr("__version__") = py::str("dev");
#endif

    return riflib.ptr();

}
