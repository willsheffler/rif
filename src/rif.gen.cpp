#include <pybind11/pybind11.h>

namespace py = pybind11;

void RIFLIB_PYBIND_sampling_orientations(py::module & m);
void RIFLIB_PYBIND_numeric_pigen(py::module & m);
void RIFLIB_PYBIND_test_numpy(py::module & m);
void RIFLIB_PYBIND_test_example_add(py::module & m);
void RIFLIB_PYBIND_test_example_sub(py::module & m);
void RIFLIB_PYBIND_test_example_mul(py::module & m);
void RIFLIB_PYBIND_test_example_dummy(py::module & m);


PYBIND11_PLUGIN(rif_cpp) {
    py::module rif("rif_cpp", R"pbdoc(
        rif docs
        -----------------------

        .. currentmodule:: rif

        .. autosummary::
           :toctree: _generate

    )pbdoc");

    py::module test__test_numpy = rif.def_submodule("test").def_submodule("test_numpy");
    py::module test__example = rif.def_submodule("test").def_submodule("example");
    py::module numeric__pigen = rif.def_submodule("numeric").def_submodule("pigen");
    py::module sampling__orientations = rif.def_submodule("sampling").def_submodule("orientations");

    RIFLIB_PYBIND_sampling_orientations(sampling__orientations);
    RIFLIB_PYBIND_numeric_pigen(numeric__pigen);
    RIFLIB_PYBIND_test_numpy(test__test_numpy);
    RIFLIB_PYBIND_test_example_add(test__example);
    RIFLIB_PYBIND_test_example_sub(test__example);
    RIFLIB_PYBIND_test_example_mul(test__example);
    RIFLIB_PYBIND_test_example_dummy(test__example);


#ifdef VERSION_INFO
    rif.attr("__version__") = py::str(VERSION_INFO);
#else
    rif.attr("__version__") = py::str("dev");
#endif

    return rif.ptr();

}