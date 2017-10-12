#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <rif/cluster/cookie_cutter.hpp>

namespace py = pybind11;
using namespace rif::cluster;
using namespace py::literals;

void RIFLIB_PYBIND_rif_cluster(py::module& m) {
  m.def("cookie_cutter_i1", &cookie_cutter<int8_t>, "data"_a, "thresh"_a,
        "metric"_a = "hamming");
  m.def("cookie_cutter_i2", &cookie_cutter<int16_t>, "data"_a, "thresh"_a,
        "metric"_a = "hamming");
  m.def("cookie_cutter_i4", &cookie_cutter<int32_t>, "data"_a, "thresh"_a,
        "metric"_a = "hamming");
  m.def("cookie_cutter_i8", &cookie_cutter<int64_t>, "data"_a, "thresh"_a,
        "metric"_a = "hamming");
  m.def("cookie_cutter_f4", &cookie_cutter<float>, "data"_a, "thresh"_a,
        "metric"_a = "L2");
  m.def("cookie_cutter_f8", &cookie_cutter<double>, "data"_a, "thresh"_a,
        "metric"_a = "L2");
  m.def("cookie_cutter_update_i1", &cookie_cutter_update<int8_t>, "data"_a,
        "thresh"_a, "centers"_a, "metric"_a = "hamming");
  m.def("cookie_cutter_update_i2", &cookie_cutter_update<int16_t>, "data"_a,
        "thresh"_a, "centers"_a, "metric"_a = "hamming");
  m.def("cookie_cutter_update_i4", &cookie_cutter_update<int32_t>, "data"_a,
        "thresh"_a, "centers"_a, "metric"_a = "hamming");
  m.def("cookie_cutter_update_i8", &cookie_cutter_update<int64_t>, "data"_a,
        "thresh"_a, "centers"_a, "metric"_a = "hamming");
  m.def("cookie_cutter_update_f4", &cookie_cutter_update<float>, "data"_a,
        "thresh"_a, "centers"_a, "metric"_a = "L2");
  m.def("cookie_cutter_update_f8", &cookie_cutter_update<double>, "data"_a,
        "thresh"_a, "centers"_a, "metric"_a = "L2");
}