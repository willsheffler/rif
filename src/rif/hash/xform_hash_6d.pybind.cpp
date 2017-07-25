#include <pybind11/pybind11.h>
#include "pyutil/pybind_numpy.hpp"
#include "rif/eigen_types.hpp"
#include "rif/hash/XformHash.hpp"
#include "rif/util/str.hpp"

namespace py = pybind11;
using namespace rif;
using namespace pybind11::literals;
using rif::util::str;

void RIFLIB_PYBIND_rif_hash_xform_hash_6d(py::module& m) {
  using XH = hash::XformHash_bt24_BCC6<X3f>;
  using F = typename XH::Scalar;
  py::class_<XH>(m, "XformHash_bt24_BCC6_X3f")
      .def(py::init<F, F, F>(), "cart_resl"_a, "ang_resl"_a,
           "cart_bound"_a = 512.0)
      .def("__repr__",
           [](XH const& xh) {
             return "XformHash_bt24_BCC6_X3f(cart_resl=" + str(xh.cart_resl_) +
                    ", ang_resl=" + str(xh.ang_resl_) +
                    ", cart_bound=" + str(xh.cart_bound_) + ")";
           })
      .def("get_center", &XH::get_center)
      .def("get_center",
           [](XH const& self, py::array_t<int64_t> keys) {
             auto f = [self](int64_t k) { return self.get_center(k); };
             return py::vectorize(f)(keys);
           })
      .def("get_key", &XH::get_key)
      .def("get_key",
           [](XH const& self, py::array_t<X3f> xforms) {
             auto f = [self](X3f x) { return self.get_key(x); };
             return py::vectorize(f)(xforms);
           })
      /**/;
}

// m.def("vectorized_func",
//     [](py::array_t<int> x, py::array_t<float> y, my_custom_type *z) {
//         auto stateful_closure = [z](int x, float y) { return my_func(x, y,
//         z); };
//         return py::vectorize(stateful_closure)(x, y);
//     }
// );
