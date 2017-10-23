#include <pybind11/pybind11.h>
#include "pyutil/pybind_numpy.hpp"
#include "rif/eigen_types.hpp"
#include "rif/hash/XformHash.hpp"
#include "rif/util/str.hpp"

namespace py = pybind11;
using namespace rif;
using namespace pybind11::literals;
using rif::util::str;

template <typename Cls>
void bind_xform_hash_common(Cls& cls) {
  using T = typename Cls::type;
  cls.def_property_readonly("phi_resl", &T::phi_resl);
  cls.def_property_readonly("cart_resl", &T::cart_resl);
  cls.def_property_readonly("ori_resl", &T::ori_resl);
  cls.def_property_readonly("cart_bound", &T::cart_bound);
  cls.def_property_readonly("ori_nside", &T::ori_nside);
  cls.def_property_readonly("grid", &T::grid);
}

void RIFLIB_PYBIND_rif_hash_xform_hash(py::module& m) {
  using XH = hash::XformHash_bt24_BCC6<X3f>;
  using XHA = hash::XformAngHash_bt24_BCC6<X3f>;
  using XHAA = hash::Xform2AngHash_bt24_BCC6<X3f>;
  using F = typename XH::Scalar;
  auto clsXH =
      py::class_<XH>(m, "XformHash_bt24_BCC6_X3f")
          .def(py::init<F, F, F>(), "cart_resl"_a, "ori_resl"_a,
               "cart_bound"_a = 512.0)
          .def("__repr__",
               [](XH const& xh) {
                 return "XformHash_bt24_BCC6_X3f(cart_resl=" +
                        str(xh.cart_resl_) + ", ori_resl=" + str(xh.ori_resl_) +
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

  auto clsXHA =
      py::class_<XHA, XH>(m, "XformAngHash_bt24_BCC6_X3f")
          .def(py::init<F, F, F, F>(), "phi_resl"_a, "cart_resl"_a,
               "ori_resl"_a, "cart_bound"_a = 512.0)
          .def("__repr__",
               [](XHA const& xha) {
                 return "XformAngHash_bt24_BCC6_X3f(phi_resl=" +
                        str(xha.phi_resl_) +
                        ", cart_resl=" + str(xha.cart_resl_) +
                        ", ori_resl=" + str(xha.ori_resl_) +
                        ", cart_bound=" + str(xha.cart_bound_) + ")";
               })
          .def("get_center", &XHA::get_center)
          // .def("get_center",
          // [](XHA const& self, py::array_t<int64_t> keys) {
          // auto f = [self](int64_t k) { return self.get_center(k); };
          // return py::vectorize(f)(keys);
          // })
          .def("get_key", (uint64_t(XHA::*)(X3f, float) const) & XHA::get_key)
          .def("get_key",
               (uint64_t(XHA::*)(std::pair<X3f, float>) const) & XHA::get_key)
          .def("get_key",
               [](XHA const& self, py::array_t<X3f> xforms,
                  py::array_t<float> phis) {
                 auto f = [self](X3f x, float phi) {
                   return self.get_key(x, phi);
                 };
                 return py::vectorize(f)(xforms, phis);
               })
      /**/;
  auto clsXHAA =
      py::class_<XHAA, XH>(m, "Xform2AngHash_bt24_BCC6_X3f")
          .def(py::init<F, F, F, F>(), "phi_resl"_a, "cart_resl"_a,
               "ori_resl"_a, "cart_bound"_a = 512.0)
          .def("__repr__",
               [](XHAA const& xha) {
                 return "Xform2AngHash_bt24_BCC6_X3f(phi_resl=" +
                        str(xha.phi_resl_) +
                        ", cart_resl=" + str(xha.cart_resl_) +
                        ", ori_resl=" + str(xha.ori_resl_) +
                        ", cart_bound=" + str(xha.cart_bound_) + ")";
               })
          .def("get_center", &XHAA::get_center)
          // .def("get_center",
          // [](XHAA const& self, py::array_t<int64_t> keys) {
          // auto f = [self](int64_t k) { return self.get_center(k); };
          // return py::vectorize(f)(keys);
          // })
          .def("get_key",
               (uint64_t(XHAA::*)(X3f, float, float) const) & XHAA::get_key)
          .def("get_key",
               (uint64_t(XHAA::*)(std::tuple<X3f, float, float>) const) &
                   XHAA::get_key)
          .def("get_key",
               [](XHAA const& self, py::array_t<X3f> xforms,
                  py::array_t<float> phis, py::array_t<float> phi2s) {
                 auto f = [self](X3f x, float phi, float phi2) {
                   return self.get_key(x, phi, phi2);
                 };
                 return py::vectorize(f)(xforms, phis, phi2s);
               })
      /**/;
  bind_xform_hash_common(clsXH);
  bind_xform_hash_common(clsXHA);
  bind_xform_hash_common(clsXHAA);
}

// m.def("vectorized_func",
//     [](py::array_t<int> x, py::array_t<float> y, my_custom_type *z) {
//         auto stateful_closure = [z](int x, float y) { return my_func(x, y,
//         z); };
//         return py::vectorize(stateful_closure)(x, y);
//     }
// );
