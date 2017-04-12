#include "pyutil/pybind_simplearray.hpp"

#include <util/str.hpp>
#include "bcc_lattice.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace rif::numeric;
using namespace rif::util;

template <int N, class F, class I>
auto bcc_center_impl(BCC<N, F, I>& self, py::array_t<I> arg) {
  auto v = [&self](I idx) { return self[idx]; };
  auto r = py::vectorize(v)(arg);
  // py::dtype dt(str(N) + short_str<F>());
  return r;
}

template <int N, class F, class I>
auto bcc_index_impl(BCC<N, F, I>& self,
                    py::array_t<typename BCC<N, F, I>::Floats> arg) {
  using Floats = typename BCC<N, F, I>::Floats;
  auto v = [&self](Floats fs) { return self[fs]; };
  auto r = py::vectorize(v)(arg);
  // py::dtype dt(str(N) + short_str<I>());
  return r;
}

template <int N, class F, class I>
void bind_bcc_N_F_I(py::module& m) {
  using T = BCC<N, F, I>;
  std::string name = "BCC" + str(N) + short_str<F>() + short_str<I>();
  // std::cout << "!!!!!!!!!!!! BCC!!!!!!!!!!!! " << name << std::endl;
  bcc_center_impl<N, F, I>;  // instantiate template
  bcc_index_impl<N, F, I>;   // instantiate template
  py::class_<T>(m, name.c_str())
      .def(py::init<typename T::Indices, typename T::Floats,
                    typename T::Floats>(),
           "nc"_a, "lb"_a, "ub"_a)
      .def_property_readonly("dim", &T::dim)
      .def("__len__", &T::size)
      .def("_center_impl", &bcc_center_impl<N, F, I>)
      .def("_index_impl", &bcc_index_impl<N, F, I>);
}

template <int N>
void bind_bcc_N(py::module& m) {
  // bind_bcc_N_F_I<N, float, int32_t>(m);
  // bind_bcc_N_F_I<N, float, int64_t>(m);
  // bind_bcc_N_F_I<N, double, int32_t>(m);
  bind_bcc_N_F_I<N, double, int64_t>(m);
}

void RIFLIB_PYBIND_numeric_bcc_lattice(py::module& m) {
  bind_bcc_N<3>(m);
  bind_bcc_N<4>(m);
  bind_bcc_N<5>(m);
  bind_bcc_N<6>(m);
  bind_bcc_N<7>(m);
  bind_bcc_N<8>(m);
  bind_bcc_N<9>(m);
  bind_bcc_N<10>(m);
}