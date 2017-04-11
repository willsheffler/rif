#include "util/pybind_simplearray.hpp"

#include "bcc_lattice.hpp"

#include <boost/lexical_cast.hpp>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace rif::numeric;

template <class T>
std::string str() {
  return "?";
}
template <>
std::string str<float>() {
  return "f4";
}
template <>
std::string str<double>() {
  return "f8";
}
template <>
std::string str<int32_t>() {
  return "i4";
}
template <>
std::string str<int64_t>() {
  return "i8";
}
template <>
std::string str<uint32_t>() {
  return "u4";
}
template <>
std::string str<uint64_t>() {
  return "u8";
}

template <int N, class F, class I>
void bind_bcc_N_F_I(py::module& m) {
  using T = BCC<N, F, I>;
  std::string name =
      "BCC" + boost::lexical_cast<std::string>(N) + str<F>() + str<I>();
  // std::cout << "!!!!!!!!!!!! BCC!!!!!!!!!!!! " << name << std::endl;
  py::class_<T>(m, name.c_str())
      .def(py::init<typename T::Indices, typename T::Floats,
                    typename T::Floats>(),
           "ncell"_a, "lb"_a, "ub"_a);
}

template <int N>
void bind_bcc_N(py::module& m) {
  bind_bcc_N_F_I<N, float, uint32_t>(m);
  bind_bcc_N_F_I<N, float, uint64_t>(m);
  bind_bcc_N_F_I<N, double, uint32_t>(m);
  bind_bcc_N_F_I<N, double, uint64_t>(m);
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