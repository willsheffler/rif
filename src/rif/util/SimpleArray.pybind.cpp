// this was mostly copied from pybind11/eigen.h

#include "SimpleArray.hpp"
#include "pyutil/pybind_simplearray.hpp"
#include "str.hpp"

struct SATest {
  rif::util::SimpleArray<2, int> m;
  rif::util::SimpleArray<2, int> &member() { return m; }
};

namespace std {
template <>
struct is_pod<SATest> : public std::integral_constant<bool, true> {};
}

template <class T>
T echo(T const &t) {
  return t;
}

using namespace rif::util;
namespace py = pybind11;

template <int N, class E>
void register_dtype_SimpleArray_N_E(py::module &m) {
  py::detail::npy_format_descriptor<SimpleArray<N, E>>::register_dtype();
  m.attr(("sa" + str(N) + short_str<E>() + "_t").c_str()) =
      py::dtype::of<SimpleArray<N, E>>();
}

template <class E>
void register_dtype_SimpleArrays_1_to_10(py::module &m) {
  register_dtype_SimpleArray_N_E<1, E>(m);
  register_dtype_SimpleArray_N_E<2, E>(m);
  register_dtype_SimpleArray_N_E<3, E>(m);
  register_dtype_SimpleArray_N_E<4, E>(m);
  register_dtype_SimpleArray_N_E<5, E>(m);
  register_dtype_SimpleArray_N_E<6, E>(m);
  register_dtype_SimpleArray_N_E<7, E>(m);
  register_dtype_SimpleArray_N_E<8, E>(m);
  register_dtype_SimpleArray_N_E<9, E>(m);
  register_dtype_SimpleArray_N_E<10, E>(m);
}

void RIFLIB_PYBIND_util_SimpleArray(py::module &m) {
  register_dtype_SimpleArrays_1_to_10<float>(m);
  register_dtype_SimpleArrays_1_to_10<uint32_t>(m);
  register_dtype_SimpleArrays_1_to_10<int32_t>(m);
  register_dtype_SimpleArrays_1_to_10<double>(m);
  register_dtype_SimpleArrays_1_to_10<uint64_t>(m);
  register_dtype_SimpleArrays_1_to_10<int64_t>(m);
  PYBIND11_NUMPY_DTYPE(SATest, m);
  m.def("echo_SimpleArray_1_int", &echo<SimpleArray<1, int>>);
  m.def("echo_SimpleArray_3_int", &echo<SimpleArray<3, int>>);
  py::class_<SATest>(m, "SATest")
      .def(py::init<>())
      .def("member", &SATest::member,
           py::return_value_policy::reference_internal);
}
