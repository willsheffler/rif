// this was mostly copied from pybind11/eigen.h
//
#include "SimpleArray.hpp"
#include "pyutil/pybind_simplearray.hpp"

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
void register_dtype_SimpleArray_N_E() {
  py::detail::npy_format_descriptor<
      rif::util::SimpleArray<N, E>>::register_dtype();
}

template <class E>
void register_dtype_SimpleArrays_1_to_10() {
  register_dtype_SimpleArray_N_E<1, E>();
  register_dtype_SimpleArray_N_E<2, E>();
  register_dtype_SimpleArray_N_E<3, E>();
  register_dtype_SimpleArray_N_E<4, E>();
  register_dtype_SimpleArray_N_E<5, E>();
  register_dtype_SimpleArray_N_E<6, E>();
  register_dtype_SimpleArray_N_E<7, E>();
  register_dtype_SimpleArray_N_E<8, E>();
  register_dtype_SimpleArray_N_E<9, E>();
  register_dtype_SimpleArray_N_E<10, E>();
}

void RIFLIB_PYBIND_util_SimpleArray(py::module &m) {
  register_dtype_SimpleArrays_1_to_10<float>();
  register_dtype_SimpleArrays_1_to_10<uint32_t>();
  register_dtype_SimpleArrays_1_to_10<int32_t>();
  register_dtype_SimpleArrays_1_to_10<double>();
  register_dtype_SimpleArrays_1_to_10<uint64_t>();
  register_dtype_SimpleArrays_1_to_10<int64_t>();
  PYBIND11_NUMPY_DTYPE(SATest, m);
  m.def("echo_SimpleArray_1_int", &echo<SimpleArray<1, int>>);
  m.def("echo_SimpleArray_3_int", &echo<SimpleArray<3, int>>);
  py::class_<SATest>(m, "SATest")
      .def(py::init<>())
      .def("member", &SATest::member,
           py::return_value_policy::reference_internal);
}
