#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

#pragma once

namespace pybind11 {
namespace detail {

template <typename Scalar, std::size_t N>
struct npy_format_descriptor<std::array<Scalar, N>> {
  using T = std::array<Scalar, N>;

  static PYBIND11_DESCR name() { return make_caster<T>::name(); }

  static pybind11::dtype dtype() {
    return reinterpret_borrow<pybind11::dtype>(dtype_ptr());
  }

  static std::string format() {
    static auto format_str =
        get_numpy_internals().get_type_info<T>(true)->format_str;
    return format_str;
  }

  static void register_dtype() {
    const std::type_info &tinfo = typeid(T);
    auto tindex = std::type_index(tinfo);

    auto &numpy_internals = get_numpy_internals();
    if (numpy_internals.get_type_info(tinfo, false)) {
      pybind11_fail("NumPy: dtype is already registered.)");
    }

    std::string scalar_fmt = format_descriptor<Scalar>::format();

    std::ostringstream oss2;
    oss2 << N << scalar_fmt;
    std::string format_str = oss2.str();

    auto dtype_ptr = pybind11::dtype(format_str).release().ptr();

    numpy_internals.registered_dtypes[tindex] = {dtype_ptr, format_str};
    get_internals().direct_conversions[tindex].push_back(direct_converter);
  }

 private:
  static PyObject *dtype_ptr() {
    static PyObject *ptr =
        get_numpy_internals().get_type_info<T>(true)->dtype_ptr;
    return ptr;
  }

  static bool direct_converter(PyObject *obj, void *&value) {
    auto &api = npy_api::get();
    if (!PyObject_TypeCheck(obj, api.PyVoidArrType_Type_)) return false;
    if (auto descr =
            reinterpret_steal<object>(api.PyArray_DescrFromScalar_(obj))) {
      if (api.PyArray_EquivTypes_(dtype_ptr(), descr.ptr())) {
        value = ((PyVoidScalarObject_Proxy *)obj)->obval;
        return true;
      }
    }
    return false;
  }
};
}  // end detail

template <class Scalar, std::size_t N>
struct format_descriptor<std::array<Scalar, N>> {
  using T = std::array<Scalar, N>;
  static std::string format() {
    return detail::npy_format_descriptor<
        typename std::remove_cv<T>::type>::format();
  }
};

}
