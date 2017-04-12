#pragma once

#include <Eigen/Geometry>
#include <iostream>
#include <util/str.hpp>
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"

namespace pybind11 {
namespace detail {

using rif::util::str;
using rif::util::cpp_repr;

template <class Scalar, int NROW, int NCOL, int OPTS>
struct npy_format_descriptor<Eigen::Matrix<Scalar, NROW, NCOL, OPTS>> {
  using T = Eigen::Matrix<Scalar, NROW, NCOL, OPTS>;
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
    if (numpy_internals.get_type_info(tinfo, false))
      pybind11_fail("NumPy: dtype is already registered:\n    Eigen::Matrix<" +
                    cpp_repr<Scalar>() + ", " + str(NROW) + ", " + str(NCOL) +
                    ", " + str(OPTS) + ">\n    hash_code: " +
                    str(tindex.hash_code()));

    std::ostringstream oss;
    std::string scalar_fmt = format_descriptor<Scalar>::format();
    if (NCOL == 1) oss << NROW << scalar_fmt;
    if (NCOL >= 2) oss << "(" << NROW << "," << NCOL << ")" << scalar_fmt;
    std::string dtype_str = oss.str();
    list names, formats, offsets;
    names.append(PYBIND11_STR_TYPE("crd"));
    formats.append(dtype_str);
    offsets.append(pybind11::int_(0));
    auto dtype_ptr =
        pybind11::dtype(names, formats, offsets, sizeof(T)).release().ptr();

    std::ostringstream oss2;
    if (NCOL == 1) oss2 << "T{" << NROW << scalar_fmt << ":crd:}";
    if (NCOL >= 2)
      oss2 << "T{(" << NROW << "," << NCOL << ")" << scalar_fmt << ":crd:}";
    std::string format_str = oss2.str();

    // std::cout << "npy_format_descriptor" << std::endl;
    // std::cout << "    size: " << sizeof(T) << std::endl;
    // std::cout << "    dtype_str: '" << dtype_str << "'" << std::endl;
    // std::cout << "    format_str: '" << format_str << "'" << std::endl;

    auto &api = npy_api::get();
    auto arr = array(buffer_info(nullptr, sizeof(T), format_str, 1));
    if (!api.PyArray_EquivTypes_(dtype_ptr, arr.dtype().ptr()))
      pybind11_fail("NumPy: invalid buffer descriptor!");

    numpy_internals.registered_dtypes[tindex] = {dtype_ptr, format_str};
    get_internals().direct_conversions[tindex].push_back(direct_converter);
  }

 private:
  static PyObject *dtype_ptr() {
    // std::cout << "eigen dtype_ptr Matrix(" << sizeof(Scalar) * 8 << "bit "
    // << NROW << " " << NCOL << ")" << std::endl;
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
template <class Scalar, int NROW, int NCOL, int OPTS>
struct format_descriptor<Eigen::Matrix<Scalar, NROW, NCOL, OPTS>> {
  using T = Eigen::Matrix<Scalar, NROW, NCOL, OPTS>;
  static std::string format() {
    return detail::npy_format_descriptor<
        typename std::remove_cv<T>::type>::format();
  }
};
}
