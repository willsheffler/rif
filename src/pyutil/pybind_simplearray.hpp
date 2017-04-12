#pragma once

#include <pybind11/numpy.h>
#include "util/SimpleArray.hpp"

#include <iostream>

// this is actually counter-productive
//    ambiguous with pybinds template
// namespace std {
// template <int N, class E, bool b>
// struct is_pod<rif::util::SimpleArray<N, E, b>>
// : public std::integral_constant<bool, std::is_pod<E>::value> {};
// }

namespace pybind11 {
namespace detail {

template <class Type>
handle simplearray_array_cast(Type const &src, handle base = handle(),
                              bool writeable = true) {
  constexpr size_t elem_size = sizeof(typename Type::Scalar);
  std::vector<size_t> shape, strides;
  shape.push_back(src.size());
  strides.push_back(elem_size);
  array a(std::move(shape), std::move(strides), (typename Type::Scalar *)&src,
          base);
  if (!writeable)
    array_proxy(a.ptr())->flags &= ~detail::npy_api::NPY_ARRAY_WRITEABLE_;
  return a.release();
}

template <typename Type>
handle simplearray_ref_array(Type &src, handle parent = none()) {
  return simplearray_array_cast(src, parent, !std::is_const<Type>::value);
}

template <typename Type>
handle simplearray_encapsulate(Type *src) {
  capsule base(src, [](void *o) { delete static_cast<Type *>(o); });
  return simplearray_ref_array(*src, base);
}

/* @brief A pybind11 type_caster for SimpleArray, taken from pybind11/eigen.h
 */
template <int DIM, class E, bool init0>
class type_caster<rif::util::SimpleArray<DIM, E, init0>> {
 public:
  using Type = ::rif::util::SimpleArray<DIM, E, init0>;

  bool load(handle src, bool) {
    auto buf = array_t<E>::ensure(src);
    if (!buf) return false;

    auto dims = buf.ndim();
    if (dims != 1) return false;

    for (std::size_t s = 0; s < dims; ++s)
      if (buf.strides(s) % sizeof(E) != 0) return false;

    if (!std::is_const<E>::value && !buf.writeable()) return false;

    size_t len = buf.shape(0);
    if (len != DIM) return false;
    size_t stride = buf.strides(0) / sizeof(E);

    for (size_t i = 0; i < DIM; ++i)
      value[i] = const_cast<E *>(buf.data())[i * stride];

    return true;
  }

 private:
  template <class CType>
  static handle cast_impl(const CType &src, return_value_policy policy,
                          handle parent) {
    switch (policy) {
      case return_value_policy::take_ownership:
      case return_value_policy::automatic:
        return simplearray_encapsulate(src);
      case return_value_policy::move:
      case return_value_policy::copy:
        return simplearray_array_cast(*src);
      case return_value_policy::reference:
      case return_value_policy::automatic_reference:
        return simplearray_ref_array(*src);
      case return_value_policy::reference_internal:
        return simplearray_ref_array(*src, parent);
      default:
        throw cast_error("unhandled return_value_policy: should not happen!");
    };
  }

 protected:
  Type value;

 public:
  static PYBIND11_DESCR name() { return type_descr(_<Type>()); }

  // Normal returned non-reference, non-const value:
  static handle cast(Type &&src, return_value_policy /* policy */,
                     handle parent) {
    return cast_impl(&src, return_value_policy::move, parent);
  }
  // If you return a non-reference const, we mark the numpy array readonly:
  static handle cast(const Type &&src, return_value_policy /* policy */,
                     handle parent) {
    return cast_impl(&src, return_value_policy::move, parent);
  }
  // lvalue reference return; default (automatic) becomes copy
  static handle cast(Type &src, return_value_policy policy, handle parent) {
    if (policy == return_value_policy::automatic ||
        policy == return_value_policy::automatic_reference)
      policy = return_value_policy::copy;
    return cast_impl(&src, policy, parent);
  }
  // const lvalue reference return; default (automatic) becomes copy
  static handle cast(const Type &src, return_value_policy policy,
                     handle parent) {
    if (policy == return_value_policy::automatic ||
        policy == return_value_policy::automatic_reference)
      policy = return_value_policy::copy;
    return cast(&src, policy, parent);
  }
  // non-const pointer return
  static handle cast(Type *src, return_value_policy policy, handle parent) {
    return cast_impl(src, policy, parent);
  }
  // const pointer return
  static handle cast(const Type *src, return_value_policy policy,
                     handle parent) {
    return cast_impl(src, policy, parent);
  }

  operator Type *() { return &value; }
  operator Type &() { return value; }
  template <typename _T>
  using cast_op_type = pybind11::detail::cast_op_type<_T>;
};

// todo: code duplication with Eigen
template <int DIM, class Scalar, bool init0>
struct npy_format_descriptor<rif::util::SimpleArray<DIM, Scalar, init0>> {
  using T = rif::util::SimpleArray<DIM, Scalar, init0>;
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
    auto &numpy_internals = get_numpy_internals();
    if (numpy_internals.get_type_info(tinfo, false))
      pybind11_fail(
          "NumPy: dtype is already registered: SimpleArray<something, "
          "something>");

    std::ostringstream oss;
    std::string scalar_fmt = format_descriptor<Scalar>::format();
    oss << DIM << scalar_fmt;
    std::string dtype_str = oss.str();
    list names, formats, offsets;
    names.append(PYBIND11_STR_TYPE("crd"));
    formats.append(dtype_str);
    offsets.append(pybind11::int_(0));
    auto dtype_ptr =
        pybind11::dtype(names, formats, offsets, sizeof(T)).release().ptr();

    std::ostringstream oss2;
    oss2 << "T{" << DIM << scalar_fmt << ":crd:}";
    std::string format_str = oss2.str();

    // std::cout << "npy_format_descriptor" << std::endl;
    // std::cout << "    size: " << sizeof(T) << std::endl;
    // std::cout << "    dtype_str: '" << dtype_str << "'" << std::endl;
    // std::cout << "    format_str: '" << format_str << "'" << std::endl;

    auto &api = npy_api::get();
    auto arr = array(buffer_info(nullptr, sizeof(T), format_str, 1));
    if (!api.PyArray_EquivTypes_(dtype_ptr, arr.dtype().ptr()))
      pybind11_fail("NumPy: invalid buffer descriptor!");

    auto tindex = std::type_index(tinfo);
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
template <int DIM, class Scalar, bool init0>
struct format_descriptor<rif::util::SimpleArray<DIM, Scalar, init0>> {
  static std::string format() {
    return detail::npy_format_descriptor<
        rif::util::SimpleArray<DIM, Scalar, init0>>::format();
  }
};
}
