#include <pybind11/numpy.h>
#include "util/SimpleArray.hpp"

namespace pybind11 {
namespace detail {

// Casts an Eigen type to numpy array.  If given a base, the numpy array
// references the src data,
// otherwise it'll make a copy.  writeable lets you turn off the writeable flag
// for the array.
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

// Takes an lvalue ref to some Eigen type and a (python) base object, creating a
// numpy array that
// reference the Eigen object's data with `base` as the python-registered base
// class (if omitted,
// the base will be set to None, and lifetime management is up to the caller).
// The numpy array is
// non-writeable if the given type is const.
template <typename Type>
handle simplearray_ref_array(Type &src, handle parent = none()) {
  // none here is to get past array's should-we-copy detection, which currently
  // always
  // copies when there is no base.  Setting the base to None should be harmless.
  return simplearray_array_cast(src, parent, !std::is_const<Type>::value);
}

// Takes a pointer to some dense, plain Eigen type, builds a capsule around it,
// then returns a numpy
// array that references the encapsulated data with a python-side reference to
// the capsule to tie
// its destruction to that of any dependent python objects.  Const-ness is
// determined by whether or
// not the Type of the pointer given is const.
template <typename Type>
handle simplearray_encapsulate(Type *src) {
  capsule base(src, [](void *o) { delete static_cast<Type *>(o); });
  return simplearray_ref_array(*src, base);
}

/* @brief A pybind11 type_caster for SimpleArray
 */
template <int DIM, class Element, bool init0>
class type_caster<rif::util::SimpleArray<DIM, Element, init0>> {
 public:
  using Type = ::rif::util::SimpleArray<DIM, Element, init0>;

  bool load(handle src, bool) {
    auto buf = array_t<Element>::ensure(src);
    if (!buf) return false;

    auto dims = buf.ndim();
    if (dims != 1) return false;

    for (std::size_t s = 0; s < dims; ++s)
      if (buf.strides(s) % sizeof(Element) != 0) return false;

    if (!std::is_const<Element>::value && !buf.writeable()) return false;

    size_t len = buf.shape(0);
    if (len != DIM) return false;
    size_t stride = buf.strides(0) / sizeof(Element);

    for (size_t i = 0; i < DIM; ++i)
      value[i] = const_cast<Element *>(buf.data())[i * stride];

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

    // bool writeable = !std::is_const<Element>::value;

    // std::vector<size_t> shape(1, DIM);
    // std::vector<size_t> strides(1, sizeof(Element));

    // handle owner = handle();

    // Element *ptr = new Element[DIM];
    // array_t<Element> result(std::move(shape), std::move(strides), ptr,
    // owner);

    // if (!std::is_const<Element>::value) {
    //   array_proxy(result.ptr())->flags &= ~npy_api::NPY_ARRAY_WRITEABLE_;
    // }

    // return result.release();
  }
  /* This part is normally created by the PYBIND11_TYPE_CASTER macro, which
   * can't be used here due to the partial specialization
   */
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
}
}
