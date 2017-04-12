#include <pybind11/pybind11.h>

#define PYBIND11_FIELD_DESCRIPTOR_OFFSET(T, Ftype, Name, Offset)               \
  ::pybind11::detail::field_descriptor {                                       \
    #Name, Offset,                                                             \
        sizeof(Ftype),                                                         \
               alignof(                                                        \
                   Ftype),                                                     \
                   ::pybind11::format_descriptor < Ftype >                     \
                       ::format(), ::pybind11::detail::npy_format_descriptor < \
                                       Ftype > ::dtype()                       \
  }

#ifdef _MSC_VER
#define PYBIND11_MAP3_LIST_NEXT1(test, next) \
  PYBIND11_EVAL0(PYBIND11_MAP_NEXT0(test, PYBIND11_MAP_COMMA next, 0))
#else
#define PYBIND11_MAP3_LIST_NEXT1(test, next) \
  PYBIND11_MAP_NEXT0(test, PYBIND11_MAP_COMMA next, 0)
#endif
#define PYBIND11_MAP3_LIST_NEXT(test, next) \
  PYBIND11_MAP3_LIST_NEXT1(PYBIND11_MAP_GET_END test, next)
#define PYBIND11_MAP3_LIST0(f, t, x1, x2, x3, peek, ...)               \
  f(t, x1, x2, x3) PYBIND11_MAP3_LIST_NEXT(peek, PYBIND11_MAP3_LIST1)( \
      f, t, peek, __VA_ARGS__)
#define PYBIND11_MAP3_LIST1(f, t, x1, x2, x3, peek, ...)               \
  f(t, x1, x2, x3) PYBIND11_MAP3_LIST_NEXT(peek, PYBIND11_MAP3_LIST0)( \
      f, t, peek, __VA_ARGS__)
#define PYBIND11_MAP3_LIST(f, t, ...) \
  PYBIND11_EVAL(PYBIND11_MAP3_LIST1(f, t, __VA_ARGS__, (), 0))

#define PYBIND11_NUMPY_DTYPE_OFFSET(Type, ...)                     \
  ::pybind11::detail::npy_format_descriptor<Type>::register_dtype( \
      {PYBIND11_MAP3_LIST(PYBIND11_FIELD_DESCRIPTOR_OFFSET, Type,  \
                          __VA_ARGS__)})
