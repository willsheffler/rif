#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace pyutil {

template <class Ary, class Val>
auto sequence_to(pybind11::object obj, size_t n) {
  Ary a;
  for (size_t i = 0; i < n; ++i)
    a[i] = pybind11::cast<Val>(obj[pybind11::int_(i)]);
  return a;
}

template <class Container>
auto to_py_buffer_info(Container const& orig) {
  using Payload = typename Container::value_type;
  Payload* raw = new Payload[orig.size()];
  std::copy(orig.begin(), orig.end(), raw);
  return pybind11::buffer_info(raw, sizeof(Payload),
                               pybind11::format_descriptor<Payload>::format(),
                               1, {orig.size()}, {sizeof(Payload)});
}

template <class Container>
auto to_py_array(Container const& orig) {
  return pybind11::array(to_py_buffer_info(orig));
}
}
