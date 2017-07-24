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
auto new_pybuffer(Container const& orig) {
  using Payload = typename Container::value_type;
  // use malloc
  Payload* raw = (Payload*)malloc(orig.size() * sizeof(Payload));
  std::copy(orig.begin(), orig.end(), raw);
  return pybind11::buffer_info(raw, sizeof(Payload),
                               pybind11::format_descriptor<Payload>::format(),
                               1, {orig.size()}, {sizeof(Payload)});
}

template <class Container>
auto new_pyarray(Container const& orig) {
  return pybind11::array(new_pybuffer(orig));
}

template <class Container>
auto to_pybuffer(Container const& orig) {
  using Payload = typename Container::value_type;
  return pybind11::buffer_info((Payload*)(&orig[0]), sizeof(Payload),
                               pybind11::format_descriptor<Payload>::format(),
                               1, {orig.size()}, {sizeof(Payload)});
}

template <class Container>
auto to_pyarray(Container const& orig) {
  // this copies
  return pybind11::array(to_pybuffer(orig));
}
}
