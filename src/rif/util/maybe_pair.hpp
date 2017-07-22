#pragma once

namespace rif {
namespace util {

template <class A>
A &get_first_if_pair(A &p) {
  return p;
}
template <class A>
A const &get_first_if_pair(A const &p) {
  return p;
}

template <class A, class B>
A &get_first_if_pair(std::pair<A, B> &p) {
  return p.first;
}

template <class A, class B>
A const &get_first_if_pair(std::pair<A, B> const &p) {
  return p.first;
}

template <class A>
A &get_second_if_pair(A &p) {
  return p;
}
template <class A>
A const &get_second_if_pair(A const &p) {
  return p;
}

template <class A, class B>
B &get_second_if_pair(std::pair<A, B> &p) {
  return p.second;
}
template <class A, class B>
B const &get_second_if_pair(std::pair<A, B> const &p) {
  return p.second;
}
}
}