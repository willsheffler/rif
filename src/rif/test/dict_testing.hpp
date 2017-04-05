#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace rif {
namespace test {

namespace py = pybind11;

static std::string print_dict(py::object obj) {
  py::dict dict = obj.cast<py::dict>();
  std::ostringstream out;
  for (auto item : dict)
    out << "key=" << std::string(py::str(item.first)) << ", "
        << "value=" << std::string(py::str(item.second)) << std::endl;
  // std::cout << sizeof(py::object) << " " << sizeof(py::dict) << std::endl;
  return out.str();
}

static auto key_names_and_len(py::object obj) {
  py::dict dict = obj.cast<py::dict>();
  std::string s;
  int totlen = 0;
  for (auto item : dict) {
    s += py::str(item.first);
    py::array a = item.second.cast<py::array>();
    std::cout << std::string(py::str(item.first)) << ", " << a.size()
              << std::endl;
    totlen += a.size();
  }
  return std::make_pair(s, totlen);
}
}
}