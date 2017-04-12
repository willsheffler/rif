#include "pyutil/pybind_simplearray.hpp"

#include <util/str.hpp>
#include "bcc_lattice.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace rif::numeric;
using namespace rif::util;

template <int N, class F, class I>
void bind_bcc_N_F_I(py::module& m) {
  using T = BCC<N, F, I>;
  std::string name = "BCC" + str(N) + short_str<F>() + short_str<I>();
  // std::cout << "!!!!!!!!!!!! BCC!!!!!!!!!!!! " << name << std::endl;
  py::class_<T>(m, name.c_str())
      .def(py::init<typename T::Indices, typename T::Floats,
                    typename T::Floats>(),
           "ncell"_a, "lb"_a, "ub"_a)
      .def("__len__", &T::size)
      .def("__getitem__", [](T& self, py::array_t<I> arg) {
        auto v = [&self](I arg) { return self[arg]; };
        return py::vectorize(v)(arg);
      });
  ;
}

template <int N>
void bind_bcc_N(py::module& m) {
  bind_bcc_N_F_I<N, float, uint32_t>(m);
  bind_bcc_N_F_I<N, float, uint64_t>(m);
  bind_bcc_N_F_I<N, double, uint32_t>(m);
  bind_bcc_N_F_I<N, double, uint64_t>(m);
}

void RIFLIB_PYBIND_numeric_bcc_lattice(py::module& m) {
  // PYBIND_NUMPY_DTYPE()

  bind_bcc_N<3>(m);
  bind_bcc_N<4>(m);
  bind_bcc_N<5>(m);
  bind_bcc_N<6>(m);
  bind_bcc_N<7>(m);
  bind_bcc_N<8>(m);
  bind_bcc_N<9>(m);
  bind_bcc_N<10>(m);
}