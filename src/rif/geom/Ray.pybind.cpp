#include "pyutil/pybind_numpy.hpp"

#include "rif/geom/Ray.hpp"
#include "rif/util/str.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace rif::geom;

namespace pybind11 {
namespace detail {
template <class F>
struct npy_format_descriptor<Ray<F>>
    : npy_format_descriptor<typename Ray<F>::M42> {};
}
template <class F>
struct format_descriptor<Ray<F>> : format_descriptor<typename Ray<F>::M42> {};
}

template <class F>
void bind_ray(py::module &m) {
  using R = Ray<F>;

  std::string name1 = "rand_ray_gaussian" + rif::util::short_str<F>();
  m.def(name1.c_str(), py::vectorize(rand_ray_gaussian<F>), "sd"_a);
  m.def(name1.c_str(),
        [](F sd, size_t size) {
          py::array_t<R> out(size);
          R *ptr = (R *)out.request().ptr;
          for (int i = 0; i < size; ++i) ptr[i] = rand_ray_gaussian<F>(sd);
          return out;
        },
        "sd"_a, "size"_a);

  std::string name2 = "Ray" + rif::util::short_str<F>();
  // std::cout << "!!!!!!!!!!!!!!!!!!!!!!" << name2 << std::endl;
  py::class_<R> cls(m, name2.c_str(), py::buffer_protocol());
  cls.def_property("orig", &R::getorig, &R::setorig)
      .def_property("dirn", &R::getdirn, &R::setdirn)
      // .def_readwrite("m42", &R::_m42)
      .def("__getitem__",
           [](R &self, int i) {
             if (0 == i) return self.getorig();
             if (1 == i) return self.getdirn();
             throw std::invalid_argument("Index out of range");
             return self.getorig();
           })
      .def("print_raw", [](R r) { std::cout << r._m42 << std::endl; })
      .def("reset_meta", [](R *r) { r->_m42.row(3) = rif::A2<F>(1, 0); })
      .def("__repr__",
           [](R r) {
             std::ostringstream oss;
             oss << "Ray(orig=[" << r.orig()[0] << ", " << r.orig()[1] << ", "
                 << r.orig()[2] << "], dirn=[" << r.dirn()[0] << ", "
                 << r.dirn()[1] << ", " << r.dirn()[2] << "])";
             return oss.str();
           })
      .def_buffer([](R &r) -> py::buffer_info {
        std::vector<size_t> shape = {4, 2};
        std::vector<size_t> stride = {sizeof(F) * 4, sizeof(F)};
        return py::buffer_info(r._m42.data(), sizeof(F),
                               py::format_descriptor<F>::format(), 2, shape,
                               stride);
      });

  // std::cout << "register M42" << std::endl;
  // py::detail::npy_format_descriptor<typename R::M42>::register_dtype();
  // std::cout << "register R" << std::endl;
  py::detail::npy_format_descriptor<R>::register_dtype();
  cls.attr("dtype") = py::dtype::of<R>();

  // using R2 = Ray<double>;
  // m.def("pyarray_Ray_test",
  // [](py::array_t<R> rays1, py::array_t<R> rays2) { return 8; });
}

void RIFLIB_PYBIND_geom_Ray(py::module &m) {
  bind_ray<float>(m);
  bind_ray<double>(m);
  m.attr("Ray") = m.attr("Rayf4");
  m.attr("rand_ray_gaussian") = m.attr("rand_ray_gaussianf4");
}
