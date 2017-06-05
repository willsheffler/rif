#include "pyutil/pybind_numpy.hpp"

#include "rif/geom/Ray.hpp"
#include "rif/util/str.hpp"

namespace py = pybind11;
using namespace rif::geom;

namespace std {
template <class F>
struct is_pod<Ray<F>> : public std::integral_constant<bool, true> {};
}

template <class F>
void bind_ray(py::module &m) {
  using R = Ray<F>;

  PYBIND11_NUMPY_DTYPE(R, orig, dirn);

  std::string name1 = "rand_ray_gaussian" + rif::util::short_str<F>();
  m.def(name1.c_str(), py::vectorize(rand_ray_gaussian<F>));

  std::string name2 = "Ray" + rif::util::short_str<F>();
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!" << name2 << std::endl;
  py::class_<R> cls(m, name2.c_str());
  cls.def_readwrite("orig", &R::orig)
      .def_readwrite("dirn", &R::dirn)
      .def("__getitem__",
           [](R &self, int i) {
             if (0 == i) return self.orig;
             if (1 == i) return self.dirn;
             throw std::invalid_argument("Index out of range");
             return self.orig;
           })
      .def("__repr__", [](R const &r) {
        std::ostringstream oss;
        oss << "Ray(orig=[" << r.orig[0] << ", " << r.orig[1] << ", "
            << r.orig[2] << "], dirn=[" << r.dirn[0] << ", " << r.dirn[1]
            << ", " << r.dirn[2] << "])";
        return oss.str();
      });
  cls.attr("dtype") = py::dtype::of<R>();
}

void RIFLIB_PYBIND_geom_Ray(py::module &m) {
  bind_ray<float>(m);
  bind_ray<double>(m);
  m.attr("Ray") = m.attr("Rayf4");
  m.attr("rand_ray_gaussian") = m.attr("rand_ray_gaussianf4");
}
