#include <Eigen/Dense>
#include <numeric/types.hpp>

namespace rif {
namespace numeric {

struct Ray {
  using V = Eigen::Vector3f;
  V origin;
  V direction;
  Ray() : origin(0, 0, 0), direction(1, 0, 0) {}
  Ray(V o, V d) : origin(o), direction(d.normalized()) {}
};

Ray operator*(X const& x, Ray const& r) {
  return Ray(x * r.origin, x.rotation() * r.direction);
}
}
}
