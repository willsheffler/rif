#include <Eigen/Geometry>

namespace rif {

using V = Eigen::Vector3f;
using M = Eigen::Matrix3f;
using X = Eigen::Transform<float, 3, Eigen::Affine>;

X xform(M m, V v) {
  X x(m);
  x.translation() = v;
  return x;
}
X xform(V v) { return xform(M::Identity(), v); }
}
