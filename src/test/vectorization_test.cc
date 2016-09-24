#include<Eigen/Core>
using namespace Eigen;
void foo(Vector4f& u, Vector4f& v, Vector4f& w)
{
  EIGEN_ASM_COMMENT("begin");
  u = v + 3*w;
  EIGEN_ASM_COMMENT("end");
}
