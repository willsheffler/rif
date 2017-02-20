#include <math.h>

namespace scheme {
namespace numeric {

// assumes we are on the +Z face!!!!!  order Z,Xform,-Y,-Xform,Y,-Z
template <class Vec>
void cube_to_sphere(Vec &v) {
  typedef typename Vec::Scalar Float;
  double const xx = v.x() * v.x();
  double const yy = v.y() * v.y();
  double const zz = v.z() * v.z();
  v.x() *= sqrt(1.0 - yy * 0.5 - zz * 0.5 + yy * zz / 3.0);
  v.y() *= sqrt(1.0 - zz * 0.5 - xx * 0.5 + zz * xx / 3.0);
  v.z() *= sqrt(1.0 - xx * 0.5 - yy * 0.5 + xx * yy / 3.0);
  v.normalize();
}

template <class Vec>
void sphere_to_cube_facenum0(Vec &v) {
  typedef typename Vec::Scalar Float;
  Float const x = v.x();
  Float const y = v.y();
  // Float const z = v.z();
  Float &Xform(v.x());
  Float &Y(v.y());
  Float &Z(v.z());
  Float const INV_SQRT_2 = 0.70710676908493042;
  Float const a2 = x * x * 2.0;
  Float const b2 = y * y * 2.0;
  Float const inner = -a2 + b2 - 3.0;
  Float const innersqrt = -sqrt((inner * inner) - 12.0 * a2);
  Xform =
      (fabs(x) < 0.000001) ? 0.0 : sqrt(innersqrt + a2 - b2 + 3.0) * INV_SQRT_2;
  Y = (fabs(y) < 0.000001) ? 0.0 : sqrt(innersqrt - a2 + b2 + 3.0) * INV_SQRT_2;
  Z = 1.0;  // not really necessary for lookups...
  Xform = x < 0.0 ? -Xform : Xform;
  Y = y < 0.0 ? -Y : Y;
  // if(Xform > 1.0) Xform = 1.0; // shouldn't be necessary
  // if(Y > 1.0) Y = 1.0; // shouldn't be necessary
}

template <class Float>
void permute_cube_face_xyz(int face, Float &x, Float &y, Float &z) {
  Float X = x, Y = y, Z = z;
  switch (face) {
    case 0:
      x = X;
      y = Y;
      z = Z;
      break;
    case 1:
      x = -Y;
      y = Z;
      z = -X;
      break;
    case 2:
      x = -Z;
      y = -Y;
      z = -X;
      break;
    case 3:
      x = Y;
      y = -Z;
      z = -X;
      break;
    case 4:
      x = Z;
      y = Y;
      z = -X;
      break;
    case 5:
      x = Y;
      y = X;
      z = -Z;
      break;
    default:
      break;
  }
}

template <class Float>
void inverse_permute_cube_face_xyz(int face, Float &x, Float &y, Float &z) {
  Float X = x, Y = y, Z = z;
  switch (face) {
    case 0:
      x = X;
      y = Y;
      z = Z;
      break;
    case 1:
      x = -Z;
      y = -X;
      z = Y;
      break;
    case 2:
      x = -Z;
      y = -Y;
      z = -X;
      break;
    case 3:
      x = -Z;
      y = X;
      z = -Y;
      break;
    case 4:
      x = -Z;
      y = Y;
      z = X;
      break;
    case 5:
      x = Y;
      y = X;
      z = -Z;
      break;
    default:
      break;
  }
}
template <class Vec>
uint64_t get_cube_facenum(Vec const &p) {
  typedef typename Vec::Scalar Float;
  Float const ax = fabs(p.x()), ay = fabs(p.y()), az = fabs(p.z());
  uint64_t facenum = 0;
  if (ax >= ay && ax >= az)
    facenum = (p.x() > 0) ? 4 : 2;
  else if (ay >= az && ay >= ax)
    facenum = (p.y() > 0) ? 1 : 3;
  else if (az >= ax && az >= ay)
    facenum = (p.z() > 0) ? 0 : 5;
  return facenum;
}

// old usage from rosetta branch sheffler/scheme/master
// void NestQSph::set_position( uint64_t index, uint64_t resl ) {
// 	debug_assert_msg(index < size(resl),"index out of bounds");
// 	debug_assert_msg(resl  < 64,"resl >= 0 and < 64");
// 	uint64_t const nhier = ONE << (TWO*resl);
// 	uint64_t const baseindex = ((index >>TWO*resl)<<resl);
// 	uint64_t const hierindex = (index & (nhier-1));
// 	uint64_t const hi0 = utility::zorder::undilate<2>(hierindex     );
// 	uint64_t const hi1 = utility::zorder::undilate<2>(hierindex>>ONE);
// 	set_position(baseindex+hi0,hi1,resl);
//  }
// void NestQSph::set_position( uint64_t index1, uint64_t index2, uint64_t resl
// ) {
// 	uint64_t facenum = index1>>resl;
// 	// cout << "get_position face: " << facenum << endl;
// 	index1 = index1 & ((ONE<<resl)-ONE);
// 	debug_assert_msg( facenum < 6,"NestQSph must have 0 <= base_index < 6");
// 	uint64_t const nside = (ONE<<resl);
// 	double const x = 2.0*((double)index1+0.5)/(double)nside-1.0;
// 	double const y = 2.0*((double)index2+0.5)/(double)nside-1.0;
// 	Vec psphere(x,y,1.0);
// 	cube_to_sphere(psphere);
// 	debug_assert_msg( -1.0 <= psphere.x() && psphere.x() <= 1.0, "chi out of
// bounds" );
// 	debug_assert_msg( -1.0 <= psphere.y() && psphere.y() <= 1.0, "psi out of
// bounds" );
// 	debug_assert_msg( -1.0 <= psphere.z() && psphere.z() <= 1.0, "psi out of
// bounds" );
//     permute_cube_face_xyz(facenum,psphere.x(),psphere.y(),psphere.z());
// 	managed_xform_ptr_->t = psphere;
//  }
// uint64_t NestQSph::get_index( Xform const & position, uint64_t resl ) const {
// 	Vec p = position.t.normalized();
// 	uint64_t const facenum = get_cube_facenum(p);
// 	// cout << "get_index face " << facenum << endl;

// 	// cout << facenum << endl;
// 	inverse_permute_cube_face_xyz(facenum,p.x(),p.y(),p.z()); // move to
// face 0
// 	sphere_to_cube_facenum0(p);

// 	uint64_t const nside = (ONE<<resl);
// 	double const delta0 = p.x()/2.0 + 0.5;
// 	double const delta1 = p.y()/2.0 + 0.5;
// 	uint64_t const index0unbound = delta0 * (double)nside;
// 	uint64_t const index1unbound = delta1 * (double)nside;
// 	uint64_t const index0 = std::min(index0unbound,nside-ONE);
// 	uint64_t const index1 = std::min(index1unbound,nside-ONE);
// 	uint64_t const hi0 = utility::zorder::dilate<2>( index0 & (nside-ONE) );
// 	uint64_t const hi1 = utility::zorder::dilate<2>( index1 & (nside-ONE) );
// 	uint64_t const index = facenum<<(TWO*resl) | hi0 | (hi1<<ONE);
// 	return index;
//  }

}  // numeric
}  // scheme
