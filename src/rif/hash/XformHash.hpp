#pragma once

#include "nest/pmap/TetracontoctachoronMap.hpp"
#include "numeric/lattice.hpp"
#include "numeric/util.hpp"
#include "util/SimpleArray.hpp"
#include "util/dilated_int.hpp"

#include <boost/utility/binary.hpp>

namespace rif {
namespace hash {

template <class A>
void clamp01(A &a) {
  for (int i = 0; i < A::N; ++i) {
    a[i] = fmin(1.0, fmax(0.0, a[i]));
  }
}

template <class _Xform>
struct XformHash_bt24_BCC6 {
  using Xform = _Xform;
  typedef uint64_t Key;
  typedef typename Xform::Scalar Float;
  typedef typename Xform::Scalar Scalar;
  typedef rif::nest::pmap::TetracontoctachoronMap<> OriMap;
  typedef rif::numeric::BCC<6, Float, uint64_t> Grid;
  typedef rif::util::SimpleArray<3, Float> F3;
  typedef rif::util::SimpleArray<3, uint64_t> I3;
  typedef rif::util::SimpleArray<6, Float> F6;
  typedef rif::util::SimpleArray<6, uint64_t> I6;

  Float grid_size_;
  Float grid_spacing_;
  OriMap ori_map_;
  Grid grid_;
  Float cart_resl_, ang_resl_, cart_bound_;

  static std::string name() { return "XformHash_bt24_BCC6"; }

  XformHash_bt24_BCC6() {}

  XformHash_bt24_BCC6(Float cart_resl, Float ang_resl,
                      Float cart_bound = 512.0) {
    init(cart_resl, ang_resl, cart_bound);
  }

  int get_ori_nside() {
    static float const covrad[64] = {
        49.66580, 25.99805, 17.48845, 13.15078, 10.48384, 8.76800, 7.48210,
        6.56491,  5.84498,  5.27430,  4.78793,  4.35932,  4.04326, 3.76735,
        3.51456,  3.29493,  3.09656,  2.92407,  2.75865,  2.62890, 2.51173,
        2.39665,  2.28840,  2.19235,  2.09949,  2.01564,  1.94154, 1.87351,
        1.80926,  1.75516,  1.69866,  1.64672,  1.59025,  1.54589, 1.50077,
        1.46216,  1.41758,  1.38146,  1.35363,  1.31630,  1.28212, 1.24864,
        1.21919,  1.20169,  1.17003,  1.14951,  1.11853,  1.09436, 1.07381,
        1.05223,  1.02896,  1.00747,  0.99457,  0.97719,  0.95703, 0.93588,
        0.92061,  0.90475,  0.89253,  0.87480,  0.86141,  0.84846, 0.83677,
        0.82164};
    int ori_nside = 1;
    while (covrad[ori_nside - 1] * 1.45 > ang_resl_ && ori_nside < 62)
      ++ori_nside;  // TODO: HACK multiplier!
    return ori_nside;
  }

  void init(Float cart_resl, Float ang_resl, Float cart_bound = 512.0) {
    cart_resl_ = cart_resl;
    ang_resl_ = ang_resl;
    cart_bound_ = cart_bound;
    init2(cart_resl, get_ori_nside(), cart_bound);
  }
  void init2(Float cart_resl, int ori_nside, Float cart_bound) {
    cart_resl /= sqrt(3.0) / 2.0;  // TODO: HACK multiplier!
    F6 lb, ub;
    I6 nside = get_bounds(cart_resl, ori_nside, cart_bound, lb, ub);
    grid_.init(nside, lb, ub);
  }

  I6 get_bounds(Float cart_resl, int ori_nside, float cart_bound, F6 &lb,
                F6 &ub) {
    I6 nside;
    if (2 * (int)(cart_bound / cart_resl) > 8192) {
      throw std::out_of_range("can have at most 8192 cart cells!");
    }
    nside[0] = nside[1] = nside[2] = 2.0 * cart_bound / cart_resl;
    nside[3] = nside[4] = nside[5] = ori_nside + 1;
    lb[0] = lb[1] = lb[2] = -cart_bound;
    ub[0] = ub[1] = ub[2] = cart_bound;
    lb[3] = lb[4] = lb[5] = -1.0 / ori_nside;
    ub[3] = ub[4] = ub[5] = 1.0;
    return nside;
  }

  XformHash_bt24_BCC6(Float cart_resl, int ori_nside, Float cart_bound) {
    // std::cout << "ori_nside c'tor" << std::endl;
    init2(cart_resl, ori_nside, cart_bound);
  }

  F6 xfrom_to_F6(Xform x, Key &cell_index) const {
    Eigen::Matrix<Float, 3, 3> rotation = x.linear();
    Eigen::Quaternion<Float> q(rotation);
    numeric::get_cell_48cell_half(q.coeffs(), cell_index);
    q = nest::pmap::hbt24_cellcen<Float>(cell_index).inverse() * q;
    q = numeric::to_half_cell(q);
    F3 params(
        params[0] = q.x() / q.w() / nest::pmap::cell_width<Float>() + 0.5,
        params[1] = q.y() / q.w() / nest::pmap::cell_width<Float>() + 0.5,
        params[2] = q.z() / q.w() / nest::pmap::cell_width<Float>() + 0.5);
    assert(cell_index < 24);
    clamp01(params);
    F6 params6;
    for (int i = 0; i < 3; ++i) {
      params6[i] = x.translation()[i];
      params6[i + 3] = params[i];
    }
    return params6;
  }

  Xform F6_to_xform(F6 params6, Key cell_index) const {
    F3 params = params6.template last<3>();
    Float const &w(nest::pmap::cell_width<Float>());
    clamp01(params);
    params = w * (params - 0.5);  // now |params| < sqrt(2)-1
    Eigen::Quaternion<Float> q(1.0, params[0], params[1], params[2]);
    q.normalize();
    q = nest::pmap::hbt24_cellcen<Float>(cell_index) * q;
    Xform center(q.matrix());
    for (int i = 0; i < 3; ++i) center.translation()[i] = params6[i];
    return center;
  }

  Key get_key(Xform const &x) const {
    Key cell_index;
    F6 params6 = xfrom_to_F6(x, cell_index);
    return cell_index << 59 | grid_[params6];
  }

  Xform get_center(Key key) const {
    Key cell_index = key >> 59;
    F6 params6 = grid_[key & (((Key)1 << 59) - (Key)1)];
    return F6_to_xform(params6, cell_index);
  }

  Key approx_size() const { return grid_.size() * 24; }

  Key approx_nori() const {
    static int const nori[18] = {192,   648,   1521,  2855,   4990,   7917,
                                 11682, 16693, 23011, 30471,  39504,  50464,
                                 62849, 77169, 93903, 112604, 133352, 157103};
    return nori[grid_.nside_[3] - 2];  // -1 for 0-index, -1 for ori_side+1
  }
};

template <class _Xform>
struct XformHash_bt24_BCC6_Phi : public XformHash_bt24_BCC6<_Xform> {
  using Key = uint64_t;
  using Xform = _Xform;
  using Float = typename Xform::Scalar;
  using Scalar = typename Xform::Scalar;
  using OriMap = rif::nest::pmap::TetracontoctachoronMap<>;
  using Grid7 = rif::numeric::BCC<7, Float, uint64_t>;
  using F3 = rif::util::SimpleArray<3, Float>;
  using I3 = rif::util::SimpleArray<3, uint64_t>;
  using F6 = rif::util::SimpleArray<6, Float>;
  using I6 = rif::util::SimpleArray<6, uint64_t>;
  using F7 = rif::util::SimpleArray<7, Float>;
  using I7 = rif::util::SimpleArray<7, uint64_t>;
  Grid7 grid7_;
  float phi_resl_;
  static std::string name() { return "XformHash_bt24_BCC6_Phi"; }
  XformHash_bt24_BCC6_Phi() {}
  XformHash_bt24_BCC6_Phi(Float phi_resl, Float cart_resl, Float ang_resl,
                          Float cart_bound = 256.0) {
    this->init(cart_resl, ang_resl, cart_bound);
    phi_resl_ = phi_resl_;
  }

  void init2(Float cart_resl, int ori_nside, Float cart_bound) {
    std::cout << "init2 subclass" << std::endl;
    cart_resl /= sqrt(3.0) / 2.0;  // TODO: HACK multiplier!
    F6 lb, ub;
    I6 ns = get_bounds(cart_resl, ori_nside, cart_bound, lb, ub);
    grid7_.init(concat(ns, uint64_t(360.0 / phi_resl_)), concat(lb, -180.0),
                concat(ub, 180.0));
  }

  Key get_key(Xform const &x, float phi) const {
    assert(phi <= 180.0);
    assert(phi >= -180.0);
    Key cell_index;
    F6 params6 = xfrom_to_F6(x, cell_index);
    F7 params7 = params6;
    params7[6] = phi;
    return cell_index << 59 | grid7_[params7];
  }

  std::pair<Xform, float> get_center(Key key) const {
    Key cell_index = key >> 59;
    F7 params7 = grid7_[key & (((Key)1 << 59) - (Key)1)];
    F6 params6 = params7;
    return std::make_pair(F6_to_xform(params6, cell_index), params7[6]);
  }
};

template <class _Xform>
struct XformHash_Quat_BCC7_Zorder {
  typedef uint64_t Key;
  typedef _Xform Xform;
  typedef typename Xform::Scalar Float;
  typedef rif::nest::pmap::TetracontoctachoronMap<> OriMap;
  typedef rif::numeric::BCC<7, Float, uint64_t> Grid;
  typedef rif::util::SimpleArray<7, Float> F7;
  typedef rif::util::SimpleArray<7, uint64_t> I7;

  static Key const ORI_MASK = ~BOOST_BINARY(
      11111111 11111111 11111000 01110000 11100001 11000011 10000111 00001110);
  static Key const ORI_MASK_NO0 = ~BOOST_BINARY(
      11111111 11111111 11111000 01110000 11100001 11000011 10000111 00001111);
  static Key const CART_MASK_NO0 = BOOST_BINARY(
      11111111 11111111 11111000 01110000 11100001 11000011 10000111 00001110);

  Grid grid_;

  bool operator==(XformHash_Quat_BCC7_Zorder<Xform> const &o) const {
    return grid_ == o.grid_;
  }
  bool operator!=(XformHash_Quat_BCC7_Zorder<Xform> const &o) const {
    return grid_ != o.grid_;
  }

  static std::string name() { return "XformHash_Quat_BCC7_Zorder"; }

  Float cart_spacing() const { return grid_.width_[0]; }

  XformHash_Quat_BCC7_Zorder() {}

  XformHash_Quat_BCC7_Zorder(Float cart_resl, Float ang_resl,
                             Float cart_bound = 512.0) {
    init(cart_resl, cart_resl, cart_bound);
  }

  XformHash_Quat_BCC7_Zorder(Float cart_resl, int ori_nside,
                             Float cart_bound = 512.0) {
    init_nside(ori_nside, cart_resl, cart_bound);
  }

  void init(Float cart_resl, Float ang_resl, Float cart_bound) {
    // bcc orientation grid covering radii
    static float const covrad[61] = {
        84.09702, 54.20621, 43.98427, 31.58683, 27.58101, 22.72314, 20.42103,
        17.58167, 16.12208, 14.44320, 13.40178, 12.15213, 11.49567, 10.53203,
        10.11448, 9.32353,  8.89083,  8.38516,  7.95147,  7.54148,  7.23572,
        6.85615,  6.63594,  6.35606,  6.13243,  5.90677,  5.72515,  5.45705,
        5.28864,  5.06335,  4.97668,  4.78774,  4.68602,  4.51794,  4.46654,
        4.28316,  4.20425,  4.08935,  3.93284,  3.84954,  3.74505,  3.70789,
        3.58776,  3.51407,  3.45023,  3.41919,  3.28658,  3.24700,  3.16814,
        3.08456,  3.02271,  2.96266,  2.91052,  2.86858,  2.85592,  2.78403,
        2.71234,  2.69544,  2.63151,  2.57503,  2.59064};
    uint64_t ori_nside = 1;
    while (covrad[ori_nside - 1] * 1.35 > ang_resl && ori_nside < 61)
      ++ori_nside;  // TODO: fix this number!
    init_nside((int)ori_nside, cart_resl, cart_bound);
  }

  void init_nside(int ori_nside, Float cart_resl, Float cart_bound) {
    cart_resl /= sqrt(3) / 2.0;
    if (2.0 * cart_bound / cart_resl > 8192.0)
      throw std::out_of_range("too many cart cells, > 8192");
    if (ori_nside > 62) throw std::out_of_range("too many ori cells");

    // cart_grid_.init(  I3(2.0*cart_bound/cart_resl), F3(-cart_bound),
    // F3(cart_bound) );
    // ori_grid_.init(  I3(ori_nside+2), F3(-1.0/ori_nside),
    // F3(1.0+1.0/ori_nside) );
    I7 nside;
    nside[0] = nside[1] = nside[2] = 2 * (int)(cart_bound / cart_resl) + 1;
    if (ori_nside % 2 == 0)
      nside[0] = nside[1] = nside[2] =
          nside[0] + 1;  // allows repr of 0,0,0... sometimes
    nside[3] = nside[4] = nside[5] = nside[6] = ori_nside + 1;
    F7 lb, ub;
    lb[0] = lb[1] = lb[2] = -cart_bound;
    ub[0] = ub[1] = ub[2] = cart_bound;
    lb[3] = lb[4] = lb[5] = lb[6] = -1.0 - 2.0 / ori_nside;
    ub[3] = ub[4] = ub[5] = ub[6] = 1.0;
    grid_.init(nside, lb, ub);
    // std::cout << "NSIDE " << nside << std::endl;
  }

  Key get_key(Xform const &x) const {
    Eigen::Matrix<Float, 3, 3> rotation = x.linear();
    Eigen::Quaternion<Float> q(rotation);
    q = numeric::to_half_cell(q);
    F7 f7;
    f7[0] = x.translation()[0];
    f7[1] = x.translation()[1];
    f7[2] = x.translation()[2];
    f7[3] = q.w();
    f7[4] = q.x();
    f7[5] = q.y();
    f7[6] = q.z();
    // std::cout << f7 << std::endl;
    bool odd;
    I7 i7 = grid_.get_indices(f7, odd);
    // std::cout << std::endl << (i7[0]>>6) << " " << (i7[1]>>6) << " " <<
    // (i7[2]>>6) << " " << i7[3] << " " << i7[4] << " "
    // << i7[5] << " " << i7[6] << " " << std::endl;
    // std::cout << std::endl << i7 << std::endl;
    Key key = odd;
    key = key | (i7[0] >> 6) << 57;
    key = key | (i7[1] >> 6) << 50;
    key = key | (i7[2] >> 6) << 43;
    key = key | util::dilate<7>(i7[0] & 63) << 1;
    key = key | util::dilate<7>(i7[1] & 63) << 2;
    key = key | util::dilate<7>(i7[2] & 63) << 3;
    key = key | util::dilate<7>(i7[3]) << 4;
    key = key | util::dilate<7>(i7[4]) << 5;
    key = key | util::dilate<7>(i7[5]) << 6;
    key = key | util::dilate<7>(i7[6]) << 7;
    return key;
  }

  I7 get_indices(Key key, bool &odd) const {
    odd = key & (Key)1;
    I7 i7;
    i7[0] = (util::undilate<7>(key >> 1) & 63) | ((key >> 57) & 127) << 6;
    i7[1] = (util::undilate<7>(key >> 2) & 63) | ((key >> 50) & 127) << 6;
    i7[2] = (util::undilate<7>(key >> 3) & 63) | ((key >> 43) & 127) << 6;
    i7[3] = util::undilate<7>(key >> 4) & 63;
    i7[4] = util::undilate<7>(key >> 5) & 63;
    i7[5] = util::undilate<7>(key >> 6) & 63;
    i7[6] = util::undilate<7>(key >> 7) & 63;
    return i7;
  }

  Xform get_center(Key key) const {
    bool odd;
    I7 i7 = get_indices(key, odd);
    // std::cout << i7 << std::endl << std::endl;
    F7 f7 = grid_.get_center(i7, odd);
    // std::cout << f7 << std::endl;
    Eigen::Quaternion<Float> q(f7[3], f7[4], f7[5], f7[6]);
    q.normalize();
    // q = numeric::to_half_cell(q); // should be un-necessary
    Xform center(q.matrix());
    center.translation()[0] = f7[0];
    center.translation()[1] = f7[1];
    center.translation()[2] = f7[2];
    return center;
  }

  Key cart_shift_key(Key key, int dx, int dy, int dz, Key d_o = 0) const {
    Key o = key % 2;
    Key x = (util::undilate<7>(key >> 1) & 63) | ((key >> 57) & 127) << 6;
    Key y = (util::undilate<7>(key >> 2) & 63) | ((key >> 50) & 127) << 6;
    Key z = (util::undilate<7>(key >> 3) & 63) | ((key >> 43) & 127) << 6;
    x += dx;
    y += dy;
    z += dz;
    o ^= d_o;
    key &= ORI_MASK;  // zero cart parts of key
    key &= ~(Key)1;   // zero even/odd
    key |= o;
    key |= (x >> 6) << 57 | util::dilate<7>(x & 63) << 1;
    key |= (y >> 6) << 50 | util::dilate<7>(y & 63) << 2;
    key |= (z >> 6) << 43 | util::dilate<7>(z & 63) << 3;
    return key;
  }

  Key approx_size() const { return grid_.size(); }

  Key approx_nori() const {
    static int const nori[63] = {
        0,      53,     53,     181,    321,    665,    874,    1642,   1997,
        2424,   3337,   4504,   5269,   6592,   8230,   10193,  11420,  14068,
        16117,  19001,  21362,  25401,  29191,  33227,  37210,  41454,  45779,
        51303,  57248,  62639,  69417,  76572,  83178,  92177,  99551,  108790,
        117666, 127850, 138032, 149535, 159922, 171989, 183625, 196557, 209596,
        226672, 239034, 253897, 271773, 288344, 306917, 324284, 342088, 364686,
        381262, 405730, 427540, 450284, 472265, 498028, 521872, 547463};
    return nori[grid_.nside_[3] - 2];  // -1 for 0-index, -1 for ori_side+1
  }

  Float cart_width() const { return grid_.width_[0]; }
  Float ang_width() const { return grid_.width_[3]; }

  Key asym_key(Key key, Key &isym) const {
    // get o,w,x,y,z for orig key
    Key o = key & (Key)1;
    Key w = util::undilate<7>(key >> 4) & 63;
    Key x = util::undilate<7>(key >> 5) & 63;
    Key y = util::undilate<7>(key >> 6) & 63;
    Key z = util::undilate<7>(key >> 7) & 63;
    // std::cout << grid_.nside_[3]-o <<" " << o << "    " << w << "\t" << x
    // <<
    // "\t" << y << "\t" << z << std::endl;
    // move to primaary
    assert(o == 0 || o == 1);
    Key nside1 = grid_.nside_[3] - o;
    isym = /*(w>nside1/2)<<3 |*/ (x <= nside1 / 2 && nside1 != 2 * x) << 2 |
           (y <= nside1 / 2 && nside1 != 2 * y) << 1 |
           (z <= nside1 / 2 && nside1 != 2 * z) << 0;
    assert(w >= nside1 / 2);
    // w = w > nside1/2 ? w : nside1-w;

    // std::cout << "FX " << (x <= nside1/2) << " " << (2*x==nside1) <<
    // std::endl;
    // std::cout << "FY " << (y <= nside1/2) << " " << (2*y==nside1) <<
    // std::endl;
    // std::cout << "FZ " << (z <= nside1/2) << " " << (2*z==nside1) <<
    // std::endl;

    x = x > nside1 / 2 ? x : nside1 - x;
    y = y > nside1 / 2 ? y : nside1 - y;
    z = z > nside1 / 2 ? z : nside1 - z;
    // make new key
    Key k = key & ~ORI_MASK;
    k = k | o;
    k = k | util::dilate<7>(w) << 4;
    k = k | util::dilate<7>(x) << 5;
    k = k | util::dilate<7>(y) << 6;
    k = k | util::dilate<7>(z) << 7;
    return k;
  }
  Key sym_key(Key key, Key isym) const {
    // get o,w,x,y,z for orig key
    Key o = key & (Key)1;
    Key w = util::undilate<7>(key >> 4) & 63;
    Key x = util::undilate<7>(key >> 5) & 63;
    Key y = util::undilate<7>(key >> 6) & 63;
    Key z = util::undilate<7>(key >> 7) & 63;
    // std::cout << grid_.nside_[3]-o <<" " << o << "    " << w << "\t" << x
    // <<
    // "\t" << y << "\t" << z << std::endl;
    // move to isym
    Key nside1 = grid_.nside_[3] - o;
    assert(w >= nside1 / 2);
    // w = (isym>>3)&1 ? w : nside1-w;
    x = (isym >> 2) & 1 ? nside1 - x : x;
    y = (isym >> 1) & 1 ? nside1 - y : y;
    z = (isym >> 0) & 1 ? nside1 - z : z;
    // make new key
    Key k = key & ~ORI_MASK;
    k = k | o;
    k = k | util::dilate<7>(w) << 4;
    k = k | util::dilate<7>(x) << 5;
    k = k | util::dilate<7>(y) << 6;
    k = k | util::dilate<7>(z) << 7;
    return k;
  }
  Key num_key_symmetries() const { return 16; }

  void print_key(Key key) const {
    bool odd;
    I7 i7 = get_indices(key, odd);
    std::cout << i7 << " " << odd << std::endl;
  }

  F7 lever_coord(Key key, Float lever_dist, F7 const &ref) const {
    bool odd;
    I7 i7 = get_indices(key, odd);
    F7 f7 = grid_.get_center(i7, odd);
    Eigen::Quaternion<Float> q(f7[3], f7[4], f7[5], f7[6]);
    q.normalize();
    Eigen::Quaternion<Float> qref(ref[3], ref[4], ref[5], ref[6]);
    bool const neg = q.dot(qref) < 0.0;
    f7[3] = (neg ? -q.w() : q.w()) * lever_dist * 2.0;
    f7[4] = (neg ? -q.x() : q.x()) * lever_dist * 2.0;
    f7[5] = (neg ? -q.y() : q.y()) * lever_dist * 2.0;
    f7[6] = (neg ? -q.z() : q.z()) * lever_dist * 2.0;
    return f7;
  }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// EVERYTHING BELOW THIS POINT MAY BE INCOMPLETE IN SOME WAY
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class _Xform>
struct XformHash_Quat_BCC7 {
  typedef uint64_t Key;
  typedef _Xform Xform;
  typedef typename Xform::Scalar Float;
  typedef rif::nest::pmap::TetracontoctachoronMap<> OriMap;
  typedef rif::numeric::BCC<7, Float, uint64_t> Grid;
  typedef rif::util::SimpleArray<7, Float> F7;
  typedef rif::util::SimpleArray<7, uint64_t> I7;

  Float grid_size_;
  Float grid_spacing_;
  OriMap ori_map_;
  Grid grid_;

  static std::string name() { return "XformHash_Quat_BCC7"; }

  XformHash_Quat_BCC7(Float cart_resl, Float ang_resl,
                      Float cart_bound = 512.0) {
    cart_resl /= sqrt(3) / 2.0;  // TODO: fix this number!
    // bcc orientation grid covering radii
    static float const covrad[99] = {
        84.09702, 54.20621, 43.98427, 31.58683, 27.58101, 22.72314, 20.42103,
        17.58167, 16.12208, 14.44320, 13.40178, 12.15213, 11.49567, 10.53203,
        10.11448, 9.32353,  8.89083,  8.38516,  7.95147,  7.54148,  7.23572,
        6.85615,  6.63594,  6.35606,  6.13243,  5.90677,  5.72515,  5.45705,
        5.28864,  5.06335,  4.97668,  4.78774,  4.68602,  4.51794,  4.46654,
        4.28316,  4.20425,  4.08935,  3.93284,  3.84954,  3.74505,  3.70789,
        3.58776,  3.51407,  3.45023,  3.41919,  3.28658,  3.24700,  3.16814,
        3.08456,  3.02271,  2.96266,  2.91052,  2.86858,  2.85592,  2.78403,
        2.71234,  2.69544,  2.63151,  2.57503,  2.59064,  2.55367,  2.48010,
        2.41046,  2.40289,  2.36125,  2.33856,  2.29815,  2.26979,  2.21838,
        2.19458,  2.17881,  2.12842,  2.14030,  2.06959,  2.05272,  2.04950,
        2.00790,  1.96385,  1.96788,  1.91474,  1.90942,  1.90965,  1.85602,
        1.83792,  1.81660,  1.80228,  1.77532,  1.76455,  1.72948,  1.72179,
        1.68324,  1.67009,  1.67239,  1.64719,  1.63832,  1.60963,  1.60093,
        1.58911};
    uint64_t ori_nside = 1;
    while (covrad[ori_nside - 1] * 1.35 > ang_resl && ori_nside < 100)
      ++ori_nside;

    if (2.0 * cart_bound / cart_resl > 8192.0)
      throw std::out_of_range("too many cart cells, > 8192");

    // cart_grid_.init(  I3(2.0*cart_bound/cart_resl), F3(-cart_bound),
    // F3(cart_bound) );
    // ori_grid_.init(  I3(ori_nside+2), F3(-1.0/ori_nside),
    // F3(1.0+1.0/ori_nside) );
    I7 nside;
    nside[0] = nside[1] = nside[2] = 2.0 * cart_bound / cart_resl;
    nside[3] = nside[4] = nside[5] = nside[6] = ori_nside + 2;
    F7 ub;
    ub[0] = ub[1] = ub[2] = cart_bound;
    ub[3] = ub[4] = ub[5] = ub[6] = 1.0 + 2.0 / ori_nside;
    grid_.init(nside, -ub, ub);
  }

  Key get_key(Xform const &x) const {
    Eigen::Matrix<Float, 3, 3> rotation = x.linear();
    Eigen::Quaternion<Float> q(rotation);
    q = numeric::to_half_cell(q);
    F7 f7;
    f7[0] = x.translation()[0];
    f7[1] = x.translation()[1];
    f7[2] = x.translation()[2];
    f7[3] = q.w();
    f7[4] = q.x();
    f7[5] = q.y();
    f7[6] = q.z();
    // std::cout << f7 << std::endl;
    Key key = grid_[f7];
    return key;
  }

  Xform get_center(Key key) const {
    F7 f7 = grid_[key];
    // std::cout << f7 << std::endl;
    Eigen::Quaternion<Float> q(f7[3], f7[4], f7[5], f7[6]);
    q.normalize();
    Xform center(q.matrix());
    center.translation()[0] = f7[0];
    center.translation()[1] = f7[1];
    center.translation()[2] = f7[2];
    return center;
  }

  Key approx_size() const { return grid_.size(); }

  Key approx_nori() const {
    static int const nori[63] = {
        0,      53,     53,     181,    321,    665,    874,    1642,   1997,
        2424,   3337,   4504,   5269,   6592,   8230,   10193,  11420,  14068,
        16117,  19001,  21362,  25401,  29191,  33227,  37210,  41454,  45779,
        51303,  57248,  62639,  69417,  76572,  83178,  92177,  99551,  108790,
        117666, 127850, 138032, 149535, 159922, 171989, 183625, 196557, 209596,
        226672, 239034, 253897, 271773, 288344, 306917, 324284, 342088, 364686,
        381262, 405730, 427540, 450284, 472265, 498028, 521872, 547463};
    return nori[grid_.nside_[3] - 2];  // -1 for 0-index, -1 for ori_side+1
  }

  Float ang_width() const { return grid_.with_[3]; }
};

template <class Xform>
struct XformHash_bt24_BCC3_Zorder {
  typedef uint64_t Key;
  typedef typename Xform::Scalar Float;
  typedef rif::nest::pmap::TetracontoctachoronMap<> OriMap;
  typedef rif::numeric::BCC<3, Float, uint64_t> Grid;
  typedef rif::util::SimpleArray<3, Float> F3;
  typedef rif::util::SimpleArray<3, uint64_t> I3;

  Float grid_size_;
  Float grid_spacing_;
  OriMap ori_map_;
  Grid cart_grid_, ori_grid_;

  static std::string name() { return "XformHash_bt24_BCC3_Zorder"; }

  XformHash_bt24_BCC3_Zorder(Float cart_resl, Float ang_resl,
                             Float cart_bound = 512.0) {
    // bcc orientation grid covering radii
    static float const covrad[64] = {
        49.66580, 25.99805, 17.48845, 13.15078, 10.48384, 8.76800, 7.48210,
        6.56491,  5.84498,  5.27430,  4.78793,  4.35932,  4.04326, 3.76735,
        3.51456,  3.29493,  3.09656,  2.92407,  2.75865,  2.62890, 2.51173,
        2.39665,  2.28840,  2.19235,  2.09949,  2.01564,  1.94154, 1.87351,
        1.80926,  1.75516,  1.69866,  1.64672,  1.59025,  1.54589, 1.50077,
        1.46216,  1.41758,  1.38146,  1.35363,  1.31630,  1.28212, 1.24864,
        1.21919,  1.20169,  1.17003,  1.14951,  1.11853,  1.09436, 1.07381,
        1.05223,  1.02896,  1.00747,  0.99457,  0.97719,  0.95703, 0.93588,
        0.92061,  0.90475,  0.89253,  0.87480,  0.86141,  0.84846, 0.83677,
        0.82164};
    uint64_t ori_nside = 1;
    while (covrad[ori_nside - 1] * 1.01 > ang_resl && ori_nside < 62)
      ++ori_nside;
    init(cart_resl, ori_nside, cart_bound);
  }

  XformHash_bt24_BCC3_Zorder(Float cart_resl, int ori_nside, Float cart_bound) {
    init(cart_resl, ori_nside, cart_bound);
  }

  void init(Float cart_resl, int ori_nside, Float cart_bound) {
    cart_resl /= 0.56;  // TODO: fix this number!
    // std::cout << "requested ang_resl: " << ang_resl << " got " <<
    // covrad[ori_nside-1] << std::endl;
    if (2 * (int)(cart_bound / cart_resl) > 8192) {
      throw std::out_of_range("can have at most 8192 cart cells!");
    }
    cart_grid_.init(I3(2.0 * cart_bound / cart_resl), F3(-cart_bound),
                    F3(cart_bound));
    ori_grid_.init(I3(ori_nside + 2), F3(-1.0 / ori_nside),
                   F3(1.0 + 1.0 / ori_nside));

    // std::cout << "cart_bcc " << cart_grid_ << std::endl;
    // std::cout << " ori_bcc " <<  ori_grid_ << std::endl;
    // std::cout << log2( (double)cart_grid_.size() / ori_grid_.size() ) <<
    // std::endl;
  }

  // Key structure
  // 5 bits bt24 cell index
  // 7 bits high order cart X bits
  // 7 bits high order cart Y bits
  // 7 bits high order cart Z bits
  // 36 bits 6*6 zorder cart/ori
  // 2 bits cart/ori even/odd
  Key get_key(Xform const &x) const {
    Eigen::Matrix3d rotation = x.linear();

    uint64_t cell_index;
    F3 params;
    ori_map_.value_to_params(rotation, 0, params, cell_index);
    assert(cell_index < 24);
    assert(0.0 <= params[0] && params[0] <= 1.0);
    assert(0.0 <= params[1] && params[1] <= 1.0);
    assert(0.0 <= params[2] && params[2] <= 1.0);

    bool ori_odd, cart_odd;
    I3 ori_indices, cart_indices;
    ori_indices = ori_grid_.get_indices(params, ori_odd);
    F3 trans(x.translation());
    cart_indices = cart_grid_.get_indices(trans, cart_odd);

    // std::cout << "get_index  " << cell_index << " " << cart_indices << " "
    // <<
    // cart_odd << " " << ori_indices << " " << ori_odd << std::endl;

    Key key;

    key = cell_index << 59;

    key = key | (cart_indices[0] >> 6) << 52;
    key = key | (cart_indices[1] >> 6) << 45;
    key = key | (cart_indices[2] >> 6) << 38;

    // 6*6 zorder
    key = key >> 2;
    key = key | (util::dilate<6>(ori_indices[0]) << 0);
    key = key | (util::dilate<6>(ori_indices[1]) << 1);
    key = key | (util::dilate<6>(ori_indices[2]) << 2);
    key = key | (util::dilate<6>((cart_indices[0] & 63)) << 3);
    key = key | (util::dilate<6>((cart_indices[1] & 63)) << 4);
    key = key | (util::dilate<6>((cart_indices[2] & 63)) << 5);
    key = key << 2;

    // lowest two bits, even/odd
    key = key | ori_odd | cart_odd << 1;

    return key;
  }

  Xform get_center(Key key) const {
    I3 cart_indices, ori_indices;

    uint64_t cell_index = key >> 59;

    cart_indices[0] =
        (((key >> 52) & 127) << 6) | (util::undilate<6>((key >> 5)) & 63);
    cart_indices[1] =
        (((key >> 45) & 127) << 6) | (util::undilate<6>((key >> 6)) & 63);
    cart_indices[2] =
        (((key >> 38) & 127) << 6) | (util::undilate<6>((key >> 7)) & 63);

    ori_indices[0] =
        util::undilate<6>((key >> 2) & (((uint64_t)1 << 36) - 1)) & 63;
    ori_indices[1] =
        util::undilate<6>((key >> 3) & (((uint64_t)1 << 36) - 1)) & 63;
    ori_indices[2] =
        util::undilate<6>((key >> 4) & (((uint64_t)1 << 36) - 1)) & 63;

    bool ori_odd = key & (Key)1;
    bool cart_odd = key & (Key)2;

    F3 trans = cart_grid_.get_center(cart_indices, cart_odd);
    F3 params = ori_grid_.get_center(ori_indices, ori_odd);
    Eigen::Matrix3d m;
    ori_map_.params_to_value(params, cell_index, 0, m);
    // std::cout << "get_center " << cell_index << " " << cart_indices << " "
    // <<
    // cart_odd << " " << ori_indices << " " << ori_odd << std::endl;

    Xform center(m);
    center.translation()[0] = trans[0];
    center.translation()[1] = trans[1];
    center.translation()[2] = trans[2];

    return center;
  }

  Key approx_size() const {
    return (ori_grid_.sizes_[0] - 1) * (ori_grid_.sizes_[1] - 1) *
           (ori_grid_.sizes_[2] - 1) * 2 * cart_grid_.size() * 24;
  }

  Key approx_nori() const { throw std::logic_error("not implemented"); }
};

template <class Xform>
struct XformHash_bt24_BCC3 {
  typedef uint64_t Key;
  typedef typename Xform::Scalar Float;
  typedef rif::nest::pmap::TetracontoctachoronMap<> OriMap;
  typedef rif::numeric::BCC<3, Float, uint64_t> Grid;
  typedef rif::util::SimpleArray<3, Float> F3;
  typedef rif::util::SimpleArray<3, uint64_t> I3;

  Float grid_size_;
  Float grid_spacing_;
  OriMap ori_map_;
  Grid cart_grid_, ori_grid_;

  static std::string name() { return "XformHash_bt24_BCC3"; }

  XformHash_bt24_BCC3(Float cart_resl, Float ang_resl,
                      Float cart_bound = 512.0) {
    cart_resl /= 0.56;  // TODO: fix this number!
    // bcc orientation grid covering radii
    static float const covrad[64] = {
        49.66580, 25.99805, 17.48845, 13.15078, 10.48384, 8.76800, 7.48210,
        6.56491,  5.84498,  5.27430,  4.78793,  4.35932,  4.04326, 3.76735,
        3.51456,  3.29493,  3.09656,  2.92407,  2.75865,  2.62890, 2.51173,
        2.39665,  2.28840,  2.19235,  2.09949,  2.01564,  1.94154, 1.87351,
        1.80926,  1.75516,  1.69866,  1.64672,  1.59025,  1.54589, 1.50077,
        1.46216,  1.41758,  1.38146,  1.35363,  1.31630,  1.28212, 1.24864,
        1.21919,  1.20169,  1.17003,  1.14951,  1.11853,  1.09436, 1.07381,
        1.05223,  1.02896,  1.00747,  0.99457,  0.97719,  0.95703, 0.93588,
        0.92061,  0.90475,  0.89253,  0.87480,  0.86141,  0.84846, 0.83677,
        0.82164};
    uint64_t ori_nside = 1;
    while (covrad[ori_nside - 1] * 1.01 > ang_resl && ori_nside < 62)
      ++ori_nside;
    // std::cout << "requested ang_resl: " << ang_resl << " got " <<
    // covrad[ori_nside-1] << std::endl;
    if (2 * (int)(cart_bound / cart_resl) > 8192) {
      throw std::out_of_range("can have at most 8192 cart cells!");
    }
    cart_grid_.init(I3(2.0 * cart_bound / cart_resl), F3(-cart_bound),
                    F3(cart_bound));
    ori_grid_.init(I3(ori_nside + 2), F3(-1.0 / ori_nside),
                   F3(1.0 + 1.0 / ori_nside));

    // std::cout << "cart_bcc " << cart_grid_ << std::endl;
    // std::cout << " ori_bcc " <<  ori_grid_ << std::endl;
    // std::cout << log2( (double)cart_grid_.size() / ori_grid_.size() ) <<
    // std::endl;
  }

  Key get_key(Xform const &x) const {
    Eigen::Matrix3d rotation = x.linear();

    uint64_t cell_index;
    F3 params;
    ori_map_.value_to_params(rotation, 0, params, cell_index);
    assert(cell_index < 24);
    assert(0.0 <= params[0] && params[0] <= 1.0);
    assert(0.0 <= params[1] && params[1] <= 1.0);
    assert(0.0 <= params[2] && params[2] <= 1.0);

    F3 trans(x.translation());
    Key ori_index, cart_index;
    ori_index = ori_grid_[params];
    cart_index = cart_grid_[trans];
    // std::cout << cart_index << " " << ori_isndex << " " << cell_index <<
    // std::endl;
    Key key;

    key = cell_index << 59;
    key = key | cart_index << 18;
    key = key | ori_index;

    return key;
  }

  Xform get_center(Key key) const {
    I3 cart_indices, ori_indices;

    Key cell_index = key >> 59;
    Key cart_index = (key << 5) >> 23;
    Key ori_index = key & (((Key)1 << 18) - 1);
    // std::cout << cart_index << " " << ori_index << " " << cell_index <<
    // std::endl;

    F3 trans = cart_grid_[cart_index];
    F3 params = ori_grid_[ori_index];
    Eigen::Matrix3d m;
    ori_map_.params_to_value(params, cell_index, 0, m);

    Xform center(m);
    center.translation()[0] = trans[0];
    center.translation()[1] = trans[1];
    center.translation()[2] = trans[2];

    return center;
  }

  Key approx_size() const {
    return (ori_grid_.sizes_[0] - 1) * (ori_grid_.sizes_[1] - 1) *
           (ori_grid_.sizes_[2] - 1) * 2 * cart_grid_.size() * 24;
  }

  Key approx_nori() const { throw std::logic_error("not implemented"); }
};

template <class Xform>
struct XformHash_bt24_Cubic_Zorder {
  typedef uint64_t Key;
  typedef typename Xform::Scalar Float;
  typedef rif::nest::pmap::TetracontoctachoronMap<> OriMap;
  typedef rif::numeric::Cubic<3, Float, uint64_t> Grid;
  typedef rif::util::SimpleArray<3, Float> F3;
  typedef rif::util::SimpleArray<3, uint64_t> I3;

  Float grid_size_;
  Float grid_spacing_;
  OriMap ori_map_;
  Grid cart_grid_, ori_grid_;

  static std::string name() { return "XformHash_bt24_Cubic_Zorder"; }

  XformHash_bt24_Cubic_Zorder(Float cart_resl, Float ang_resl,
                              Float cart_bound = 512.0) {
    cart_resl /= 0.867;  // TODO: fix this number!
    // bcc orientation grid covering radii
    static float const covrad[64] = {
        62.71876, 39.26276, 26.61019, 20.06358, 16.20437, 13.45733, 11.58808,
        10.10294, 9.00817,  8.12656,  7.37295,  6.74856,  6.23527,  5.77090,
        5.38323,  5.07305,  4.76208,  4.50967,  4.25113,  4.04065,  3.88241,
        3.68300,  3.53376,  3.36904,  3.22018,  3.13437,  2.99565,  2.89568,
        2.78295,  2.70731,  2.61762,  2.52821,  2.45660,  2.37996,  2.31057,
        2.25207,  2.18726,  2.13725,  2.08080,  2.02489,  1.97903,  1.92123,
        1.88348,  1.83759,  1.79917,  1.76493,  1.72408,  1.68516,  1.64581,
        1.62274,  1.57909,  1.55846,  1.52323,  1.50846,  1.47719,  1.44242,
        1.42865,  1.39023,  1.37749,  1.34783,  1.32588,  1.31959,  1.29872,
        1.26796};
    uint64_t ori_nside = 1;
    while (covrad[ori_nside - 1] * 1.01 > ang_resl && ori_nside < 62)
      ++ori_nside;
    // std::cout << "requested ang_resl: " << ang_resl << " got " <<
    // covrad[ori_nside-1] << std::endl;
    if (2 * (int)(cart_bound / cart_resl) > 8192) {
      throw std::out_of_range("can have at most 8192 cart cells!");
    }
    cart_grid_.init(I3(2.0 * cart_bound / cart_resl), F3(-cart_bound),
                    F3(cart_bound));
    ori_grid_.init(I3(ori_nside), F3(0.0), F3(1.0));

    // std::cout << "cart_bcc " << cart_grid_ << std::endl;
    // std::cout << " ori_bcc " <<  ori_grid_ << std::endl;
    // std::cout << log2( (double)cart_grid_.size() / ori_grid_.size() ) <<
    // std::endl;
  }

  Key get_key(Xform const &x) const {
    Eigen::Matrix3d rotation = x.linear();

    uint64_t cell_index;
    F3 params;
    ori_map_.value_to_params(rotation, 0, params, cell_index);
    assert(cell_index < 24);
    assert(0.0 <= params[0] && params[0] <= 1.0);
    assert(0.0 <= params[1] && params[1] <= 1.0);
    assert(0.0 <= params[2] && params[2] <= 1.0);

    I3 ori_indices, cart_indices;
    ori_indices = ori_grid_.get_indices(params);
    F3 trans(x.translation());
    cart_indices = cart_grid_.get_indices(trans);

    // std::cout << "get_index  " << cell_index << " " << cart_indices << " "
    // <<
    // ori_indices << std::endl;

    Key key;

    key = cell_index << 59;

    key = key | (cart_indices[0] >> 6) << 52;
    key = key | (cart_indices[1] >> 6) << 45;
    key = key | (cart_indices[2] >> 6) << 38;

    // 6*6 zorder
    key = key >> 2;
    key = key | (util::dilate<6>(ori_indices[0]) << 0);
    key = key | (util::dilate<6>(ori_indices[1]) << 1);
    key = key | (util::dilate<6>(ori_indices[2]) << 2);
    key = key | (util::dilate<6>((cart_indices[0] & 63)) << 3);
    key = key | (util::dilate<6>((cart_indices[1] & 63)) << 4);
    key = key | (util::dilate<6>((cart_indices[2] & 63)) << 5);
    key = key << 2;

    return key;
  }

  Xform get_center(Key key) const {
    I3 cart_indices, ori_indices;

    uint64_t cell_index = key >> 59;

    cart_indices[0] =
        (((key >> 52) & 127) << 6) | (util::undilate<6>((key >> 5)) & 63);
    cart_indices[1] =
        (((key >> 45) & 127) << 6) | (util::undilate<6>((key >> 6)) & 63);
    cart_indices[2] =
        (((key >> 38) & 127) << 6) | (util::undilate<6>((key >> 7)) & 63);

    ori_indices[0] = util::undilate<6>((key >> 2) & (((Key)1 << 36) - 1)) & 63;
    ori_indices[1] = util::undilate<6>((key >> 3) & (((Key)1 << 36) - 1)) & 63;
    ori_indices[2] = util::undilate<6>((key >> 4) & (((Key)1 << 36) - 1)) & 63;

    F3 trans = cart_grid_.get_center(cart_indices);
    F3 params = ori_grid_.get_center(ori_indices);
    // std::cout << "get_center " << cell_index << " " << cart_indices << " "
    // <<
    // ori_indices << std::endl;
    Eigen::Matrix3d m;
    ori_map_.params_to_value(params, cell_index, 0, m);

    Xform center(m);
    center.translation()[0] = trans[0];
    center.translation()[1] = trans[1];
    center.translation()[2] = trans[2];

    return center;
  }

  Key approx_size() const {
    return (ori_grid_.sizes_[0]) * (ori_grid_.sizes_[1]) *
           (ori_grid_.sizes_[2]) * cart_grid_.size() * 24;
  }

  Key approx_nori() const { throw std::logic_error("not implemented"); }
};

template <class Xform>
struct XformHash_Quatgrid_Cubic {
  typedef uint64_t Key;
  typedef typename Xform::Scalar Float;
  typedef rif::nest::pmap::TetracontoctachoronMap<> OriMap;
  typedef rif::numeric::Cubic<3, Float, uint64_t> Grid;
  typedef rif::util::SimpleArray<3, Float> F3;
  typedef rif::util::SimpleArray<3, uint64_t> I3;

  Float grid_size_;
  Float grid_spacing_;
  OriMap ori_map_;
  Grid cart_grid_, ori_grid_;

  static std::string name() { return "XformHash_Quatgrid_Cubic"; }

  XformHash_Quatgrid_Cubic(Float cart_resl, Float ang_resl,
                           Float cart_bound = 512.0) {
    cart_resl /= 0.56;  // TODO: fix this number!
    // bcc orientation grid covering radii
    static float const covrad[64] = {
        49.66580, 25.99805, 17.48845, 13.15078, 10.48384, 8.76800, 7.48210,
        6.56491,  5.84498,  5.27430,  4.78793,  4.35932,  4.04326, 3.76735,
        3.51456,  3.29493,  3.09656,  2.92407,  2.75865,  2.62890, 2.51173,
        2.39665,  2.28840,  2.19235,  2.09949,  2.01564,  1.94154, 1.87351,
        1.80926,  1.75516,  1.69866,  1.64672,  1.59025,  1.54589, 1.50077,
        1.46216,  1.41758,  1.38146,  1.35363,  1.31630,  1.28212, 1.24864,
        1.21919,  1.20169,  1.17003,  1.14951,  1.11853,  1.09436, 1.07381,
        1.05223,  1.02896,  1.00747,  0.99457,  0.97719,  0.95703, 0.93588,
        0.92061,  0.90475,  0.89253,  0.87480,  0.86141,  0.84846, 0.83677,
        0.82164};
    uint64_t ori_nside = 1;
    while (covrad[ori_nside - 1] > ang_resl && ori_nside < 62) ++ori_nside;
    // std::cout << "requested ang_resl: " << ang_resl << " got " <<
    // covrad[ori_nside-1] << std::endl;
    if (2 * (int)(cart_bound / cart_resl) > 8192) {
      throw std::out_of_range("can have at most 8192 cart cells!");
    }
    cart_grid_.init(I3(2.0 * cart_bound / cart_resl), F3(-cart_bound),
                    F3(cart_bound));
    ori_grid_.init(I3(ori_nside + 2), F3(-1.0 / ori_nside),
                   F3(1.0 + 1.0 / ori_nside));

    // std::cout << "cart_bcc " << cart_grid_ << std::endl;
    // std::cout << " ori_bcc " <<  ori_grid_ << std::endl;
    // std::cout << log2( (double)cart_grid_.size() / ori_grid_.size() ) <<
    // std::endl;
  }

  Key get_key(Xform const &x) const {
    Eigen::Matrix3d rotation = x.linear();

    uint64_t cell_index;
    F3 params;
    ori_map_.value_to_params(rotation, 0, params, cell_index);
    assert(cell_index < 24);
    assert(0.0 <= params[0] && params[0] <= 1.0);
    assert(0.0 <= params[1] && params[1] <= 1.0);
    assert(0.0 <= params[2] && params[2] <= 1.0);

    bool ori_odd, cart_odd;
    I3 ori_indices, cart_indices;
    ori_indices = ori_grid_.get_indices(params, ori_odd);
    F3 trans(x.translation());
    cart_indices = cart_grid_.get_indices(trans, cart_odd);

    // std::cout << "get_index  " << cell_index << " " << cart_indices << " "
    // <<
    // cart_odd << " " << ori_indices << " " << ori_odd << std::endl;

    Key key;

    key = cell_index << 59;

    key = key | (cart_indices[0] >> 6) << 52;
    key = key | (cart_indices[1] >> 6) << 45;
    key = key | (cart_indices[2] >> 6) << 38;

    // 6*6 zorder
    key = key >> 2;
    key = key | (util::dilate<6>(ori_indices[0]) << 0);
    key = key | (util::dilate<6>(ori_indices[1]) << 1);
    key = key | (util::dilate<6>(ori_indices[2]) << 2);
    key = key | (util::dilate<6>((cart_indices[0] & 63)) << 3);
    key = key | (util::dilate<6>((cart_indices[1] & 63)) << 4);
    key = key | (util::dilate<6>((cart_indices[2] & 63)) << 5);
    key = key << 2;

    // lowest two bits, even/odd
    key = key | ori_odd | cart_odd << 1;

    return key;
  }

  Xform get_center(Key key) const {
    I3 cart_indices, ori_indices;

    uint64_t cell_index = key >> 59;

    cart_indices[0] =
        (((key >> 52) & 127) << 6) | (util::undilate<6>((key >> 5)) & 63);
    cart_indices[1] =
        (((key >> 45) & 127) << 6) | (util::undilate<6>((key >> 6)) & 63);
    cart_indices[2] =
        (((key >> 38) & 127) << 6) | (util::undilate<6>((key >> 7)) & 63);

    ori_indices[0] = util::undilate<6>((key >> 2) & (((Key)1 << 36) - 1)) & 63;
    ori_indices[1] = util::undilate<6>((key >> 3) & (((Key)1 << 36) - 1)) & 63;
    ori_indices[2] = util::undilate<6>((key >> 4) & (((Key)1 << 36) - 1)) & 63;

    bool ori_odd = key & (Key)1;
    bool cart_odd = key & (Key)2;

    F3 trans = cart_grid_.get_center(cart_indices, cart_odd);
    F3 params = ori_grid_.get_center(ori_indices, ori_odd);
    Eigen::Matrix3d m;
    ori_map_.params_to_value(params, cell_index, 0, m);
    // std::cout << "get_center " << cell_index << " " << cart_indices << " "
    // <<
    // cart_odd << " " << ori_indices << " " << ori_odd << std::endl;

    Xform center(m);
    center.translation()[0] = trans[0];
    center.translation()[1] = trans[1];
    center.translation()[2] = trans[2];

    return center;
  }

  Key approx_size() const { return ori_grid_.size() * cart_grid_.size() * 24; }

  Key approx_nori() const { throw std::logic_error("not implemented"); }
};
}
}
