#ifndef INCLUDED_scheme_nest_maps_HecatonicosachoronMap_HH
#define INCLUDED_scheme_nest_maps_HecatonicosachoronMap_HH

#include "numeric/geom_4d.hpp"

#include <Eigen/Dense>
#include "util/SimpleArray.hpp"

#include <boost/static_assert.hpp>
#include <iostream>
#include <vector>

namespace scheme {
namespace nest {
namespace pmap {

using namespace Eigen;
using std::cout;
using std::endl;

template <class Float>
Float cell_width() {
  return 2.0 * sqrt(2.0) - 2.0;
}

template <class Float, class Index>
Eigen::Map<Eigen::Quaternion<Float> const> hbt24_cellcen(Index const &i) {
  // Float const * tmp = numeric::get_raw_48cell_half<Float>() + 4*i;
  // std::cout << "   raw hbt24_cellcen " << tmp[0] << " " << tmp[1] << " " <<
  // tmp[2] << " " << tmp[3] << std::endl;
  return Eigen::Map<Eigen::Quaternion<Float> const>(
      numeric::get_raw_48cell_half<Float>() + 4 * i);
}

template <int DIM = 3, class Value = Eigen::Matrix3d, class Index = uint64_t,
          class Float = typename Value::Scalar>
struct TetracontoctachoronMap {
  BOOST_STATIC_ASSERT_MSG(DIM == 3, "TetracontoctachoronMap DIM must be == 3");

  static std::string pmap_name() { return "TetracontoctachoronMap"; }

  static int const DIMENSION = DIM;
  typedef Value ValueType;
  typedef Float FloatType;
  typedef Index IndexType;
  typedef util::SimpleArray<DIM, Index> Indices;
  typedef util::SimpleArray<DIM, Float> Params;

  Index nside_;
  Float one_over_nside_;

  TetracontoctachoronMap() { init(1); }
  TetracontoctachoronMap(Index nside) { init(nside); }

  void init(uint64_t nside) {
    nside_ = nside;
    one_over_nside_ = 1.0 / nside_;
  }

  ///@brief sets value to parameters without change
  ///@return false iff invalid parameters
  bool params_to_value(Params const &params, Index cell_index, Index resl,
                       Value &value) const {
    // // cout << "        set p0 " << params << endl;
    Float const &w(cell_width<Float>());

    Params p = params * one_over_nside_;

    Index h48_cell_index = cell_index / (nside_ * nside_ * nside_);
    cell_index = cell_index % (nside_ * nside_ * nside_);
    p[0] += one_over_nside_ * (Float)(cell_index % nside_);
    p[1] += one_over_nside_ * (Float)(cell_index / nside_ % nside_);
    p[2] += one_over_nside_ * (Float)(cell_index / (nside_ * nside_) % nside_);

    // if( !( p[0] >= 0.0 && p[0] <= 1.0 ) ) cout << "BAD param val: " << p[0]
    // << endl;
    // if( !( p[1] >= 0.0 && p[1] <= 1.0 ) ) cout << "BAD param val: " << p[1]
    // << endl;
    // if( !( p[2] >= 0.0 && p[2] <= 1.0 ) ) cout << "BAD param val: " << p[2]
    // << endl;

    assert(p[0] >= -0.00001 && p[0] <= 1.00001);
    assert(p[1] >= -0.00001 && p[1] <= 1.00001);
    assert(p[2] >= -0.00001 && p[2] <= 1.00001);
    p[0] = fmax(0.0, p[0]);
    p[1] = fmax(0.0, p[1]);
    p[2] = fmax(0.0, p[2]);
    p[0] = fmin(1.0, p[0]);
    p[1] = fmin(1.0, p[1]);
    p[2] = fmin(1.0, p[2]);

    // std::cout << cell_index << " " << p << " " << p << std::endl;
    // static int count = 0; if( ++count > 30 ) std::exit(-1);

    p = w * (p - 0.5);  // now |p| < sqrt(2)-1

    // if( resl > 3 ){
    Float corner_dist = fabs(p[0]) + fabs(p[1]) + fabs(p[2]);
    Float delta = sqrt(3.0) / 2.0 / w / (Float)(1 << resl);
    // // static int count = 0;
    // //          std::cout << corner_dist << "    " << p << " " << p <<
    // std::endl;
    //          if(++count > 100) std::exit(-1);
    if (corner_dist - delta > 1.0)
      return false;  // TODO make this check more rigerous???
    // }

    // Eigen::Quaternion<Float> q( sqrt(1.0-p.squaredNorm()), p[0], p[1], p[2]
    // );
    // assert( fabs(q.squaredNorm()-1.0) < 0.000001 );
    Eigen::Quaternion<Float> q(1.0, p[0], p[1], p[2]);
    q.normalize();

    q = hbt24_cellcen<Float>(h48_cell_index) * q;

    value = q.matrix();

    return true;
  }

  ///@brief sets params/cell_index from value
  ///@note necessary for value lookup and neighbor lookup
  bool value_to_params(Value const &value, Index /*resl*/, Params &params,
                       Index &cell_index) const {
    Quaternion<Float> q(value);
    // q = numeric::to_half_cell(q);
    // // cout << "get q  " << q.coeffs().transpose() << endl;

    Index h48_cell_index;
    numeric::get_cell_48cell_half(q.coeffs(), h48_cell_index);

    q = hbt24_cellcen<Float>(h48_cell_index).inverse() * q;

    q = numeric::to_half_cell(q);

    // // cout << q.w() << endl;
    // assert( q.w() > 0.7 );

    // // cout << "    get q0 " << q.coeffs().transpose() << endl;

    // // cout << "      get p  " << q.x() << " " << q.y() << " " << q.z() <<
    // endl;

    params[0] = q.x() / q.w() / cell_width<Float>() + 0.5;
    params[1] = q.y() / q.w() / cell_width<Float>() + 0.5;
    params[2] = q.z() / q.w() / cell_width<Float>() + 0.5;

    Indices ci = params * nside_;
    cell_index = ci[0] + ci[1] * nside_ + ci[2] * nside_ * nside_;

    params = params * nside_ - ci.template cast<Float>();
    assert(params[0] >= 0.0 && params[0] <= 1.0);
    assert(params[1] >= 0.0 && params[1] <= 1.0);
    assert(params[2] >= 0.0 && params[2] <= 1.0);

    // just for testing...
    // params[0] = fmax(0.000001,params[0]);
    // params[1] = fmax(0.000001,params[1]);
    // params[2] = fmax(0.000001,params[2]);
    // params[0] = fmin(0.999999,params[0]);
    // params[1] = fmin(0.999999,params[1]);
    // params[2] = fmin(0.999999,params[2]);

    cell_index += h48_cell_index * nside_ * nside_ * nside_;

    // cout << "        get p0 " << params << endl;

    return true;
  }

  ///@brief get parameter space repr of Value for particular cell
  ///@note necessary only for neighbor lookup
  void value_to_params_for_cell(Value const &value, Params &params) const {
    std::cerr << "Not Implemented" << std::endl;
    std::exit(-1);
  }

  ///@brief return the cell_index of neighboring cells within radius of value
  ///@note delta parameter is in "Parameter Space"
  template <class OutIter>
  void get_neighboring_cells(Value const &value, Float radius,
                             OutIter out) const {
    std::cerr << "Not Implemented" << std::endl;
    std::exit(-1);
  }

  static Float const *get_covrad_data() {
    static Float const covrad[25] = {
        62.76235,  // 1
        38.63604,  // 2
        26.71264,  // 3
        20.62393,  // 4
        17.02567,  // 5
        14.25487,  // 6
        12.42992,  // 7
        11.02897,  // 8
        9.62588,   // 9
        8.70544,   // 10
        7.82964,   // 11
        7.28521,   // 12
        6.62071,   // 13
        6.13243,   // 14
        5.81918,   // 15
        5.44871,   // 16
        5.14951,   // 17
        4.82331,   // 18
        4.52938,   // 19
        4.31905,   // 20
        4.07469,   // 21
        3.93772,   // 22
        3.77275,   // 23
        3.64786,   // 24
        3.44081    // 25
    };
    return covrad;
  }

  static int get_nside_for_rot_resl_deg(Float rot_resl_deg) {
    static Float const *covrad = get_covrad_data();
    int nside = 0;
    while (covrad[nside] > rot_resl_deg && nside < 23) {
      // std::cout << nside << " " << covrad[nside] << std::endl;
      ++nside;
    }
    return nside + 1;
  }

  ///@brief aka covering radius max distance from bin center to any value within
  /// bin
  Float bin_circumradius(Index resl) const {
    BOOST_VERIFY(resl < 6);
    if (resl == 0) {
      static Float const *covrad = get_covrad_data();
      if (nside_ > 25) {
        std::cerr
            << "TetracontoctachoronMap::bin_circumradius > 25 not implemented"
            << std::endl;
        std::exit(-1);
      }
      return covrad[nside_ - 1] * M_PI / 180.0;
    }
    assert(nside_ == 1);
    static Float const covrad[6] = {62.76235, 37.95720, 20.53126,
                                    11.00814, 5.31355,  2.66953};
    return covrad[resl] * M_PI / 180.0;
  }

  ///@brief maximum distance from the bin center which must be within the bin
  Float bin_inradius(Index resl) const {
    // double const delta = 1.0/(double)(1ul<<resl);
    // Vec pworst = Vec(1,1,1);// - delta*Vec(2.0,0,0);
    // Vec pNest0 = Vec(1,1,1) - delta*Vec(1.0,1.0,0);
    // cube_to_sphere(pworst);
    // cube_to_sphere(pNest0);
    // pworst.normalize();
    // pNest0.normalize();
    // return pworst.distance(pNest0) * 0.5999; // should be half of
    // curcumradius based on geometry
  }

  ///@brief cell size
  Index num_cells() const { return 24 * nside_ * nside_ * nside_; }
};

template <int DIM, class Value, class Index, class Float>
std::ostream &operator<<(
    std::ostream &out,
    TetracontoctachoronMap<DIM, Value, Index, Float> const &tm) {
  out << "TetracontoctachoronMap nside = " << tm.nside_
      << " covrad0 = " << tm.bin_circumradius(0) * 180.0 / M_PI;
  return out;
}
}
}
}

#endif
