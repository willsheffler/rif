#include <gtest/gtest.h>
#include <boost/assign/std/vector.hpp>  // for 'operator+=()'
#include <boost/foreach.hpp>
#include <random>
#include "nest/NEST.hpp"
#include "nest/NEST_test_util.hpp"
#include "nest/pmap/DiscreteChoiceMap.hpp"
#include "nest/pmap/ScaleMap.hpp"
#include "nest/pmap/UnitMap.hpp"

namespace rif {
namespace nest {
namespace pmap {

using std::cout;
using std::endl;

using rif::nest::StorePointer;

TEST(NEST_ScaleMap, particular_values) {
  {
    typedef util::SimpleArray<1, double> VAL;
    VAL lb, ub;
    util::SimpleArray<1, size_t> bs;
    lb.fill(-16.0);
    ub.fill(16.0);
    bs.fill(1);
    NEST<1, VAL, ScaleMap> nest1(lb, ub, bs);
    bs.fill(2);
    NEST<1, VAL, ScaleMap> nest2(lb, ub, bs);
    ASSERT_EQ(nest1.set_and_get(0, 0)[0], 0.0);
    ASSERT_EQ(nest1.set_and_get(0, 1)[0], -8.0);
    ASSERT_EQ(nest1.set_and_get(1, 1)[0], 8.0);

    ASSERT_EQ(nest2.set_and_get(0, 0)[0], -8.0);
    ASSERT_EQ(nest2.set_and_get(1, 0)[0], 8.0);
    ASSERT_EQ(nest2.set_and_get(0, 1)[0], -12.0);
    ASSERT_EQ(nest2.set_and_get(1, 1)[0], -4.0);
    ASSERT_EQ(nest2.set_and_get(2, 1)[0], 4.0);
    ASSERT_EQ(nest2.set_and_get(3, 1)[0], 12.0);
  }
  {
    typedef util::SimpleArray<2, double> VAL;
    util::SimpleArray<2, double> lb(-16.0, -24), ub(16.0, 24);
    util::SimpleArray<2, size_t> bs(2, 3);
    NEST<2, VAL, ScaleMap> nest(lb, ub, bs);
    ASSERT_EQ(nest.set_and_get(0, 0), VAL(-8, -16));
    ASSERT_EQ(nest.set_and_get(1, 0), VAL(8, -16));
    ASSERT_EQ(nest.set_and_get(2, 0), VAL(-8, 0));
    ASSERT_EQ(nest.set_and_get(3, 0), VAL(8, 0));
    ASSERT_EQ(nest.set_and_get(4, 0), VAL(-8, 16));
    ASSERT_EQ(nest.set_and_get(5, 0), VAL(8, 16));
    ASSERT_EQ(nest.set_and_get(0, 1), VAL(-12, -20));
    ASSERT_EQ(nest.set_and_get(1, 1), VAL(-4, -20));
    ASSERT_EQ(nest.set_and_get(2, 1), VAL(-12, -12));
    ASSERT_EQ(nest.set_and_get(3, 1), VAL(-4, -12));
    ASSERT_EQ(nest.set_and_get(0, 2), VAL(-14, -22));
    ASSERT_EQ(nest.set_and_get(1, 2), VAL(-10, -22));
    ASSERT_EQ(nest.set_and_get(2, 2), VAL(-14, -18));
    ASSERT_EQ(nest.set_and_get(3, 2), VAL(-10, -18));
    ASSERT_EQ(nest.set_and_get(0, 3), VAL(-15, -23));
    ASSERT_EQ(nest.set_and_get(1, 3), VAL(-13, -23));
    ASSERT_EQ(nest.set_and_get(2, 3), VAL(-15, -21));
    ASSERT_EQ(nest.set_and_get(3, 3), VAL(-13, -21));
    ASSERT_EQ(nest.set_and_get(5 * 64 + 60, 3), VAL(13, 21));
    ASSERT_EQ(nest.set_and_get(5 * 64 + 61, 3), VAL(15, 21));
    ASSERT_EQ(nest.set_and_get(5 * 64 + 62, 3), VAL(13, 23));
    ASSERT_EQ(nest.set_and_get(5 * 64 + 63, 3), VAL(15, 23));
    // cout << VAL(0,0) << endl;
  }
}

TEST(UnitMap, value_to_params_for_cell) {
  typedef UnitMap<2> MapType;
  MapType umap(3);
  MapType::Params params;

  umap.value_to_params_for_cell(MapType::ValueType(1.5, 0.5), 0, params, 0);
  ASSERT_EQ(params[0], 1.5);
  ASSERT_EQ(params[1], 0.5);

  umap.value_to_params_for_cell(MapType::ValueType(1.5, 0.5), 0, params, 1);
  ASSERT_EQ(params[0], 0.5);
  ASSERT_EQ(params[1], 0.5);

  umap.value_to_params_for_cell(MapType::ValueType(1.5, 0.5), 0, params, 2);
  ASSERT_EQ(params[0], -0.5);
  ASSERT_EQ(params[1], 0.5);

  umap.value_to_params_for_cell(MapType::ValueType(0.2, 0.5), 0, params, 2);
  ASSERT_EQ(params[0], -1.8);
  ASSERT_EQ(params[1], 0.5);
}

TEST(ScaleMap, value_to_params_for_cell) {
  {
    typedef ScaleMap<2> MapType;
    MapType smap(MapType::Params(0, 0), MapType::Params(4, 4),
                 MapType::Indices(4, 4));

    MapType::Params params;

    smap.value_to_params_for_cell(
        MapType::ValueType(1.5, 0.5), 0, params,
        smap.indices_to_cellindex(MapType::Indices(0, 0)));
    ASSERT_EQ(params[0], 1.5);
    ASSERT_EQ(params[1], 0.5);

    smap.value_to_params_for_cell(
        MapType::ValueType(1.5, 0.5), 0, params,
        smap.indices_to_cellindex(MapType::Indices(1, 0)));
    ASSERT_EQ(params[0], 0.5);
    ASSERT_EQ(params[1], 0.5);

    smap.value_to_params_for_cell(
        MapType::ValueType(1.5, 0.5), 0, params,
        smap.indices_to_cellindex(MapType::Indices(0, 1)));
    ASSERT_EQ(params[0], 1.5);
    ASSERT_EQ(params[1], -0.5);

    smap.value_to_params_for_cell(
        MapType::ValueType(1.5, 0.5), 0, params,
        smap.indices_to_cellindex(MapType::Indices(2, 0)));
    ASSERT_EQ(params[0], -0.5);
    ASSERT_EQ(params[1], 0.5);

    smap.value_to_params_for_cell(
        MapType::ValueType(1.5, 0.5), 0, params,
        smap.indices_to_cellindex(MapType::Indices(1, 2)));
    ASSERT_EQ(params[0], 0.5);
    ASSERT_EQ(params[1], -1.5);

    smap.value_to_params_for_cell(
        MapType::ValueType(1.5, 0.5), 0, params,
        smap.indices_to_cellindex(MapType::Indices(3, 3)));
    ASSERT_EQ(params[0], -1.5);
    ASSERT_EQ(params[1], -2.5);
  }
  {
    typedef ScaleMap<2> MapType;
    double scale = 2.0;
    MapType::Params shift(0.597, 1.1243);
    MapType smap(MapType::Params(0.0, 0.0) * scale - shift,
                 MapType::Params(4.0, 4.0) * scale - shift,
                 MapType::Indices(4, 4));

    MapType::Params params;

    smap.value_to_params_for_cell(
        MapType::ValueType(1.5 * scale - shift[0], 0.5 * scale - shift[1]), 0,
        params, smap.indices_to_cellindex(MapType::Indices(0, 0)));
    ASSERT_EQ(params[0], 1.5);
    ASSERT_EQ(params[1], 0.5);

    smap.value_to_params_for_cell(
        MapType::ValueType(1.5 * scale - shift[0], 0.5 * scale - shift[1]), 0,
        params, smap.indices_to_cellindex(MapType::Indices(1, 0)));
    ASSERT_EQ(params[0], 0.5);
    ASSERT_EQ(params[1], 0.5);

    smap.value_to_params_for_cell(
        MapType::ValueType(1.5 * scale - shift[0], 0.5 * scale - shift[1]), 0,
        params, smap.indices_to_cellindex(MapType::Indices(0, 1)));
    ASSERT_EQ(params[0], 1.5);
    ASSERT_EQ(params[1], -0.5);

    smap.value_to_params_for_cell(
        MapType::ValueType(1.5 * scale - shift[0], 0.5 * scale - shift[1]), 0,
        params, smap.indices_to_cellindex(MapType::Indices(2, 0)));
    ASSERT_EQ(params[0], -0.5);
    ASSERT_EQ(params[1], 0.5);

    smap.value_to_params_for_cell(
        MapType::ValueType(1.5 * scale - shift[0], 0.5 * scale - shift[1]), 0,
        params, smap.indices_to_cellindex(MapType::Indices(1, 2)));
    ASSERT_EQ(params[0], 0.5);
    ASSERT_EQ(params[1], -1.5);

    smap.value_to_params_for_cell(
        MapType::ValueType(1.5 * scale - shift[0], 0.5 * scale - shift[1]), 0,
        params, smap.indices_to_cellindex(MapType::Indices(3, 3)));
    ASSERT_EQ(params[0], -1.5);
    ASSERT_EQ(params[1], -2.5);
  }
}

template <int DIM>
void test_index_lookup_scaled() {
  typedef util::SimpleArray<DIM, double> VAL;
  util::SimpleArray<10, double> lb0;
  lb0 << -1.3, 2.2, 0, -3, -5, -9.9, 1.3, 44, -13.3, 99;
  util::SimpleArray<10, double> ub0;
  ub0 << 1.3, 4.2, 1, 10, -3, 9.9, 4.3, 44, 13.3, 199;
  util::SimpleArray<10, size_t> bs0;
  bs0 << 2, 1, 2, 3, 4, 5, 6, 7, 8, 9;

  typedef util::SimpleArray<DIM, double> VAL;
  VAL lb, ub;
  util::SimpleArray<DIM, size_t> bs;
  for (size_t i = 0; i < DIM; ++i) {
    lb[i] = lb0[i];
    ub[i] = ub0[i];
    bs[i] = bs0[i];
  }
  NEST<DIM, VAL, ScaleMap, StoreValue> nest(lb, ub, bs);
  size_t rmax = 6 / DIM;
  for (size_t r = 0; r <= rmax; ++r) {
    for (size_t i = 0; i < nest.size(r); ++i) {
      ASSERT_TRUE(nest.set_state(i, r));
      VAL value = nest.value();
      size_t index = nest.get_index(value, r);
      ASSERT_EQ(i, index);
      for (size_t r2 = 0; r2 <= r; ++r2) {
        size_t index2 = nest.get_index(value, r2);
        ASSERT_EQ(i >> (DIM * (r - r2)), index2);
      }
    }
  }
}

TEST(NEST_ScaleMap, index_lookup_scaled) {
  test_index_lookup_scaled<1>();
  test_index_lookup_scaled<2>();
  test_index_lookup_scaled<3>();
  test_index_lookup_scaled<4>();
  test_index_lookup_scaled<5>();
  test_index_lookup_scaled<6>();
}

template <int DIM>
void test_map_scale_bounds() {
  BOOST_STATIC_ASSERT((DIM < 10));
  util::SimpleArray<10, double> lb0, ub0;
  lb0 << -1.3, 2.2, 0, -3, -5, -9.9, 1.3, 44, -13.3, 99;
  ub0 << 1.3, 4.2, 1, 10, -3, 9.9, 4.3, 44, 13.3, 199;
  util::SimpleArray<10, size_t> bs0;
  bs0 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

  typedef util::SimpleArray<DIM, double> VAL;
  VAL lb, ub;
  util::SimpleArray<DIM, size_t> bs;
  for (size_t i = 0; i < DIM; ++i) {
    lb[i] = lb0[i];
    ub[i] = ub0[i];
    bs[i] = bs0[i];
  }
  NEST<DIM, VAL, ScaleMap, StoreValue> nest(lb, ub, bs);

  size_t resl = 8 / DIM;
  for (size_t i = 0; i < nest.size(resl); ++i) {
    ASSERT_TRUE(nest.set_state(i, resl));
    for (size_t j = 0; j < DIM; ++j) {
      ASSERT_LT(lb[j], nest.value()[j]);
      ASSERT_LT(nest.value()[j], ub[j]);
    }
    ASSERT_FALSE(nest.set_state(i + nest.size(resl), resl));
    for (size_t j = 0; j < DIM; ++j) {
      ASSERT_LT(lb[j], nest.value()[j]);
      ASSERT_LT(nest.value()[j], ub[j]);
    }
  }
}

TEST(NEST_ScaleMap, map_scale) {
  test_map_scale_bounds<1>();
  test_map_scale_bounds<2>();
  test_map_scale_bounds<3>();
  test_map_scale_bounds<4>();
  test_map_scale_bounds<5>();
  test_map_scale_bounds<6>();
}

template <class NEST>
void test_bin_circumradius(NEST nest, typename NEST::ValueType lb,
                           typename NEST::ValueType ub) {
  size_t NITER = 10 * 1000;
#ifdef SCHEME_BENCHMARK
  NITER *= 50;
#endif
  std::mt19937 rng((unsigned int)time(0));
  std::uniform_real_distribution<> uniform;
  typename NEST::ValueType randpt;
  for (size_t r = 0; r <= std::min((size_t)10, (size_t)NEST::MAX_RESL_ONE_CELL);
       ++r) {
    double maxdis = 0;
    for (size_t iter = 0; iter < NITER / NEST::DIMENSION; ++iter) {
      for (size_t i = 0; i < NEST::DIMENSION; ++i)
        randpt[i] = uniform(rng) * (ub[i] - lb[i]) + lb[i];
      nest.set_state(nest.get_index(randpt, r), r);
      double dist = (randpt - nest.value()).norm();
      // cout << randpt.transpose() << " " << nest.value().transpose() << endl;
      maxdis = fmax(maxdis, dist);
      ASSERT_LE(dist, nest.bin_circumradius(r));
    }
    // covering radius should be reasonably tigth
    // cout << NEST::DIMENSION << " " << r << " " << nest.bin_circumradius(r) <<
    // " " << maxdis << endl;
    if (nest.bin_circumradius(r) * (1.0 - (double)NEST::DIMENSION / 20.0) >
        maxdis) {
      cout << "WARNING(PROBABILISTIC): covering radius may be too loose DIM="
           << NEST::DIMENSION << " resl=" << r << endl;
      cout << "                        covering radius "
           << nest.bin_circumradius(r) << " max observerd: " << maxdis << endl;
      cout << "                        if you see a handful of these, don't "
              "worry. if you see lots, then worry"
           << endl;
    }
    // cout << DIM << " " << maxdis << " " << nest.bin_circumradius(r) <<
    // std::endl;
  }
}

TEST(NEST_UnitMap, NEST_bin_circumradius_unitmap) {
  {
    NEST<1>::ValueType lb, ub;
    lb.fill(0);
    ub.fill(1);
    test_bin_circumradius(NEST<1>(), lb, ub);
  }
  {
    NEST<2>::ValueType lb, ub;
    lb.fill(0);
    ub.fill(1);
    test_bin_circumradius(NEST<2>(), lb, ub);
  }
  {
    NEST<3>::ValueType lb, ub;
    lb.fill(0);
    ub.fill(1);
    test_bin_circumradius(NEST<3>(), lb, ub);
  }
  {
    NEST<4>::ValueType lb, ub;
    lb.fill(0);
    ub.fill(1);
    test_bin_circumradius(NEST<4>(), lb, ub);
  }
  {
    NEST<5>::ValueType lb, ub;
    lb.fill(0);
    ub.fill(1);
    test_bin_circumradius(NEST<5>(), lb, ub);
  }
  {
    NEST<6>::ValueType lb, ub;
    lb.fill(0);
    ub.fill(1);
    test_bin_circumradius(NEST<6>(), lb, ub);
  }
}

template <int DIM>
void test_bin_circumradius_scalemap() {
  typedef NEST<DIM, util::SimpleArray<DIM, double>, ScaleMap> NestType;
  typename NestType::ValueType lb, ub;
  for (size_t i = 0; i < DIM; ++i) {
    lb[i] = -1;
    ub[i] = 1 + 2 * i;
  }
  NestType nest(lb, ub);
  test_bin_circumradius<NestType>(nest, lb, ub);
}

TEST(NEST_ScaleMap, NEST_bin_circumradius_scalemap) {
  test_bin_circumradius_scalemap<1>();
  test_bin_circumradius_scalemap<2>();
  test_bin_circumradius_scalemap<3>();
  test_bin_circumradius_scalemap<4>();
  test_bin_circumradius_scalemap<5>();
  test_bin_circumradius_scalemap<6>();
}

template <class NestType>
void test_coverage_random_ScaleMap() {
  size_t NITER = 1 * 1000;
  double FUDGE = 0.06;
#ifdef SCHEME_BENCHMARK
  NITER *= 50;
  FUDGE = 0.03;
#endif

  std::mt19937 rng((unsigned int)time(0));
  std::uniform_real_distribution<> uniform;

  // set up random bounds
  typename NestType::Params lb, ub;
  typename NestType::Indices cs;
  for (size_t i = 0; i < NestType::DIMENSION; ++i) {
    double b1 = uniform(rng) * 20.0 - 10.0;
    double b2 = uniform(rng) * 20.0 - 10.0;
    lb[i] = std::min(b1, b2);
    ub[i] = std::max(b1, b2);
    cs[i] = (typename NestType::IndexType)(uniform(rng) * 3.999) + 1;
    assert(lb[i] < ub[i]);
    assert(cs[i] > 0);
  }
  // cout << "LB " << lb.transpose() << endl;
  // cout << "UB " << ub.transpose() << endl;
  // cout << "CS " << cs.transpose() << endl;

  NestType nest(lb, ub, cs);
  size_t max_resl = std::min(
      (size_t)10, (size_t)NestType::MAX_RESL_ONE_CELL -
                      2);  // -2 because we set cs up to 4 per dimension
  std::vector<double> largest_d2_for_r(max_resl + 1, 0.0);
  for (size_t i = 0; i < NITER; ++i) {
    // set up random value within bounds
    typename NestType::ValueType val;
    for (size_t j = 0; j < NestType::DIMENSION; ++j) {
      val[j] = (ub[j] - lb[j]) * uniform(rng) + lb[j];
      assert(lb[j] <= val[j]);
      assert(ub[j] >= val[j]);
    }

    // run generic coverage test
    generic_test_coverage_of_value(nest, val, largest_d2_for_r, max_resl);
  }
  for (size_t r = 0; r <= max_resl; ++r) {
    // cout << r << " " << largest_d2_for_r[r] << endl;
    // factor of 1.0+0.03*DIM is a *TOTAL* hack, errer does not scale this way
    // by dimension
    // but it's enough to make sure the curcumradius is reasonably tight
    ASSERT_LT(nest.bin_circumradius(r),
              (1.0 + FUDGE * (double)NestType::DIMENSION) *
                  sqrt(largest_d2_for_r[r]));
  }
}

TEST(NEST_ScaleMap, test_coverage_DIM_1_to_9) {
  test_coverage_random_ScaleMap<
      NEST<1, util::SimpleArray<1, double>, ScaleMap>>();
  test_coverage_random_ScaleMap<
      NEST<2, util::SimpleArray<2, double>, ScaleMap>>();
  test_coverage_random_ScaleMap<
      NEST<3, util::SimpleArray<3, double>, ScaleMap>>();
  test_coverage_random_ScaleMap<
      NEST<4, util::SimpleArray<4, double>, ScaleMap>>();
  test_coverage_random_ScaleMap<
      NEST<5, util::SimpleArray<5, double>, ScaleMap>>();
  test_coverage_random_ScaleMap<
      NEST<6, util::SimpleArray<6, double>, ScaleMap>>();
  test_coverage_random_ScaleMap<
      NEST<7, util::SimpleArray<7, double>, ScaleMap>>();
  test_coverage_random_ScaleMap<
      NEST<8, util::SimpleArray<8, double>, ScaleMap>>();
  test_coverage_random_ScaleMap<
      NEST<9, util::SimpleArray<9, double>, ScaleMap>>();
}
}
}
}
