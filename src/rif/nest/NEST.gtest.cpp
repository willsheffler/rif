#include <gtest/gtest.h>

#include <boost/assign/std/vector.hpp>  // for 'operator+=()'
#include <boost/foreach.hpp>
#include <iterator>
#include <random>
#include "nest/NEST.hpp"
#include "nest/NEST_concepts.hpp"
#include "nest/NEST_test_util.hpp"
#include "nest/pmap/DiscreteChoiceMap.hpp"
#include "nest/pmap/UnitMap.hpp"
#include "util/meta/util.hpp"

// #include <Eigen/Core>

namespace rif {
namespace nest {

using std::cout;
using std::endl;

using namespace pmap;

TEST(NEST, concpets) {
  NEST<1, concept::ValueArchitype, concept::ParamMapArchitype>().set_state(1,
                                                                           0);
  // NEST< 1, double, concept::ParamMapArchitype >().get_index(1.0,0); // this
  // should not compile!
  NEST<1, concept::ValueArchitype, concept::ParamMapInvertableArchitype>()
      .set_state(1, 0);
  NEST<1, concept::ValueArchitype, concept::ParamMapInvertableArchitype>()
      .get_index(concept::ValueArchitype(), 0);
  // NEST< 1, concept::ValueArchitype, UnitMap >().set_state(1,0);
  NEST<1, concept::ArrayValueArchitype<1>, UnitMap>().set_state(1, 0);
  NEST<1, concept::ArrayValueArchitype<1>, UnitMap>().get_index(
      concept::ArrayValueArchitype<1>(), 0);
  SUCCEED();
}

TEST(NEST, particular_cases) {
  NEST<1> nest;
  ASSERT_EQ(nest.set_and_get(0, 0)[0], 0.5);
  ASSERT_EQ(nest.set_and_get(0, 1)[0], 0.25);
  ASSERT_EQ(nest.set_and_get(1, 1)[0], 0.75);
  ASSERT_EQ(nest.set_and_get(0, 2)[0], 0.125);
  ASSERT_EQ(nest.set_and_get(1, 2)[0], 0.375);
  ASSERT_EQ(nest.set_and_get(2, 2)[0], 0.625);
  ASSERT_EQ(nest.set_and_get(3, 2)[0], 0.875);
  NEST<2> nest2;
  ASSERT_EQ(nest2.set_and_get(0, 0), NEST<2>::ValueType(0.5, 0.5));
  ASSERT_EQ(nest2.set_and_get(0, 1), NEST<2>::ValueType(0.25, 0.25));
  ASSERT_EQ(nest2.set_and_get(1, 1), NEST<2>::ValueType(0.75, 0.25));
  ASSERT_EQ(nest2.set_and_get(2, 1), NEST<2>::ValueType(0.25, 0.75));
  ASSERT_EQ(nest2.set_and_get(3, 1), NEST<2>::ValueType(0.75, 0.75));
}

// // keep this out because #include <Eigen/Core> takes almost a second to
// compile
// TEST(NEST,works_with_eigen){
//  typedef Eigen::Matrix<float,2,1> MAT;
//  typedef NEST<2,MAT> Nest;
//  Nest nest;
//  MAT m; m << 0.5,0.5;
//  ASSERT_EQ( nest.set_and_get(0,0), m );
// }

TEST(NEST, map_discrete) {
  using std::cout;
  using std::endl;

  using namespace boost::assign;
  std::vector<double> choices;
  choices += 42.0, 152.345, 8049782.83402;
  NEST<0, double, pmap::DiscreteChoiceMap, StoreValue> nest0(choices);
  ASSERT_EQ(choices.size(), nest0.size());
  ASSERT_EQ(choices.size(), nest0.size(0));
  ASSERT_EQ(choices.size(), nest0.size(3));
  ASSERT_EQ(choices.size(), nest0.size(4));
  for (size_t i = 0; i < nest0.size(); ++i) {
    ASSERT_TRUE(nest0.set_state(i));
    ASSERT_EQ(choices[i], nest0.set_and_get(i));
    ASSERT_EQ(choices[i], nest0.set_and_get(i, 0));
    ASSERT_EQ(choices[i], nest0.set_and_get(i, 7));
    ASSERT_EQ(choices[i], nest0.set_and_get(i, 3));
    ASSERT_EQ(choices[i], nest0.set_and_get(i, 8));
  }
  ASSERT_FALSE(nest0.set_state(5));
}

template <int DIM>
void test_coverage_unit1cell() {
  std::mt19937 rng((unsigned int)time(0));
  std::uniform_real_distribution<> uniform;
  NEST<DIM> nest;
  size_t rmax = 9 / DIM + 1;
  for (size_t r = 0; r <= rmax; ++r) {
    double cellradius = sqrt(DIM) * 0.5 / (1 << r);
    for (size_t iter = 0; iter < 100000 / DIM; ++iter) {
      util::SimpleArray<DIM, double> tgt;
      for (size_t i = 0; i < DIM; ++i) tgt[i] = uniform(rng);
      size_t index = nest.get_index(tgt, r);
      ASSERT_LT(index, nest.size(r));
      util::SimpleArray<DIM, double> val = nest.set_and_get(index, r);
      ASSERT_LE((tgt - val).norm(), cellradius);
    }
    util::SimpleArray<DIM, double> zeros;
    zeros.fill(0);
    ASSERT_LE((zeros - nest.set_and_get(nest.get_index(zeros, r), r)).norm(),
              cellradius);
    util::SimpleArray<DIM, double> ones;
    ones.fill(0.999999);
    ASSERT_LE((ones - nest.set_and_get(nest.get_index(ones, r), r)).norm(),
              cellradius);
  }
}

TEST(NEST, coverage_cube) {
  test_coverage_unit1cell<1>();
  test_coverage_unit1cell<2>();
  test_coverage_unit1cell<3>();
  test_coverage_unit1cell<4>();
  test_coverage_unit1cell<5>();
  test_coverage_unit1cell<6>();
}

TEST(NEST, index_nesting) {
  generic_test_index_nesting_of_bincenters(NEST<1>());
  generic_test_index_nesting_of_bincenters(NEST<2>());
  generic_test_index_nesting_of_bincenters(NEST<3>());
  generic_test_index_nesting_of_bincenters(NEST<4>());
  generic_test_index_nesting_of_bincenters(NEST<5>());
  generic_test_index_nesting_of_bincenters(NEST<6>());
}

template <int DIM>
void test_store_pointer_virtual() {
  typedef util::SimpleArray<DIM, double> VAL;
  NEST<DIM, VAL, UnitMap, StoreValue, size_t, double, true> nest_val(2);
  NEST<DIM, VAL, UnitMap, StorePointer, size_t, double, true> nest_ptr(2);
  NestBase<size_t> *nest_val_virtual = &nest_val;
  NestBase<size_t> *nest_ptr_virtual = &nest_ptr;
  StorePointer<VAL> *ptr_store =
      dynamic_cast<StorePointer<VAL> *>(nest_ptr_virtual);
  StorePointer<VAL> *ptr_store_fail =
      dynamic_cast<StorePointer<VAL> *>(nest_val_virtual);
  ASSERT_TRUE(ptr_store);
  ASSERT_FALSE(ptr_store_fail);
  VAL val_pointed_to;
  ptr_store->set_pointer(&val_pointed_to);
  size_t rmax = 9 / DIM + 1;
  for (size_t r = 0; r <= rmax; ++r) {
    for (size_t i = 0; i < nest_val.size(r); ++i) {
      VAL tmp;
      boost::any a = &tmp;
      ASSERT_TRUE(nest_val_virtual->virtual_get_state(i, r, a));
      ASSERT_TRUE(nest_ptr_virtual->virtual_get_state(i, r, a));
      ASSERT_EQ(nest_val.value(), val_pointed_to);
      ASSERT_EQ(nest_val.value(), tmp);
      ASSERT_EQ(tmp, val_pointed_to);
    }
  }
}
TEST(NEST, store_pointer_virtual) {
  test_store_pointer_virtual<1>();
  test_store_pointer_virtual<2>();
  test_store_pointer_virtual<3>();
  test_store_pointer_virtual<4>();
  test_store_pointer_virtual<5>();
  test_store_pointer_virtual<6>();
}

TEST(NEST, store_nothing) {
  typedef util::SimpleArray<2, double> VAL;
  NEST<2, VAL, UnitMap, util::StoreNothing, size_t, double, true> nest;
  VAL tmp(-10, -10);
  boost::any a = &tmp;
  ASSERT_TRUE(nest.virtual_get_state(0, 0, a));
  ASSERT_EQ(tmp, VAL(0.5, 0.5));

  std::vector<size_t> vi(2);
  vi[0] = 1;
  vi[1] = 0;
  size_t ii = 0;
  ASSERT_TRUE(nest.virtual_get_state(vi, 0, ii, 1, a));
  ASSERT_EQ(tmp, VAL(0.75, 0.25));
  ASSERT_EQ(ii, 2);
}

template <int DIM>
void test_uniformity() {
  typedef util::SimpleArray<DIM, double> VAL;
  NEST<DIM, VAL, UnitMap, StoreValue> nest(1);
  size_t rmax = 9 / DIM + 1;
  for (size_t r = 0; r <= rmax; ++r) {
    double scale = 1 << r;
    std::vector<int> counts(1 << r, 0);
    for (size_t i = 0; i < nest.size(r); ++i) {
      ASSERT_TRUE(nest.set_state(i, r));
      for (size_t j = 0; j < DIM; ++j) {
        counts[nest.value()[j] * scale]++;
      }
    }
    ASSERT_EQ((1 << ((DIM - 1) * r)) * DIM,
              *std::min_element(counts.begin(), counts.end()));
    ASSERT_EQ((1 << ((DIM - 1) * r)) * DIM,
              *std::max_element(counts.begin(), counts.end()));
  }
}

TEST(NEST, uniformity) {
  test_uniformity<1>();
  test_uniformity<2>();
  test_uniformity<3>();
  test_uniformity<4>();
  test_uniformity<5>();
  test_uniformity<6>();
}

TEST(NEST, bounds) {
  NEST<1> nest;
  ASSERT_TRUE(nest.set_state(0, 0));
  ASSERT_FALSE(nest.set_state(1, 0));
#ifndef NDEBUG
#ifndef CXX14
  ASSERT_DEATH(nest.set_and_get(1, 0), "Assertion.*failed");
#endif
#endif
}

TEST(NEST, virtual_get_index) {
  std::mt19937 rng((unsigned int)time(0));
  std::uniform_real_distribution<> uniform;

  typedef util::SimpleArray<2, double> VAL;
  NEST<2, VAL, UnitMap, util::StoreNothing, uint64_t, double, true> nest;
  NestBase<uint64_t> *nestp = &nest;

  boost::any a;
  for (int resl = 0; resl < 6; ++resl) {
    for (int iter = 0; iter < 10000; ++iter) {
      VAL v(uniform(rng), uniform(rng));
      a = &v;
      ASSERT_EQ(nest.get_index(v, resl), nestp->virtual_get_index(a, resl));
    }
  }

  NEST<2, VAL, nest::concept::ParamMapArchitype> nest2;
#ifndef CXX14
  ASSERT_DEATH(nest2.virtual_get_index(a, 0), ".*");
#endif
}
}
}
