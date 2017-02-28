#include "nest/NEST.hpp"
#include <gtest/gtest.h>
#include <boost/assign/std/vector.hpp>  // for 'operator+=()'
// // #include <random>
#include <boost/foreach.hpp>
#include <iterator>
#include "nest/pmap/DiscreteChoiceMap.hpp"
#include "nest/pmap/ScaleMap.hpp"
#include "nest/pmap/UnitMap.hpp"

namespace rif {
namespace nest {

using std::cout;
using std::endl;

using namespace pmap;

TEST(NEST_NEIGHBOR, dim2_test_case) {
  NEST<2> nest;
  std::vector<size_t> neighbors;
  NEST<2>::ValueType val(0.85, 0.15);
  size_t r = 1;
  std::back_insert_iterator<std::vector<size_t>> back_it(neighbors);
  nest.get_neighbors(val, r, back_it);
  ASSERT_EQ(neighbors.size(), 4);
  ASSERT_EQ(neighbors[0], 0);
  ASSERT_EQ(neighbors[1], 1);
  ASSERT_EQ(neighbors[2], 2);
  ASSERT_EQ(neighbors[3], 3);
}

TEST(NEST_NEIGHBOR, dim3_test_case) {
  NEST<3> nest;
  typedef NEST<3>::ValueType VAL;
  VAL val(0.85, 0.15, 0.5);
  // val[0] = 0.85; val[1] = 0.15; val[2] = 0.5;
  size_t r = 3;
  std::vector<size_t> neighbors;
  std::back_insert_iterator<std::vector<size_t>> back_it(neighbors);
  nest.get_neighbors(val, r, back_it);
  ASSERT_EQ(neighbors.size(), 27);
  EXPECT_EQ(nest.set_and_get(neighbors[0], r), VAL(0.6875, 0.0625, 0.4375));
  EXPECT_EQ(nest.set_and_get(neighbors[1], r), VAL(0.8125, 0.0625, 0.4375));
  EXPECT_EQ(nest.set_and_get(neighbors[2], r), VAL(0.9375, 0.0625, 0.4375));
  EXPECT_EQ(nest.set_and_get(neighbors[3], r), VAL(0.6875, 0.1875, 0.4375));
  EXPECT_EQ(nest.set_and_get(neighbors[4], r), VAL(0.8125, 0.1875, 0.4375));
  EXPECT_EQ(nest.set_and_get(neighbors[5], r), VAL(0.9375, 0.1875, 0.4375));
  EXPECT_EQ(nest.set_and_get(neighbors[6], r), VAL(0.6875, 0.3125, 0.4375));
  EXPECT_EQ(nest.set_and_get(neighbors[7], r), VAL(0.8125, 0.3125, 0.4375));
  EXPECT_EQ(nest.set_and_get(neighbors[8], r), VAL(0.9375, 0.3125, 0.4375));
  EXPECT_EQ(nest.set_and_get(neighbors[9], r), VAL(0.6875, 0.0625, 0.5625));
  EXPECT_EQ(nest.set_and_get(neighbors[10], r), VAL(0.8125, 0.0625, 0.5625));
  EXPECT_EQ(nest.set_and_get(neighbors[11], r), VAL(0.9375, 0.0625, 0.5625));
  EXPECT_EQ(nest.set_and_get(neighbors[12], r), VAL(0.6875, 0.1875, 0.5625));
  EXPECT_EQ(nest.set_and_get(neighbors[13], r), VAL(0.8125, 0.1875, 0.5625));
  EXPECT_EQ(nest.set_and_get(neighbors[14], r), VAL(0.9375, 0.1875, 0.5625));
  EXPECT_EQ(nest.set_and_get(neighbors[15], r), VAL(0.6875, 0.3125, 0.5625));
  EXPECT_EQ(nest.set_and_get(neighbors[16], r), VAL(0.8125, 0.3125, 0.5625));
  EXPECT_EQ(nest.set_and_get(neighbors[17], r), VAL(0.9375, 0.3125, 0.5625));
  EXPECT_EQ(nest.set_and_get(neighbors[18], r), VAL(0.6875, 0.0625, 0.6875));
  EXPECT_EQ(nest.set_and_get(neighbors[19], r), VAL(0.8125, 0.0625, 0.6875));
  EXPECT_EQ(nest.set_and_get(neighbors[20], r), VAL(0.9375, 0.0625, 0.6875));
  EXPECT_EQ(nest.set_and_get(neighbors[21], r), VAL(0.6875, 0.1875, 0.6875));
  EXPECT_EQ(nest.set_and_get(neighbors[22], r), VAL(0.8125, 0.1875, 0.6875));
  EXPECT_EQ(nest.set_and_get(neighbors[23], r), VAL(0.9375, 0.1875, 0.6875));
  EXPECT_EQ(nest.set_and_get(neighbors[24], r), VAL(0.6875, 0.3125, 0.6875));
  EXPECT_EQ(nest.set_and_get(neighbors[25], r), VAL(0.8125, 0.3125, 0.6875));
  EXPECT_EQ(nest.set_and_get(neighbors[26], r), VAL(0.9375, 0.3125, 0.6875));
}

TEST(NEST_NEIGHBOR, unit_1d_boundary_1cell) {
  NEST<1> nest;
  std::vector<size_t> neighbors;
  std::back_insert_iterator<std::vector<size_t>> back_it(neighbors);
  size_t r = 2;
  NEST<1>::ValueType val;

  neighbors.clear();
  val = 1.1;
  nest.get_neighbors_for_cell(val, r, 0, back_it);
  // BOOST_FOREACH(size_t i,neighbors){ cout << i << " " ; } cout << endl;
  ASSERT_EQ(neighbors.size(), 1);
  ASSERT_EQ(neighbors[0], 3);

  neighbors.clear();
  val = 0.99;
  nest.get_neighbors_for_cell(val, r, 0, back_it);
  ASSERT_EQ(neighbors.size(), 2);
  ASSERT_EQ(neighbors[0], 2);
  ASSERT_EQ(neighbors[1], 3);

  neighbors.clear();
  val = 0.6;
  nest.get_neighbors_for_cell(val, r, 0, back_it);
  ASSERT_EQ(neighbors.size(), 3);
  ASSERT_EQ(neighbors[0], 1);
  ASSERT_EQ(neighbors[1], 2);
  ASSERT_EQ(neighbors[2], 3);
}

TEST(NEST_NEIGHBOR, unit_2d_boundary_1cell) {
  typedef NEST<2> NestType;
  typedef NestType::ValueType VAL;

  NestType nest;
  std::vector<size_t> neighbors;
  std::back_insert_iterator<std::vector<size_t>> back_it(neighbors);
  size_t r = 2;

  neighbors.clear();
  nest.get_neighbors_for_cell(NestType::ValueType(1.5, 1.5), r, 0, back_it);
  ASSERT_EQ(neighbors.size(), 0);
  neighbors.clear();
  nest.get_neighbors_for_cell(NestType::ValueType(-1.5, 1.5), r, 0, back_it);
  ASSERT_EQ(neighbors.size(), 0);
  neighbors.clear();
  nest.get_neighbors_for_cell(NestType::ValueType(1.5, -1.5), r, 0, back_it);
  ASSERT_EQ(neighbors.size(), 0);
  neighbors.clear();
  nest.get_neighbors_for_cell(NestType::ValueType(-1.5, -1.5), r, 0, back_it);
  ASSERT_EQ(neighbors.size(), 0);

  neighbors.clear();
  nest.get_neighbors_for_cell(NestType::ValueType(1.1, 1.1), r, 0, back_it);
  // BOOST_FOREACH(size_t i,neighbors){ cout << i << "(" <<
  // nest.set_and_get(i,r)<<") "; } cout << endl;
  ASSERT_EQ(neighbors.size(), 1);
  ASSERT_EQ(nest.set_and_get(neighbors[0], r), VAL(0.875, 0.875));

  neighbors.clear();
  nest.get_neighbors_for_cell(NestType::ValueType(0.9, 1.1), r, 0, back_it);
  // BOOST_FOREACH(size_t i,neighbors){ cout << i << "(" <<
  // nest.set_and_get(i,r)<<") "; } cout << endl;
  ASSERT_EQ(neighbors.size(), 2);
  ASSERT_EQ(nest.set_and_get(neighbors[0], r), VAL(0.625, 0.875));
  ASSERT_EQ(nest.set_and_get(neighbors[1], r), VAL(0.875, 0.875));

  neighbors.clear();
  nest.get_neighbors_for_cell(NestType::ValueType(1.1, 0.6), r, 0, back_it);
  // BOOST_FOREACH(size_t i,neighbors){ cout << i << "(" <<
  // nest.set_and_get(i,r)<<") "; } cout << endl;
  ASSERT_EQ(neighbors.size(), 3);
  ASSERT_EQ(nest.set_and_get(neighbors[0], r), VAL(0.875, 0.375));
  ASSERT_EQ(nest.set_and_get(neighbors[1], r), VAL(0.875, 0.625));
  ASSERT_EQ(nest.set_and_get(neighbors[2], r), VAL(0.875, 0.875));

  neighbors.clear();
  nest.get_neighbors_for_cell(NestType::ValueType(-0.249, -0.249), r, 0,
                              back_it);
  // BOOST_FOREACH(size_t i,neighbors){ cout << i << "(" <<
  // nest.set_and_get(i,r)<<") "; } cout << endl;
  ASSERT_EQ(neighbors.size(), 1);
  ASSERT_EQ(nest.set_and_get(neighbors[0], r), VAL(0.125, 0.125));

  neighbors.clear();
  nest.get_neighbors_for_cell(NestType::ValueType(-0.251, -0.251), r, 0,
                              back_it);
  // BOOST_FOREACH(size_t i,neighbors){ cout << i << "(" <<
  // nest.set_and_get(i,r)<<") "; } cout << endl;
  ASSERT_EQ(neighbors.size(), 0);
}

TEST(NEST_NEIGHBOR, unitmap_neighbors_2d_boundary_2cell) {
  typedef NEST<2> NestType;
  typedef NestType::ValueType VAL;
  NestType nest(2);
  std::vector<size_t> neighbors;
  std::back_insert_iterator<std::vector<size_t>> back_it(neighbors);
  size_t r = 2;
  int i;

  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(0, 0), r, back_it);
  ASSERT_EQ(neighbors.size(), 4);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.125, 0.125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.375, 0.125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.125, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.375, 0.375));

  // cout << "size: " << neighbors.size() << endl;
  // BOOST_FOREACH(size_t i,neighbors) cout << nest.set_and_get(i,r) << " " << i
  // << endl;

  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(0.999, 0.499), r, back_it);
  ASSERT_EQ(neighbors.size(), 9);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.625, 0.125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.875, 0.125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.625, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.875, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.625, 0.625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.875, 0.625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.125, 0.125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.125, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.125, 0.625));

  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(1.000, 0.500), r, back_it);
  ASSERT_EQ(neighbors.size(), 9);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.875, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.875, 0.625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.875, 0.875));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.125, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.375, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.125, 0.625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.375, 0.625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.125, 0.875));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.375, 0.875));
}

TEST(NEST_NEIGHBOR, unitmap_neighbors_3d_boundary_2cell) {
  typedef NEST<3> NestType;
  typedef NestType::ValueType VAL;
  NestType nest(2);
  std::vector<size_t> neighbors;
  std::back_insert_iterator<std::vector<size_t>> back_it(neighbors);
  size_t r = 4;
  int i;

  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(0, 0, 0), r, back_it);
  ASSERT_EQ(neighbors.size(), 8);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.03125, 0.03125, 0.03125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.09375, 0.03125, 0.03125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.03125, 0.09375, 0.03125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.09375, 0.09375, 0.03125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.03125, 0.03125, 0.09375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.09375, 0.03125, 0.09375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.03125, 0.09375, 0.09375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.09375, 0.09375, 0.09375));

  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(0, 1, 0), r, back_it);
  ASSERT_EQ(neighbors.size(), 4);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.03125, 0.96875, 0.03125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.09375, 0.96875, 0.03125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.03125, 0.96875, 0.09375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.09375, 0.96875, 0.09375));
  // cout << "size: " << neighbors.size() << endl;
  // BOOST_FOREACH(size_t i,neighbors) cout << nest.set_and_get(i,r) << " " << i
  // << endl;

  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(1, 0, 0), r, back_it);
  ASSERT_EQ(neighbors.size(), 12);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.96875, 0.03125, 0.03125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.96875, 0.09375, 0.03125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.96875, 0.03125, 0.09375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(0.96875, 0.09375, 0.09375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(1.03125, 0.03125, 0.03125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(1.09375, 0.03125, 0.03125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(1.03125, 0.09375, 0.03125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(1.09375, 0.09375, 0.03125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(1.03125, 0.03125, 0.09375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(1.09375, 0.03125, 0.09375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(1.03125, 0.09375, 0.09375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r),
            VAL(1.09375, 0.09375, 0.09375));
  // cout << "size: " << neighbors.size() << endl;
  // BOOST_FOREACH(size_t i,neighbors) cout << nest.set_and_get(i,r) << " " << i
  // << endl;
}

TEST(NEST_NEIGHBOR, unitmap_neighbors_2d_boundary_20cell) {
  typedef NEST<2> NestType;
  typedef NestType::ValueType VAL;
  NestType nest(20);
  std::vector<size_t> neighbors;
  std::back_insert_iterator<std::vector<size_t>> back_it(neighbors);
  size_t r = 2;
  int i;

  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(10.999, 0.499), r, back_it);
  ASSERT_EQ(neighbors.size(), 9);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(10.625, 0.125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(10.875, 0.125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(10.625, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(10.875, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(10.625, 0.625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(10.875, 0.625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(11.125, 0.125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(11.125, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(11.125, 0.625));

  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(11.000, 0.500), r, back_it);
  ASSERT_EQ(neighbors.size(), 9);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(10.875, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(10.875, 0.625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(10.875, 0.875));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(11.125, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(11.375, 0.375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(11.125, 0.625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(11.375, 0.625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(11.125, 0.875));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(11.375, 0.875));
}

TEST(NEST_NEIGHBOR, scalemap_neighbors_2d_boundary) {
  typedef NEST<2, util::SimpleArray<2, double>, ScaleMap> NestType;
  typedef NestType::ValueType VAL;
  NestType nest(NestType::Params(0, 0), NestType::Params(4, 4),
                NestType::Indices(4, 4));

  std::vector<size_t> neighbors;
  std::back_insert_iterator<std::vector<size_t>> back_it(neighbors);
  size_t r = 3;
  int i;

  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(0.999, 0.499), r, back_it);
  ASSERT_EQ(neighbors.size(), 9);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.8125, 0.3125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.9375, 0.3125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.8125, 0.4375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.9375, 0.4375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.8125, 0.5625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.9375, 0.5625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.0625, 0.3125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.0625, 0.4375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.0625, 0.5625));
  // cout << "size: " << neighbors.size() << endl;
  // BOOST_FOREACH(size_t i,neighbors) cout << nest.set_and_get(i,r) << " " << i
  // << endl;

  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(0.5, 0.999), r, back_it);
  ASSERT_EQ(neighbors.size(), 9);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.4375, 0.8125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.5625, 0.8125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.6875, 0.8125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.4375, 0.9375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.5625, 0.9375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.6875, 0.9375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.4375, 1.0625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.5625, 1.0625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.6875, 1.0625));
  // cout << "size: " << neighbors.size() << endl;
  // BOOST_FOREACH(size_t i,neighbors) cout << nest.set_and_get(i,r) << " " << i
  // << endl;

  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(1.0, 0.999), r, back_it);
  ASSERT_EQ(neighbors.size(), 9);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.9375, 0.8125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.9375, 0.9375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.0625, 0.8125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.1875, 0.8125));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.0625, 0.9375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.1875, 0.9375));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.9375, 1.0625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.0625, 1.0625));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.1875, 1.0625));
  // cout << "size: " << neighbors.size() << endl;
  // BOOST_FOREACH(size_t i,neighbors) cout << nest.set_and_get(i,r) << " " << i
  // << endl;

  r = 0;
  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(1.0, 0.999), r, back_it);
  ASSERT_EQ(neighbors.size(), 6);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.5, 0.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.5, 0.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(2.5, 0.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.5, 1.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.5, 1.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(2.5, 1.5));
  // cout << "size: " << neighbors.size() << endl;
  // BOOST_FOREACH(size_t i,neighbors) cout << nest.set_and_get(i,r) << " " << i
  // << endl;

  r = 0;
  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(1.0, 1.0), r, back_it);
  ASSERT_EQ(neighbors.size(), 9);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.5, 0.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.5, 0.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(2.5, 0.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.5, 1.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.5, 1.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(2.5, 1.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(0.5, 2.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(1.5, 2.5));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(2.5, 2.5));

  r = 1;
  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(3.5, 2.9), r, back_it);
  ASSERT_EQ(neighbors.size(), 6);
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(3.25, 2.25));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(3.75, 2.25));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(3.25, 2.75));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(3.75, 2.75));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(3.25, 3.25));
  EXPECT_EQ(nest.set_and_get(neighbors[i++], r), VAL(3.75, 3.25));
  // cout << "size: " << neighbors.size() << endl;
  // BOOST_FOREACH(size_t i,neighbors) cout << nest.set_and_get(i,r) << " " << i
  // << endl;

  r = 1;
  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(4.5, 2.9), r, back_it);
  ASSERT_EQ(neighbors.size(), 0);

  r = 1;
  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(4.4999, 2.9), r, back_it);
  ASSERT_EQ(neighbors.size(), 3);

  r = 3;
  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(4.124999, 2.9), r, back_it);
  ASSERT_EQ(neighbors.size(), 3);

  r = 3;
  neighbors.clear();
  i = 0;
  nest.get_neighbors(NestType::ValueType(4.125, 2.9), r, back_it);
  ASSERT_EQ(neighbors.size(), 0);
}

/*
template<class NestType>
void generic_test_neighbors(
        NestType nest,
        typename NestType::Params lb,
        typename NestType::Params ub
){
        typedef typename NestType::IndexType Index;
        typedef std::vector<typename NestType::IndexType> IndexVec;
        typedef std::set<typename NestType::IndexType> IndexSet;
        std::mt19937 rng((unsigned int)time(0));
        std::uniform_real_distribution<> uniform;
        typename NestType::ValueType randpt;
        for(size_t r = 0; r <= 9; ++r){
                // typename NestType::FloatType maxdis = 0;
                for(size_t iter=0; iter < 10/NestType::DIMENSION; ++iter){
                        for(size_t i = 0; i < NestType::DIMENSION; ++i){
                                randpt[i] = (ub[i]-lb[i])*uniform(rng)+lb[i];
                        }
                        IndexSet nbset;
                        Index center_index = nest.get_index(randpt,r);
                        nest.get_neighbors( randpt, r,
std::inserter(nbset,nbset.end()) );
                        BOOST_FOREACH(Index i,nbset){
                                cout << r << " " << center_index << " " << i <<
endl;
                        }
                        cout << endl;

                        // double dist = (randpt-nest.value()).norm();
                        // // cout << randpt.transpose() << " " <<
nest.value().transpose() << endl;
                        // maxdis = fmax(maxdis,dist);
                        // ASSERT_LE( dist, nest.bin_circumradius(r) );
                }
                // covering radius should be reasonably tigth
                // EXPECT_LT( nest.bin_circumradius(r)*0.75 , maxdis );
                // cout << DIM << " " << maxdis << " " <<
nest.bin_circumradius(r) << std::endl;
        }

}

template<int DIM>
void generic_test_neighbors_unit(){
        typedef NEST<DIM> NestType;
        typename NestType::Params lb,ub;
        lb.fill(0);
        ub.fill(1);
        generic_test_neighbors( NestType(), lb, ub );
}


TEST(NEST_NEIGHBOR,bin_inradius_unit){
        generic_test_neighbors_unit<2>();
}
*/
}
}
