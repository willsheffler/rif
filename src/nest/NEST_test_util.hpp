#ifndef INCLUDED_scheme_nest_NEST_test_util_HH
#define INCLUDED_scheme_nest_NEST_test_util_HH

#include <gtest/gtest.h>
#include "nest/NEST.hpp"

namespace scheme {
namespace nest {

template <typename Nest>
void generic_test_index_nesting_of_value(
    Nest nest, typename Nest::ValueType value,
    typename Nest::IndexType rmax = 12 / Nest::DIMENSION + 1) {
  typedef typename Nest::ValueType Value;
  typedef typename Nest::IndexType Index;
  int const DIM = Nest::DIMENSION;
  for (Index r = 0; r <= rmax; ++r) {
    Index index = nest.get_index(value, r);
    for (Index r2 = 0; r2 <= r; ++r2) {
      Index index2 = nest.get_index(value, r2);
      ASSERT_EQ(index >> (DIM * (r - r2)), index2);
    }
  }
}

template <typename Nest>
void generic_test_index_nesting_of_bincenters(
    Nest nest, typename Nest::IndexType rmax = 12 / Nest::DIMENSION + 1) {
  typedef typename Nest::ValueType Value;
  typedef typename Nest::IndexType Index;
  int const DIM = Nest::DIMENSION;
  for (Index r = 0; r <= rmax; ++r) {
    for (Index i = 0; i < nest.size(r); ++i) {
      ASSERT_TRUE(nest.set_state(i, r));
      Value value = nest.value();
      Index index = nest.get_index(value, r);
      ASSERT_EQ(i, index);
      for (Index r2 = 0; r2 <= r; ++r2) {
        Index index2 = nest.get_index(value, r2);
        ASSERT_EQ(i >> (DIM * (r - r2)), index2);
      }
    }
  }
}

template <typename Nest>
void generic_test_coverage_of_value(
    Nest nest, typename Nest::ValueType value,
    std::vector<double> &largest_d2_for_r,
    typename Nest::IndexType rmax = Nest::MAX_RESL_ONE_CELL - 2) {
  typedef typename Nest::ValueType Value;
  typedef typename Nest::IndexType Index;
  int const DIM = Nest::DIMENSION;
  for (Index r = 0; r <= rmax; ++r) {
    Index index = nest.get_index(value, r);
    if (index >= nest.size(r)) {
      std::cout << value << " " << index << " " << r << std::endl;
    }
    assert(index < nest.size(r));
    Value bincen = nest.set_and_get(index, r);
    double d2 = 0.0;
    for (size_t i = 0; i < DIM; ++i)
      d2 += (bincen[i] - value[i]) * (bincen[i] - value[i]);
    largest_d2_for_r[r] = fmax(largest_d2_for_r[r], d2);
    double const cr = nest.bin_circumradius(r);
    ASSERT_LE(d2, cr * cr);
  }
}
}
}

#endif
