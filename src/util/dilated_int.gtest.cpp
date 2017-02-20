/// @brief  Z-order or Morton style indexing utilities for arbitrary dimension
/// @author will sheffler

// inspired by code from here: http://www.marcusbannerman.co.uk/dynamo
// see "Converting to and from Dilated Integers"(doi:10.1109/TC.2007.70814)

#include "util/dilated_int.hpp"
#include <gtest/gtest.h>

namespace scheme {
namespace util {

template <uint64_t D>
bool test_dilated_int() {
  uint64_t maxval = util::impl::MAX_DILATABLE<D>::VAL;
  maxval = std::min(maxval, uint64_t(2097151));
  for (uint64_t i(0); i < maxval; ++i) {
    uint64_t dilated = util::dilate<D>(i);
    uint64_t undilated = util::undilate<D>(dilated);
    EXPECT_EQ(undilated, i);
  }
  return true;
}

TEST(DILATED_INT, test_64bit) {
  util::test_dilated_int<1>();
  util::test_dilated_int<2>();
  util::test_dilated_int<3>();
  util::test_dilated_int<4>();
  util::test_dilated_int<5>();
  util::test_dilated_int<6>();
  util::test_dilated_int<7>();
  util::test_dilated_int<8>();
  util::test_dilated_int<9>();
  util::test_dilated_int<10>();
  util::test_dilated_int<11>();
  util::test_dilated_int<12>();
  util::test_dilated_int<13>();
  util::test_dilated_int<14>();
  util::test_dilated_int<15>();
  util::test_dilated_int<16>();
  util::test_dilated_int<17>();
  util::test_dilated_int<18>();
  util::test_dilated_int<19>();
  util::test_dilated_int<20>();
  util::test_dilated_int<21>();
  util::test_dilated_int<22>();
  util::test_dilated_int<23>();
  util::test_dilated_int<24>();
  util::test_dilated_int<25>();
  util::test_dilated_int<26>();
  util::test_dilated_int<27>();
  util::test_dilated_int<28>();
  util::test_dilated_int<29>();
  util::test_dilated_int<30>();
  util::test_dilated_int<31>();
  util::test_dilated_int<32>();
  // above 32 makes no sense for 64 bit
}
}
}
