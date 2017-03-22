#include "util/SimpleArray.hpp"
#include <gtest/gtest.h>
#include <boost/foreach.hpp>

namespace rif {
namespace util {

using std::cout;
using std::endl;

TEST(SimpleArrayLegacy, bounds_check) {
  SimpleArrayLegacy<3, int> a;
  a[3];  // non-bounds checked
#ifndef NDEBUG
#ifndef CXX14
  ASSERT_DEATH(a.at(3), ".*");
#endif
#endif
}

TEST(SimpleArrayLegacy, iteration) {
  SimpleArrayLegacy<3, int> a;
  int v;
  v = 0;
  BOOST_FOREACH (int &i, std::make_pair(a.begin(), a.end()))
    i = ++v;
  v = 0;
  BOOST_FOREACH (int i, std::make_pair(a.begin(), a.end()))
    ASSERT_EQ(++v, i);
  v = 0;
  BOOST_FOREACH (int i, a)
    ASSERT_EQ(++v, i);
  SimpleArrayLegacy<3, int> const &r = a;
  v = 0;
  BOOST_FOREACH (int i, std::make_pair(r.begin(), r.end()))
    ASSERT_EQ(++v, i);
  v = 0;
  BOOST_FOREACH (int const &i, r)
    ASSERT_EQ(++v, i);
}
}
}
