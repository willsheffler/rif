#include <gtest/gtest.h>

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/mpl/quote.hpp>
#include <boost/mpl/vector.hpp>
#include "util/meta/InstanceMap.hpp"
#include "util/meta/print_type.hpp"

#include "io/cache.hpp"
#include "util/SimpleArray.hpp"

namespace scheme {
namespace util {
namespace meta {

using std::cout;
using std::endl;

namespace bf = boost::fusion;
namespace mpl = boost::mpl;

namespace dummy {
struct T {};
struct U {};
struct V {};
}
TEST(NumericInstanceMap, basic_test) {
  using namespace dummy;
  typedef NumericInstanceMap<m::vector<T, U, V>, m::always<double>> NMAP;
  NMAP a(0, 0, 0);
  ASSERT_EQ(a.sum(), 0);

  ASSERT_EQ(NMAP(1, 2, 3), NMAP(1, 2, 3));
  ASSERT_NE(NMAP(1, 2, 3), NMAP(0, 2, 3));
  ASSERT_NE(NMAP(1, 2, 3), NMAP(1, 0, 3));
  ASSERT_NE(NMAP(1, 2, 3), NMAP(1, 2, 0));

  NMAP x;
  x.get<T>() = 1;
  x.get<U>() = 2;
  x.get<V>() = 3;
  ASSERT_EQ(NMAP(1, 2, 3), x);

  typedef NumericInstanceMap<m::vector<U, V, T>, m::always<double>> NMAP2;
  NMAP2 y;
  y.get<U>() = 1;
  y.get<V>() = 2;
  y.get<T>() = 3;
  ASSERT_EQ(NMAP2(1, 2, 3), y);
}

TEST(NumericInstanceMap, serialization) {
#ifdef CEREAL
  using namespace dummy;
  typedef NumericInstanceMap<m::vector<T, U, V>, m::always<double>> NMAP;
  NMAP imap(1, 2, 3);
  ASSERT_EQ(imap, io::test_serialization(imap));
#endif
}
}
}
}
