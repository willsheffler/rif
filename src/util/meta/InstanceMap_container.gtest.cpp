#include <gtest/gtest.h>

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/mpl/quote.hpp>
#include <boost/mpl/vector.hpp>
#include "util/meta/InstanceMap.hpp"
#include "util/meta/print_type.hpp"

#include "io/cache.hpp"
#include "util/SimpleArray.hpp"

#include <set>

namespace scheme {
namespace util {
namespace meta {

using std::cout;
using std::endl;

namespace bf = boost::fusion;
namespace mpl = boost::mpl;

TEST(ContainerInstanceMap, basic_test) {
  ContainerInstanceMap<m::vector<std::vector<int>, std::vector<char> > > cmap;
  BOOST_STATIC_ASSERT(
      (boost::is_same<std::vector<int>::value_type, int>::value));
  BOOST_STATIC_ASSERT(
      (boost::is_same<impl::get_value_type_void<int>::type, void>::value));
  cmap.get<int>().push_back(1);
  cmap.get<char>().push_back('c');
  cmap.get<char>().push_back('h');
  ASSERT_EQ(cmap.get<int>().at(0), 1);
  ASSERT_EQ(cmap.get<char>().at(0), 'c');
  ASSERT_EQ(cmap.get<char>().at(1), 'h');
}

TEST(ContainerInstanceMap, vector_default_test) {
  ContainerInstanceMap<m::vector<int,  // defaults to vector<int>
                                 char  // defaults to vector<char>
                                 > >
      cmap;
  cmap.get<int>().push_back(1);
  cmap.get<char>().push_back('c');
  cmap.get<char>().push_back('h');
  ASSERT_EQ(cmap.get<int>().at(0), 1);
  ASSERT_EQ(cmap.get<char>().at(0), 'c');
  ASSERT_EQ(cmap.get<char>().at(1), 'h');
}

TEST(ContainerInstanceMap, vector_default_set_test) {
  ContainerInstanceMap<m::vector<int,  // defaults to vector<int>
                                 std::set<char> > >
      cmap;
  cmap.get<int>().push_back(1);
  cmap.get<char>().insert('c');
  cmap.get<char>().insert('h');
  ASSERT_EQ(cmap.get<int>().at(0), 1);
  std::set<char>& tmp = cmap.get<char>();
  ASSERT_EQ(*tmp.find('c'), 'c');
  ASSERT_EQ(*tmp.find('h'), 'h');
  ASSERT_EQ(tmp.find('b'), tmp.end());
}

TEST(ContainerInstanceMap, simple_array_test) {
  ContainerInstanceMap<m::vector<SimpleArray<2, int>, SimpleArray<2, size_t> > >
      cmap;
  cmap.get<int>().at(0) = -1;
  cmap.get<size_t>()[0] = 1;
  cmap.get<size_t>()[1] = 2;
  ASSERT_EQ(cmap.get<int>().at(0), -1);
  SimpleArray<2, size_t>& tmp = cmap.get<size_t>();
  ASSERT_EQ(tmp.at(0), (size_t)1);
  ASSERT_EQ(tmp.at(1), (size_t)2);
}

TEST(ContainerInstanceMap, serialization_SimpleArray) {
#ifdef CEREAL
  ContainerInstanceMap<m::vector<SimpleArray<2, int>, SimpleArray<2, size_t> > >
      cmap;
  cmap.get<int>().at(0) = -1;
  cmap.get<int>().at(1) = -2;
  cmap.get<size_t>()[0] = 1;
  cmap.get<size_t>()[1] = 2;
  ASSERT_EQ(cmap, io::test_serialization(cmap));
#endif
}

TEST(ContainerInstanceMap, serialization_vector) {
#ifdef CEREAL
  ContainerInstanceMap<m::vector<std::vector<int>, std::vector<size_t> > > cmap;
  cmap.get<int>().push_back(-1);
  cmap.get<int>().push_back(6);
  cmap.get<int>().push_back(7);
  cmap.get<size_t>().push_back(3);
  cmap.get<size_t>().push_back(924);
  ASSERT_TRUE(cmap == io::test_serialization(cmap));
#endif
}
}
}
}
