#include <gtest/gtest.h>

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/mpl/quote.hpp>
#include <boost/mpl/vector.hpp>
#include "util/meta/InstanceMap.hpp"
#include "util/meta/print_type.hpp"

#include "io/cache.hpp"
#include "util/SimpleArray.hpp"

namespace rif {
namespace util {
namespace meta {

using std::cout;
using std::endl;

namespace bf = boost::fusion;
namespace mpl = boost::mpl;

TEST(InstanceMap, fusion_map_test) {
  {
    typedef bf::map<bf::pair<int, char>, bf::pair<double, std::string>>
        map_type;

    map_type m(bf::make_pair<int>('X'), bf::make_pair<double>("Men"));

    ASSERT_EQ('X', bf::at_key<int>(m));
    ASSERT_EQ("Men", bf::at_key<double>(m));
  }
}

TEST(InstanceMap, holds_types) {
  using mpl::_1;
  using mpl::_2;

  typedef mpl::vector<int, char, float> Types;
  typedef mpl::vector<char, float, int> Types2;

  {
    // print_type< fusion_map_pairs<Types,_1>::type >();
    typedef InstanceMap<Types, _1> TEST;
    TEST imap;
    BOOST_STATIC_ASSERT((mpl::size<TEST>::value == mpl::size<Types>::value));
    BOOST_STATIC_ASSERT((bf::result_of::has_key<TEST, int>::value));
    ASSERT_TRUE(bf::has_key<int>(imap));
    BOOST_STATIC_ASSERT(
        (boost::is_same<int,
                        bf::result_of::value_at_key<TEST, int>::type>::value));
    BOOST_STATIC_ASSERT(
        (boost::is_same<char,
                        bf::result_of::value_at_key<TEST, char>::type>::value));
    imap.get<int>() = 1;
    imap.get<char>() = 'C';
    imap.get<float>() = 1.2345f;
    ASSERT_EQ(imap.get<int>(), 1);
    ASSERT_EQ(imap.get<char>(), 'C');
    ASSERT_EQ(imap.get<float>(), 1.2345f);
    // bf::for_each((TEST::Base&)imap,PrintInstanceType());
  }

  {
    InstanceMap<Types> imap;
    imap.get<int>() = 1;
    imap.get<char>() = 'C';
    imap.get<float>() = 1.2345f;
    ASSERT_EQ(imap.get<int>(), 1);
    ASSERT_EQ(imap.get<char>(), 'C');
    ASSERT_EQ(imap.get<float>(), 1.2345f);
  }
  {
    InstanceMap<Types, Types> imap;
    imap.get<int>() = 1;
    imap.get<char>() = 'C';
    imap.get<float>() = 1.2345f;
    ASSERT_EQ(imap.get<int>(), 1);
    ASSERT_EQ(imap.get<char>(), 'C');
    ASSERT_EQ(imap.get<float>(), 1.2345f);
  }
  {
    InstanceMap<Types, Types2> imap;
    imap.get<float>() = 1;
    imap.get<int>() = 'C';
    imap.get<char>() = 1.2345f;
    ASSERT_EQ(imap.get<float>(), 1);
    ASSERT_EQ(imap.get<int>(), 'C');
    ASSERT_EQ(imap.get<char>(), 1.2345f);
  }
  {
    InstanceMap<Types, _1> imap;
    imap.get<int>() = 1;
    imap.get<char>() = 'C';
    imap.get<float>() = 1.2345f;
    ASSERT_EQ(imap.get<int>(), 1);
    ASSERT_EQ(imap.get<char>(), 'C');
    ASSERT_EQ(imap.get<float>(), 1.2345f);
  }
  {
    InstanceMap<Types, std::vector<_1>> imap;
    imap.get<int>().push_back(1);
    imap.get<char>().push_back('C');
    imap.get<float>().push_back(1.2345f);
    ASSERT_EQ(imap.get<int>()[0], 1);
    ASSERT_EQ(imap.get<char>()[0], 'C');
    ASSERT_EQ(imap.get<float>()[0], 1.2345f);
    imap.get<float>().push_back(1.2345f);
    imap.get<float>().push_back(1.2345f);
    ASSERT_EQ(imap.get<float>().size(), 3);

    InstanceMap<Types, std::vector<_1>> imap2 = imap;
    ASSERT_TRUE(imap2 == imap);
  }
  bf::map<bf::pair<int, int>, bf::pair<char, char>, bf::pair<float, float>>
      test;
  // bf::for_each( test, PrintInstanceType() );
}

TEST(InstanceMap, can_use_fusion_pairs_directly) {
  // "usual" way
  InstanceMap<m::vector<int, char>, m::vector<char, float>> zip_imap;
  zip_imap.get<int>() = 'a';
  zip_imap.get<char>() = 1.234f;
  ASSERT_EQ(zip_imap.get<int>(), 'a');
  ASSERT_EQ(zip_imap.get<char>(), 1.234f);

  // make fusion map directly
  f::result_of::as_map<
      m::vector<f::pair<int, char>, f::pair<char, float>>>::type fmap;
  f::at_key<int>(fmap) = 'a';
  f::at_key<char>(fmap) = 1.234f;
  ASSERT_EQ(f::at_key<int>(fmap), 'a');
  ASSERT_EQ(f::at_key<char>(fmap), 1.234f);

  InstanceMap<m::vector<f::pair<int, char>, f::pair<char, float>>, FUSION_PAIRS>
      imap;
  imap.get<int>() = 'a';
  imap.get<char>() = 1.234f;
  ASSERT_EQ(imap.get<int>(), 'a');
  ASSERT_EQ(imap.get<char>(), 1.234f);
}

TEST(InstanceMap, serialization) {
#ifdef CEREAL
  InstanceMap<m::vector<int, char, float>> imap;
  imap.get<int>() = 1;
  imap.get<char>() = 'C';
  imap.get<float>() = 1.2345f;
  ASSERT_EQ(imap, io::test_serialization(imap));
  InstanceMap<m::vector<int, char, float>> const &cimap = imap;
  ASSERT_EQ(imap, io::test_serialization(cimap));
#endif
}

TEST(InstanceMap, subtyping) {
  typedef mpl::vector<int, char, float> Types;
  typedef mpl::vector<int> Types_int;

  {
    typedef InstanceMap<Types, mpl::_1> TEST;
    TEST imap;
    BOOST_STATIC_ASSERT((mpl::size<TEST>::value == mpl::size<Types>::value));
    BOOST_STATIC_ASSERT((bf::result_of::has_key<TEST, int>::value));
    ASSERT_TRUE(bf::has_key<int>(imap));
    BOOST_STATIC_ASSERT(
        (boost::is_same<int,
                        bf::result_of::value_at_key<TEST, int>::type>::value));
    BOOST_STATIC_ASSERT(
        (boost::is_same<char,
                        bf::result_of::value_at_key<TEST, char>::type>::value));
    imap.get<int>() = 1;
    imap.get<char>() = 'C';
    imap.get<float>() = 1.2345f;
    ASSERT_EQ(imap.get<int>(), 1);
    ASSERT_EQ(imap.get<char>(), 'C');
    ASSERT_EQ(imap.get<float>(), 1.2345f);
    // bf::for_each((TEST::Base&)imap,PrintInstanceType());

    //    typename TEST::Base & test1 = static_cast< typename TEST::Base &
    //>(
    // imap );

    // typename boost::fusion::detail::map_impl<0> & test2 = static_cast<
    // boost::fusion::detail::map_impl<0> & >( imap );
  }
}
}
}
}
