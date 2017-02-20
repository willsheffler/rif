#include <gtest/gtest.h>

#include <stdint.h>
#include "util/container/ContainerInteractions.hpp"

namespace scheme {
namespace util {
namespace container {

TEST(util_container_ContainerInteractions, default_on_vector) {
  typedef uint32_t Index;
  typedef ContainerInteractions<double, std::vector<int>, std::vector<int>,
                                Index>
      CI;
  typedef std::pair<Index, Index> Index2;

  BOOST_STATIC_ASSERT(
      (boost::is_same<std::pair<ContainerInteractionsIter<Index>,
                                ContainerInteractionsIter<Index> >,
                      CI::Range>::value));

  std::vector<int> a;
  std::vector<int> b;
  get_citer<CI::Range>::type beg;
  get_citer<CI::Range>::type end;
  CI::Range r;

  CI::get_interaction_range(0, a, b, r);
  beg = get_cbegin(r);
  end = get_cend(r);
  ASSERT_EQ(beg, end);

  a.push_back(0);  // 1,0

  CI::get_interaction_range(0, a, b, r);
  beg = get_cbegin(r);
  end = get_cend(r);
  ASSERT_EQ(beg, end);

  a.pop_back();
  b.push_back(0);  // 0,1

  CI::get_interaction_range(0, a, b, r);
  beg = get_cbegin(r);
  end = get_cend(r);
  ASSERT_EQ(beg, end);

  a.push_back(0);  // 1,1

  CI::get_interaction_range(0, a, b, r);
  beg = get_cbegin(r);
  end = get_cend(r);
  ASSERT_EQ(*beg++, Index2(0, 0));
  ASSERT_EQ(beg, end);

  a.push_back(0);
  b.push_back(0);
  b.push_back(0);

  CI::get_interaction_range(0, a, b, r);
  beg = get_cbegin(r);
  end = get_cend(r);
  ASSERT_EQ(*beg++, Index2(0, 0));
  ASSERT_EQ(*beg++, Index2(0, 1));
  ASSERT_EQ(*beg++, Index2(0, 2));
  ASSERT_EQ(*beg++, Index2(1, 0));
  ASSERT_EQ(*beg++, Index2(1, 1));
  ASSERT_EQ(*beg++, Index2(1, 2));
  ASSERT_EQ(beg, end);
}
}
}
}
namespace scheme {
namespace util {
namespace container {

template <class T>
struct myvector : std::vector<T> {
  myvector(size_t N, T def) : std::vector<T>(N, def) {}
};
template <class Xform, class A, class Container2, class Index>
struct ContainerInteractions<Xform, myvector<A>, Container2, Index> {
  typedef std::vector<std::pair<Index, Index> > Range;
  static void get_interaction_range(Xform const &, myvector<A> const &c1,
                                    Container2 const &c2, Range &r) {
    // std::cout << "ContainerInteractions specialization myvector10" <<
    // std::endl;
    for (Index i1 = 0; i1 < c1.size(); ++i1)
      for (Index i2 = 0; i2 < c2.size(); ++i2)
        r.push_back(std::make_pair(i1, i2));
  }
};
}
}
}

namespace scheme {
namespace util {
namespace container {

TEST(util_container_ContainerInteractions, specialization_on_myvector) {
  typedef uint16_t Index;
  typedef ContainerInteractions<double, myvector<int>, std::vector<int>, Index>
      CI;
  typedef std::pair<Index, Index> Index2;

  BOOST_STATIC_ASSERT((boost::is_same<std::vector<Index2>, CI::Range>::value));

  myvector<int> a(2, 0);
  std::vector<int> b(3, 1);

  CI::Range r;
  CI::get_interaction_range(0, a, b, r);
  // BOOST_FOREACH( Index2 p, r ) cout << p << endl;
  get_citer<CI::Range>::type beg = get_cbegin(r);
  get_citer<CI::Range>::type end = get_cend(r);
  ASSERT_EQ(*beg++, Index2(0, 0));
  ASSERT_EQ(*beg++, Index2(0, 1));
  ASSERT_EQ(*beg++, Index2(0, 2));
  ASSERT_EQ(*beg++, Index2(1, 0));
  ASSERT_EQ(*beg++, Index2(1, 1));
  ASSERT_EQ(*beg++, Index2(1, 2));
  ASSERT_EQ(beg, end);

  ASSERT_EQ(r.size(), 6);
}
}
}
}
