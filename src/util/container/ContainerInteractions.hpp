#ifndef INCLUDED_util_container_ContainerInteractions_HH
#define INCLUDED_util_container_ContainerInteractions_HH

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>
#include "types.hpp"
#include "util/StoragePolicy.hpp"
#include "util/meta/InstanceMap.hpp"
#include "util/meta/print_type.hpp"
#include "util/meta/util.hpp"

#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/vector.hpp>

#include <boost/iterator/iterator_facade.hpp>

namespace scheme {
namespace util {
namespace container {

////////////////// default get container interactions ////////////////

template <class Index>
struct ContainerInteractionsIter
    : boost::iterator_facade<
          ContainerInteractionsIter<Index>, std::pair<Index, Index> const,
          boost::forward_traversal_tag, std::pair<Index, Index> const> {
  typedef ContainerInteractionsIter<Index> THIS;
  ContainerInteractionsIter() : i1(0), i2(0), s2(0) {}
  ContainerInteractionsIter(Index ia, Index ib, Index sb)
      : i1(ia), i2(ib), s2(sb) {}

 private:
  friend class boost::iterator_core_access;
  void increment() {
    ++i2;
    if (i2 == s2) {
      i2 = 0;
      ++i1;
    }
  }
  std::pair<Index, Index> const dereference() const {
    return std::make_pair(i1, i2);
  }
  bool equal(THIS const &o) const {
    // std::cout << "ContainerInteractionsIter EQUAL: " << i1 << "-" << i2 <<
    // "-" << s1 << " / " << o.i1 << "-" << o.i2 << "-" << o.s1 << std::endl;
    assert(s2 == o.s2);
    return i1 == o.i1 && i2 == o.i2;
    // return reinterpret_cast<uint64_t const&>(*this) ==
    // reinterpret_cast<uint64_t const&>(o);
  }
  Index i1, i2, s2;
};

///@brief Default implementation of ContainerInteractions, double loop over all
/// by []/size()
template <class Xform, class Container1, class Container2, class Index>
struct ContainerInteractions {
  typedef std::pair<ContainerInteractionsIter<Index>,
                    ContainerInteractionsIter<Index>>
      Range;
  static void get_interaction_range(Xform const &,  // X * A = B, X = B * ~A
                                    Container1 const &c1, Container2 const &c2,
                                    Range &r) {
    // std::cout << "ContainerInteractions default" << std::endl;
    typedef ContainerInteractionsIter<Index> Iter;
    r = std::make_pair(
        Iter((Index)0, (Index)0, (Index)c2.size()),
        Iter((Index)(c2.size() ? c1.size() : 0), (Index)0, (Index)c2.size()));
  }
};
template <class T>
typename T::const_iterator get_cbegin(T const &t) {
  return t.begin();
}
template <class T>
typename T::const_iterator get_cend(T const &t) {
  return t.end();
}
template <class I>
I const &get_cbegin(std::pair<I, I> const &p) {
  return p.first;
}
template <class I>
I const &get_cend(std::pair<I, I> const &p) {
  return p.second;
}
template <class T>
struct get_citer {
  typedef typename T::const_iterator type;
};
template <class I>
struct get_citer<std::pair<I, I>> {
  typedef I type;
};
}
}
}

#endif
