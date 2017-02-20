#ifndef INCLUDED_util_meta_intersect_HH
#define INCLUDED_util_meta_intersect_HH

#include <boost/mpl/copy.hpp>
#include <boost/mpl/copy_if.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/set.hpp>
// #include <boost/mpl/inserter.hpp>

namespace scheme {
namespace util {
namespace meta {

namespace m = boost::mpl;

template <typename SeqA, typename SeqB>
struct intersect {
  typedef typename m::copy<
      SeqB, m::inserter<m::set<>, m::insert<m::_1, m::_2>>>::type set_seq2;
  typedef typename m::copy_if<SeqA, m::has_key<set_seq2, m::_1>>::type type;
};
}
}
}

#endif
