#ifndef INCLUDED_util_meta_InstanceMap_HH
#define INCLUDED_util_meta_InstanceMap_HH

#include <boost/fusion/include/as_map.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/value_at_key.hpp>
#include <boost/mpl/always.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/zip_view.hpp>
#include <functional>
#include <vector>
#include "util/meta/util.hpp"

// #include <boost/serialization/access.hpp>
#ifdef CEREAL
#include <cereal/access.hpp>
#endif

namespace scheme {
namespace util {
namespace meta {

namespace m = boost::mpl;
namespace f = boost::fusion;

// struct FUSION_PAIRS { template<class T> struct apply { typedef void type; };
// };
struct FUSION_PAIRS {
  typedef void type;
};

namespace impl {

template <class Keys, class Arg2>
struct fusion_map_pairs {
  typedef typename m::eval_if<
      boost::is_same<FUSION_PAIRS, Arg2>,
      Keys,  // already fusion pairs
      m::transform<
          Keys,
          typename m::eval_if<
              // m::is_sequence<Arg2>,
              m::or_<boost::is_same<FUSION_PAIRS, Arg2>, m::is_sequence<Arg2> >,
              Arg2,  // It seems one of these two it getting instantiated even
              m::transform<Keys, Arg2>  // when m::eval_if<
                                        // boost::is_same<FUSION_PAIRS,Arg2> is
                                        // true ***
              >::type,
          f::pair<m::_1, m::_2> > >::type type;
};

///@brief convenience function to make a boost::fusion::map
///@detail fusion_map<Keys> will be same as tuple
///@detail fusion_map<Keys,Values> will map Key to Value
///@detail fusion_map<Keys,MetaFunc> will map Key to apply<MetaFunc,Key>::type
template <class Keys, class Arg2>
struct fusion_map {
  typedef typename f::result_of::as_map<
      typename fusion_map_pairs<Keys, Arg2>::type>::type type;
};

/// *** this compiles, so it seems like the upper eval_if should protect WTF?
///     typedef m::eval_if_c< true, m::identity<int>,
///     m::eval_if_c<true,void,void> >::type TEST;

template <class Imap, class Archive>
struct SerializeVisitor {
  Imap& fmap;
  Archive& archive;
  SerializeVisitor(Imap& f, Archive& ar) : fmap(f), archive(ar) {}
  template <class T>
  void operator()(type2type<T>) {
    archive& fmap.template get<T>();
  }
};

template <class Imap>
struct EqualsVisitor {
  bool is_equal;
  Imap const& imap1;
  Imap const& imap2;
  EqualsVisitor(Imap const& a, Imap const& b)
      : is_equal(true), imap1(a), imap2(b) {}
  template <class T>
  void operator()(type2type<T>) {
    is_equal &= imap1.template get<T>() == imap2.template get<T>();
  }
};
}

///@brief meta-container holding instances for any sequence of types
///@tparam _Keys sequence of Key types
///@tparam Arg2 sequence of Value types OR metafunction class OR placeholder
///expression
///@detail if Arg2 is a metafunction class, the values that func applied to the
///_Keys
template <typename _Keys, typename Arg2 = _Keys>
struct InstanceMap : impl::fusion_map<_Keys, Arg2>::type {
  typedef typename impl::fusion_map_pairs<_Keys, Arg2>::type Pairs;
  typedef InstanceMap<_Keys, Arg2> THIS;
  typedef typename impl::fusion_map<_Keys, Arg2>::type Base;
  typedef Base FusionType;
  // typedef typename m::eval_if< boost::is_same<FUSION_PAIRS,Arg2>,
  // m::identity<void>, _Keys >::type Keys;
  typedef
      typename m::transform<Pairs, util::meta::first_type<m::_1> >::type Keys;
  typedef typename m::transform<Pairs, util::meta::second_type<m::_1> >::type
      Values;

  // variadic ctors
  InstanceMap() {}
  template <class A>
  InstanceMap(A const& a) : Base(a) {}
  template <class A, class B>
  InstanceMap(A const& a, B const& b) : Base(a, b) {}
  template <class A, class B, class C>
  InstanceMap(A const& a, B const& b, C const& c) : Base(a, b, c) {}
  template <class A, class B, class C, class D>
  InstanceMap(A const& a, B const& b, C const& c, D const& d)
      : Base(a, b, c, d) {}
  template <class A, class B, class C, class D, class E>
  InstanceMap(A const& a, B const& b, C const& c, D const& d, E const& e)
      : Base(a, b, c, d, e) {}
  template <class A, class B, class C, class D, class E, class F>
  InstanceMap(A const& a, B const& b, C const& c, D const& d, E const& e,
              F const& f)
      : Base(a, b, c, d, e, f) {}
  template <class A, class B, class C, class D, class E, class F, class G>
  InstanceMap(A const& a, B const& b, C const& c, D const& d, E const& e,
              F const& f, G const& g)
      : Base(a, b, c, d, e, f, g) {}
  template <class A, class B, class C, class D, class E, class F, class G,
            class H>
  InstanceMap(A const& a, B const& b, C const& c, D const& d, E const& e,
              F const& f, G const& g, H const& h)
      : Base(a, b, c, d, e, f, g, h) {}
  template <class A, class B, class C, class D, class E, class F, class G,
            class H, class I>
  InstanceMap(A const& a, B const& b, C const& c, D const& d, E const& e,
              F const& f, G const& g, H const& h, I const& i)
      : Base(a, b, c, d, e, f, g, h, i) {}
  template <class A, class B, class C, class D, class E, class F, class G,
            class H, class I, class J>
  InstanceMap(A const& a, B const& b, C const& c, D const& d, E const& e,
              F const& f, G const& g, H const& h, I const& i, J const& j)
      : Base(a, b, c, d, e, f, g, h, i, j) {}

  ///@brief get reference to the instance of type associated with key Key
  ///@tparam Key input
  template <typename Key>
  typename f::result_of::value_at_key<Base, Key>::type& get() {
    return f::at_key<Key>((Base&)*this);
  }

  ///@brief get const reference to the instance of type associated with key Key
  ///@tparam Key input
  template <typename Key>
  typename f::result_of::value_at_key<Base, Key>::type const& get() const {
    return f::at_key<Key>((Base&)*this);
  }

  bool operator==(THIS const& that) const {
    impl::EqualsVisitor<THIS> eqv(*this, that);
    m::for_each<Keys, type2type<m::_1> >(eqv);
    return eqv.is_equal;
  }

// friend class boost::serialization::access;
#ifdef CEREAL
  friend class cereal::access;
#endif
  template <class Archive>
  void serialize(Archive& ar, const unsigned int) {
    impl::SerializeVisitor<THIS, Archive> serializer(*this, ar);
    m::for_each<Keys, type2type<m::_1> >(serializer);
  }
};

namespace impl {
template <class T>
struct EQUAL {
  T const& rhs;
  bool& is_equal;
  EQUAL(T const& r, bool& b) : rhs(r), is_equal(b) {}
  template <class X>
  void operator()(X const& x) const {
    is_equal &= (rhs.template get<typename X::first_type>() == x.second);
  }
};
template <class Float>
struct SUM {
  Float& sum;
  SUM(Float& s) : sum(s) {}
  template <class T>
  void operator()(T const& x) const {
    sum += x.second;
  }
};
template <class Float>
struct SETVAL {
  Float val;
  template <class T>
  void operator()(T& x) const {
    x.second = val;
  }
};
template <class T>
struct ADD {
  T& sink;
  ADD(T& s) : sink(s) {}
  template <class X>
  void operator()(X const& x) const {
    sink.template get<typename X::first_type>() += x.second;
  }
};
template <class T>
struct MUL {
  T& sink;
  MUL(T& s) : sink(s) {}
  template <class X>
  void operator()(X const& x) const {
    sink.template get<typename X::first_type>() *= x.second;
  }
};
template <class T, class OP>
struct BINARY_OP_EQUALS {
  T& sink;
  BINARY_OP_EQUALS(T& s) : sink(s) {}
  template <class X>
  void operator()(X const& x) const {
    sink.template get<typename X::first_type>() =
        OP()(sink.template get<typename X::first_type>(), x.second);
  }
};
template <class Float>
struct VEC {
  std::vector<Float>& vec;
  VEC(std::vector<Float>& s) : vec(s) {}
  template <class T>
  void operator()(T const& x) const {
    vec.push_back(x.second);
  }
};
}
template <class A, class B>
std::ostream& operator<<(std::ostream& out, InstanceMap<A, B> const& m) {
  return out << (typename InstanceMap<A, B>::FusionType&)m;
}

template <class T>
struct is_InstanceMap : m::false_ {};
template <class A, class B>
struct is_InstanceMap<InstanceMap<A, B> > : m::true_ {};

///@brief an InstanceMap where all value types are numeric, along with some
///element-wise binary ops
///@note values must be convertable into Float
template <typename Keys, typename Arg2 = Keys, class Float = double>
struct NumericInstanceMap : InstanceMap<Keys, Arg2> {
  typedef NumericInstanceMap<Keys, Arg2, Float> THIS;
  typedef InstanceMap<Keys, Arg2> BASE;
  typedef typename BASE::Base FusionType;
  typedef Float F;
  NumericInstanceMap(F f = 0) { setall(f); }
  NumericInstanceMap(F a, F b) : BASE(a, b) {}
  NumericInstanceMap(F a, F b, F c) : BASE(a, b, c) {}
  NumericInstanceMap(F a, F b, F c, F d) : BASE(a, b, c, d) {}
  NumericInstanceMap(F a, F b, F c, F d, F e) : BASE(a, b, c, d, e) {}
  NumericInstanceMap(F a, F b, F c, F d, F e, F f) : BASE(a, b, c, d, e, f) {}
  NumericInstanceMap(F a, F b, F c, F d, F e, F f, F g)
      : BASE(a, b, c, d, e, f, g) {}
  NumericInstanceMap(F a, F b, F c, F d, F e, F f, F g, F h)
      : BASE(a, b, c, d, e, f, g, h) {}
  NumericInstanceMap(F a, F b, F c, F d, F e, F f, F g, F h, F i)
      : BASE(a, b, c, d, e, f, g, h, i) {}
  NumericInstanceMap(F a, F b, F c, F d, F e, F f, F g, F h, F i, F j)
      : BASE(a, b, c, d, e, f, g, h, i, j) {}
  ///@briew set value of all instances
  void setall(Float val) {
    impl::SETVAL<Float> set;
    set.val = val;
    f::for_each((FusionType&)*this, set);
  }
  ///@briew sum of instance values
  Float sum() const {
    Float sum = 0;
    impl::SUM<Float> s(sum);
    f::for_each((FusionType&)*this, s);
    return sum;
  }
  template <class F2>
  void vector(std::vector<F2>& vec) const {
    impl::VEC<F2> s(vec);
    f::for_each((FusionType&)*this, s);
  }
  ///@brief test equality element by element
  bool operator==(THIS const& o) const {
    bool is_equal = true;
    impl::EQUAL<THIS> e(o, is_equal);
    f::for_each((FusionType&)*this, e);
    return is_equal;
    // return true;
  }
  bool operator!=(THIS const& o) const { return !(*this == o); }
  // /@brief convertable to Float as sum of elements
  // operator Float() const { return sum(); }
  void operator+=(THIS const& o) {
    impl::BINARY_OP_EQUALS<THIS, std::plus<Float> > add(*this);
    f::for_each((FusionType&)o, add);
  }
  void operator-=(THIS const& o) {
    impl::BINARY_OP_EQUALS<THIS, std::minus<Float> > add(*this);
    f::for_each((FusionType&)o, add);
  }
  void operator/=(THIS const& o) {
    impl::BINARY_OP_EQUALS<THIS, std::divides<Float> > add(*this);
    f::for_each((FusionType&)o, add);
  }
  void operator*=(THIS const& o) {
    impl::BINARY_OP_EQUALS<THIS, std::multiplies<Float> > add(*this);
    f::for_each((FusionType&)o, add);
  }
};
template <class A, class B, class C>
std::ostream& operator<<(std::ostream& out,
                         NumericInstanceMap<A, B, C> const& m) {
  return out << (typename NumericInstanceMap<A, B, C>::FusionType&)m;
}

template <class A, class B, class C, class O>
NumericInstanceMap<A, B, C> operator*(NumericInstanceMap<A, B, C> const& a,
                                      NumericInstanceMap<A, O, C> const& b) {
  NumericInstanceMap<A, B, C> result = a;
  result *= b;
  return result;
}
template <class A, class B, class C, class O>
NumericInstanceMap<A, B, C> operator+(NumericInstanceMap<A, B, C> const& a,
                                      NumericInstanceMap<A, O, C> const& b) {
  NumericInstanceMap<A, B, C> result = a;
  result += b;
  return result;
}
template <class A, class B, class C, class O>
NumericInstanceMap<A, B, C> operator-(NumericInstanceMap<A, B, C> const& a,
                                      NumericInstanceMap<A, O, C> const& b) {
  NumericInstanceMap<A, B, C> result = a;
  result -= b;
  return result;
}
template <class A, class B, class C, class O>
NumericInstanceMap<A, B, C> operator/(NumericInstanceMap<A, B, C> const& a,
                                      NumericInstanceMap<A, O, C> const& b) {
  NumericInstanceMap<A, B, C> result = a;
  result /= b;
  return result;
}

namespace impl {
SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(value_type, void)
}

struct std_vector_mfc {
  template <class T>
  struct apply {
    typedef std::vector<T> type;
  };
};

template <class DefaultCFC = std_vector_mfc>
struct make_container_pair {
  template <class T>
  struct apply {
    typedef typename impl::get_value_type_void<T>::type VALUE;
    typedef typename m::eval_if<
        boost::is_same<void, VALUE>,
        m::identity<f::pair<T, typename m::apply<DefaultCFC, T>::type> >,
        m::identity<f::pair<VALUE, T> > >::type type;
  };
};

template <class Containers>
struct ContainerInstanceMap
    : InstanceMap<
          typename m::transform<Containers, make_container_pair<> >::type,
          FUSION_PAIRS> {};

template <class A>
std::ostream& operator<<(std::ostream& out, ContainerInstanceMap<A> const& m) {
  return out << (typename ContainerInstanceMap<A>::FusionType&)m;
}
}
}
}

#endif
