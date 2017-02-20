#ifndef INCLUDED_util_meta_util_HH
#define INCLUDED_util_meta_util_HH

#include <boost/mpl/bool.hpp>

namespace scheme {
namespace util {
namespace meta {

namespace m = boost::mpl;
// namespace f = boost::fusion;

template <class T>
struct type2type {};

template <class T>
struct showclass;
template <int T>
struct showint;

template <class T>
struct is_pair : m::false_ {};
template <class A, class B>
struct is_pair<std::pair<A, B>> : m::true_ {};

template <class T>
struct is_homo_pair : m::false_ {};
template <class A>
struct is_homo_pair<std::pair<A, A>> : m::true_ {};

template <class T>
struct first_type {
  typedef typename T::first_type type;
};
template <class T>
struct second_type {
  typedef typename T::second_type type;
};

// utility to get MEMBER if it exists, else double
struct __EMPTY_TYPE_UTILITY__ {};
template <typename T>
struct __TO_EMPTY_TYPE_UTILITY__ {
  typedef __EMPTY_TYPE_UTILITY__ type;
};

#define SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(MEMBER, DEFAULT)              \
  template <typename T,                                                   \
            typename Enable = scheme::util::meta::__EMPTY_TYPE_UTILITY__> \
  struct get_##MEMBER##_##DEFAULT {                                       \
    typedef DEFAULT type;                                                 \
  };                                                                      \
  template <typename T>                                                   \
  struct get_##MEMBER##_##DEFAULT<                                        \
      T, typename scheme::util::meta::__TO_EMPTY_TYPE_UTILITY__<          \
             typename T::MEMBER>::type> {                                 \
    typedef typename T::MEMBER type;                                      \
  };

#define SCHEME_MEMBER_TYPE_DEFAULT_SELF_TEMPLATE(MEMBER)                  \
  template <typename T,                                                   \
            typename Enable = scheme::util::meta::__EMPTY_TYPE_UTILITY__> \
  struct get_##MEMBER##_SELF {                                            \
    typedef T type;                                                       \
  };                                                                      \
  template <typename T>                                                   \
  struct get_##MEMBER##_SELF<                                             \
      T, typename scheme::util::meta::__TO_EMPTY_TYPE_UTILITY__<          \
             typename T::MEMBER>::type> {                                 \
    typedef typename T::MEMBER type;                                      \
  };

#define SCHEME_HAS_MEMBER_TYPE(MEMBER)                                    \
  template <typename T,                                                   \
            typename Enable = scheme::util::meta::__EMPTY_TYPE_UTILITY__> \
  struct has_type_##MEMBER : boost::mpl::false_ {};                       \
  template <typename T>                                                   \
  struct has_type_##MEMBER<                                               \
      T, typename scheme::util::meta::__TO_EMPTY_TYPE_UTILITY__<          \
             typename T::MEMBER>::type> : boost::mpl::true_ {};

template <typename T, class R, class A, class B, class C>
struct has_const_call_oper_3 {
  template <typename U, R (U::*)(A, B, C) const>
  struct SFINAE {};
  template <typename U>
  static char Test(SFINAE<U, &U::operator()> *);
  template <typename U>
  static int Test(...);
  static const bool value = sizeof(Test<T>(0)) == sizeof(char);
  typedef boost::mpl::bool_<value> type;
};

template <typename T, class R, class A>
struct has_subscript_oper {
  template <typename U, R (U::*)(A)>
  struct SFINAE {};
  template <typename U>
  static char Test(SFINAE < U, &U::operator[]> *);
  template <typename U>
  static int Test(...);
  static const bool value = sizeof(Test<T>(0)) == sizeof(char);
  typedef boost::mpl::bool_<value> type;
};

template <typename T, class R, class A>
struct has_const_subscript_oper {
  template <typename U, R (U::*)(A) const>
  struct SFINAE {};
  template <typename U>
  static char Test(SFINAE < U, &U::operator[]> *);
  template <typename U>
  static int Test(...);
  static const bool value = sizeof(Test<T>(0)) == sizeof(char);
  typedef boost::mpl::bool_<value> type;
};

#define SCHEME_HAS_MEMBER_FUNCTION_3(MEMBER)                      \
  template <typename T, class R, class A, class B, class C>       \
  struct has_member_fun_##MEMBER {                                \
    template <typename U, R (U::*)(A, B, C)>                      \
    struct SFINAE {};                                             \
    template <typename U>                                         \
    static char Test(SFINAE<U, &U::MEMBER> *);                    \
    template <typename U>                                         \
    static int Test(...);                                         \
    static const bool value = sizeof(Test<T>(0)) == sizeof(char); \
    typedef boost::mpl::bool_<value> type;                        \
  };
#define SCHEME_HAS_CONST_MEMBER_FUNCTION_0(MEMBER)                \
  template <typename T, class R>                                  \
  struct has_const_member_fun_##MEMBER {                          \
    template <typename U, R (U::*)() const>                       \
    struct SFINAE {};                                             \
    template <typename U>                                         \
    static char Test(SFINAE<U, &U::MEMBER> *);                    \
    template <typename U>                                         \
    static int Test(...);                                         \
    static const bool value = sizeof(Test<T>(0)) == sizeof(char); \
    typedef boost::mpl::bool_<value> type;                        \
  };
#define SCHEME_HAS_CONST_MEMBER_FUNCTION_2(MEMBER)                \
  template <typename T, class R, class A, class B>                \
  struct has_const_member_fun_##MEMBER {                          \
    template <typename U, R (U::*)(A, B) const>                   \
    struct SFINAE {};                                             \
    template <typename U>                                         \
    static char Test(SFINAE<U, &U::MEMBER> *);                    \
    template <typename U>                                         \
    static int Test(...);                                         \
    static const bool value = sizeof(Test<T>(0)) == sizeof(char); \
    typedef boost::mpl::bool_<value> type;                        \
  };
#define SCHEME_HAS_CONST_MEMBER_FUNCTION_3(MEMBER)                \
  template <typename T, class R, class A, class B, class C>       \
  struct has_const_member_fun_##MEMBER {                          \
    template <typename U, R (U::*)(A, B, C) const>                \
    struct SFINAE {};                                             \
    template <typename U>                                         \
    static char Test(SFINAE<U, &U::MEMBER> *);                    \
    template <typename U>                                         \
    static int Test(...);                                         \
    static const bool value = sizeof(Test<T>(0)) == sizeof(char); \
    typedef boost::mpl::bool_<value> type;                        \
  };
#define SCHEME_HAS_CONST_MEMBER_FUNCTION_4(MEMBER)                   \
  template <typename T, class R, class A, class B, class C, class D> \
  struct has_const_member_fun_##MEMBER {                             \
    template <typename U, R (U::*)(A, B, C, D) const>                \
    struct SFINAE {};                                                \
    template <typename U>                                            \
    static char Test(SFINAE<U, &U::MEMBER> *);                       \
    template <typename U>                                            \
    static int Test(...);                                            \
    static const bool value = sizeof(Test<T>(0)) == sizeof(char);    \
    typedef boost::mpl::bool_<value> type;                           \
  };

/////////////// in progress.... ////////////////

// template< typename SeqOfSeq >
// struct flatten {
// 	typedef void type;
// };

// template< typename S, template<typename,typename> class pair = std::pair >
// struct product1 {
// 	template< typename T >
// 	struct apply {
// 		typedef typename
// 			m::transform< S, pair<T,m::_1> >::type
// 		type;
// 	};
// };

// template< typename S, typename T=S, template<typename,typename> class pair =
// std::pair >
// struct product {
// 	typedef typename m::transform< S, product1<T,pair> >::type SeqOfSeq;
// 	typedef typename flatten< SeqOfSeq >::type type;
// };
}
}
}

#endif
