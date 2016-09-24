#ifndef INCLUDED_util_meta_ref_wrap_HH
#define INCLUDED_util_meta_ref_wrap_HH

// #include <iostream>
// #include <boost/fusion/include/pair.hpp>
// #include <boost/fusion/include/for_each.hpp>
// #include <boost/fusion/include/mpl.hpp>
// #include <boost/fusion/include/is_sequence.hpp>
// #include <boost/mpl/copy.hpp>
// #include <boost/mpl/copy_if.hpp>
// #include <boost/mpl/set.hpp>
// #include <boost/mpl/insert.hpp>
// #include <boost/mpl/inserter.hpp>
// #include <boost/mpl/pair.hpp>
// #include <boost/mpl/transform.hpp>
// // #include <boost/mpl/iter_fold.hpp>
// #include <boost/mpl/for_each.hpp>
// #include <boost/mpl/vector.hpp>
// #include <boost/mpl/is_sequence.hpp>

// #include <boost/type_traits/remove_const.hpp>
#include <boost/tuple/tuple.hpp>

namespace scheme {
namespace util {
namespace meta {

namespace m = boost::mpl;
// namespace f = boost::fusion;


template<class T> struct remove_refwrap { typedef typename boost::remove_const<T>::type type; };
template<class T> struct remove_refwrap<boost::reference_wrapper<T> > { typedef typename boost::remove_const<T>::type type; };

///@brief works with pair and tuple
template<class T> struct recursive_remove_refwrap { typedef typename remove_refwrap<T>::type type; };
template<class A,class B>
	struct recursive_remove_refwrap<std::pair<A,B> > { typedef 
		std::pair< typename remove_refwrap<A>::type ,
		           typename remove_refwrap<B>::type > type; };
template<class A> 
	struct recursive_remove_refwrap<boost::tuple<A> > { typedef 
		boost::tuple< typename remove_refwrap<A>::type > type; };
template<class A,class B>
	struct recursive_remove_refwrap<boost::tuple<A,B> > { typedef 
		boost::tuple< typename remove_refwrap<A>::type ,
		              typename remove_refwrap<B>::type > type; };
template<class A,class B,class C>
	struct recursive_remove_refwrap<boost::tuple<A,B,C> > { typedef
		boost::tuple< typename remove_refwrap<A>::type ,
		              typename remove_refwrap<B>::type ,
		              typename remove_refwrap<C>::type > type; };
template<class A,class B,class C,class D>
	struct recursive_remove_refwrap<boost::tuple<A,B,C,D> > { typedef
		boost::tuple< typename remove_refwrap<A>::type ,
		              typename remove_refwrap<B>::type ,
		              typename remove_refwrap<C>::type ,
		              typename remove_refwrap<D>::type > type; };

// struct recursive_remove_refwrap_MFC {
// 	template<class T> struct apply { typedef typename recursive_remove_refwrap<T>::type type; }; };




}
}
}


#endif
