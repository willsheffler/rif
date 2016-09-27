#ifndef INCLUDED_util_meta_print_type_HH
#define INCLUDED_util_meta_print_type_HH

#include "riflib/util/meta/util.hpp"

#include <iostream>
#include <boost/fusion/include/pair.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/is_sequence.hpp>
// #include <boost/mpl/copy.hpp>
// #include <boost/mpl/copy_if.hpp>
// #include <boost/mpl/set.hpp>
// #include <boost/mpl/insert.hpp>
// #include <boost/mpl/inserter.hpp>
#include <boost/mpl/pair.hpp>
// #include <boost/mpl/transform.hpp>
// // #include <boost/mpl/iter_fold.hpp>
#include <boost/mpl/for_each.hpp>
// #include <boost/mpl/vector.hpp>
// #include <boost/mpl/is_sequence.hpp>

// #include <boost/type_traits/remove_const.hpp>
// #include <boost/tuple/tuple.hpp>

namespace scheme {
namespace util {
namespace meta {

namespace m = boost::mpl;
namespace f = boost::fusion;

struct PrintType {
	std::ostream & out;
	std::string indent;
	bool newline;


	PrintType(
		std::ostream & _out=std::cout,
		std::string _indent="",
		bool n=true
	) : 
		out(_out),indent(_indent),newline(n) {}

	template <typename T>
	typename boost::disable_if_c<f::traits::is_sequence<T>::value>::type
	operator()( T const & ) const {
		out << indent << typeid(T).name();
		if(newline) out << std::endl; else out << " ";
	}
	template <typename T>
	typename boost::enable_if_c<f::traits::is_sequence<T>::value>::type
	operator()(T const &) const {
		out << indent << "TypeSeq<";
		if(newline) out << std::endl; else out << " ";		
		m::for_each<T>( PrintType(out,newline?indent+"    ":"",newline) );
		out << (newline?(indent+"  "):"") << ">";
		if(newline) out << std::endl; else out << " ";		
	}
	template <typename T>
	void operator()( type2type<T> const & = type2type<T>() ) const {
		out << indent << typeid(T).name();
		if(newline) out << std::endl; else out << " ";
	}
	template <typename A,typename B>
	void operator()(m::pair<A,B> const &) const {
		out << indent << "mpl::pair< " << typeid(A).name() << ", " << typeid(B).name() << " > ";
		if(newline) out << std::endl; else out << " ";
	}
	template <typename A,typename B>
	void operator()(std::pair<A,B> const &) const {
		out << indent << "std::pair< " << typeid(A).name() << ", " << typeid(B).name() << " > ";
		if(newline) out << std::endl; else out << " ";		
	}
	template <typename A,typename B>
	void operator()(f::pair<A,B> const &) const {
		out << indent << "fusion::pair< " << typeid(A).name() << ", " << typeid(B).name() << " > ";
		if(newline) out << std::endl; else out << " ";
	}



	// template <typename A,typename B>
	// void operator()(m::vector<A,B> const &) const {
	// 	PrintType p(out,indent,false);
	// 	out << indent << "mpl::vector< ";
	// 	p(A()); out << ", "; p(B()); out << " >";
	// 	if(newline) out << std::endl;
	// }
};

template<class T>
void print_type(
	std::ostream & out=std::cout,
	std::string indent="",
	bool newline=true)
{
	PrintType p(out,indent,newline);
	// this is OK because PrintType does nothing with the arg
	T *t(0);
	p.template operator()<T>(*t);
}

struct PrintInstanceType {
	std::ostream & out;
	std::string indent;
	PrintInstanceType(std::ostream & _out=std::cout,std::string _indent="") : out(_out),indent(_indent) {}
	template <typename T>
	void operator()(T const & x) const {
		out << indent << x << " (" << typeid(T).name() << ")" << std::endl;
	}
};

struct PrintBFMapofVec {
	std::ostream & out;
	std::string indent;
	PrintBFMapofVec(std::ostream & _out=std::cout,std::string _indent="") : out(_out),indent(_indent) {}
	template <typename T>
	void operator()(T const & x) const {
		out << indent << "KeyType: " << typeid(typename T::first_type).name() << std::endl;
		f::for_each( x.second, PrintInstanceType(out,indent+"    ") );
	}
};





}
}
}


namespace std {
	template<class A,class B>
	std::ostream & operator<<(std::ostream & out, std::pair<A,B> const & p){ 
		// return out << "std::pair< " << typeid(A).name() << ", " 
		//<< typeid(B).name() << " >( " << p.first << ", " << p.second << " )";
		return out << p.first << ", " << p.second;		
	}
}



#endif
