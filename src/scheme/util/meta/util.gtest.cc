#include <gtest/gtest.h>

#include "scheme/util/meta/util.hh"
#include "scheme/util/meta/print_type.hh"
#include "scheme/util/meta/ref_wrap.hh"

#include <boost/tuple/tuple.hpp>
#include <boost/mpl/assert.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>

namespace scheme {
namespace util {
namespace meta {

// TEST(META,meta_product){
// 	typedef m::vector<char,int> Types;
// 	typedef m::vector<float,double> Types2;	
// 	m::for_each<Types>(PrintType());
// 	std::cout << "=========" << std::endl;
// 	m::for_each< product1<Types2>::apply<char>::type >(PrintType());
// 	m::for_each< product1<Types2>::apply<int >::type >(PrintType());
// 	std::cout << "=========" << std::endl;
// 	// typedef m::at< product<Types,Types2>::SeqOfSeq, m::int_<0> > PROD;
// 	// m::for_each< PROD >( PrintType() );

// }


using std::cout;
using std::endl;

struct EMPTY {};
struct DOUBLE { typedef double FOO; };
struct CHAR { typedef char FOO; };

struct VAL { typedef int value_type; };

SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(FOO,int)
SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(FOO,float)
SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(value_type,void)


SCHEME_MEMBER_TYPE_DEFAULT_SELF_TEMPLATE(FOO)


TEST(META_UTIL,member_type_default){

	BOOST_STATIC_ASSERT(( boost::is_same< int    , get_FOO_int<EMPTY  >::type >::value ));
	BOOST_STATIC_ASSERT(( boost::is_same< double , get_FOO_int<DOUBLE >::type >::value ));
	BOOST_STATIC_ASSERT(( boost::is_same< char   , get_FOO_int<CHAR   >::type >::value ));

	BOOST_STATIC_ASSERT(( boost::is_same< float  , get_FOO_float<EMPTY  >::type >::value ));

	BOOST_STATIC_ASSERT(( boost::is_same< EMPTY  , get_FOO_SELF<EMPTY  >::type >::value ));
	BOOST_STATIC_ASSERT(( boost::is_same< double , get_FOO_SELF<DOUBLE >::type >::value ));
	BOOST_STATIC_ASSERT(( boost::is_same< char   , get_FOO_SELF<CHAR   >::type >::value ));

	BOOST_STATIC_ASSERT(( boost::is_same< int    , get_FOO_SELF<int>::type >::value ));

	BOOST_STATIC_ASSERT(( boost::is_same< void   , get_value_type_void<EMPTY>::type >::value ));
	BOOST_STATIC_ASSERT(( boost::is_same< int    , get_value_type_void<VAL>::type >::value ));

}

using m::false_;
SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(TF,false_)
struct TRUE_  { typedef m::true_ TF; };
struct FALSE_ { typedef m::false_ TF; };

TEST(META_UTIL,member_type_default_mpl_TF){

	BOOST_STATIC_ASSERT(( !get_TF_false_<EMPTY >::type::value ));
	BOOST_STATIC_ASSERT((  get_TF_false_<TRUE_ >::type::value ));
	BOOST_STATIC_ASSERT(( !get_TF_false_<FALSE_>::type::value ));
	
}



SCHEME_HAS_MEMBER_TYPE(FOO)

TEST(META_UTIL,has_member_type){

	BOOST_STATIC_ASSERT(( !has_type_FOO<EMPTY  >::value ));
	BOOST_STATIC_ASSERT((  has_type_FOO<DOUBLE >::value ));
	BOOST_STATIC_ASSERT((  has_type_FOO<CHAR   >::value ));

}

TEST(META_UTIL,PrintType){
	m::vector<int,char> m;

	std::ostringstream out;
	PrintType(out).operator()(m);

	PrintType(out,"",false)(m); out << std::endl;

	// std::cout << out.str() << std::endl;
}

TEST(META_UTIL,remove_ref){
	using boost::is_same;
	using boost::tuple;
	using boost::reference_wrapper;
	using std::pair;

	BOOST_MPL_ASSERT(( is_same< int, remove_refwrap<int>::type > ));
	BOOST_MPL_ASSERT(( is_same< int, remove_refwrap<reference_wrapper<int> >::type > ));
	BOOST_MPL_ASSERT_NOT(( is_same< int, reference_wrapper<int> > ));

	BOOST_MPL_ASSERT(( is_same< int, remove_refwrap<const int>::type > ));
	BOOST_MPL_ASSERT(( is_same< int, remove_refwrap<reference_wrapper<const int> >::type > ));
	BOOST_MPL_ASSERT_NOT(( is_same< int const, reference_wrapper<int const> > ));

	BOOST_MPL_ASSERT(( is_same< int, recursive_remove_refwrap<int>::type > ));
	BOOST_MPL_ASSERT(( is_same< int, recursive_remove_refwrap<reference_wrapper<int> >::type > ));
	BOOST_MPL_ASSERT(( is_same< pair<int,char>, recursive_remove_refwrap<pair<int,char> >::type > ));
	BOOST_MPL_ASSERT(( is_same< pair<int,char>, recursive_remove_refwrap<pair<reference_wrapper<int>,char> >::type > ));
	BOOST_MPL_ASSERT(( is_same< pair<int,char>, recursive_remove_refwrap<
		pair<  reference_wrapper<int>,  reference_wrapper<char>   > >::type > ));

	BOOST_MPL_ASSERT(( is_same<
		tuple<                   int        >, recursive_remove_refwrap<
		tuple< reference_wrapper<int>       > >::type > ));
	BOOST_MPL_ASSERT(( is_same< 
		tuple<                   int ,                  char  >  ,   recursive_remove_refwrap<
		tuple< reference_wrapper<int>,reference_wrapper<char> >      >::type > ));
	BOOST_MPL_ASSERT(( is_same< 
		tuple<                   int ,                  char ,                   float  >  ,   recursive_remove_refwrap<
		tuple< reference_wrapper<int>,reference_wrapper<char>, reference_wrapper<float> >      >::type > ));
	BOOST_MPL_ASSERT(( is_same< 
		tuple<                   int ,                  char ,                   float ,                   int  >  ,   recursive_remove_refwrap<
		tuple< reference_wrapper<int>,reference_wrapper<char>, reference_wrapper<float>, reference_wrapper<int> >      >::type > ));

	BOOST_MPL_ASSERT(( is_same<
		tuple<                   int        >, recursive_remove_refwrap<
		tuple< reference_wrapper<int>       > >::type > ));
	BOOST_MPL_ASSERT(( is_same< 
		tuple<                   int ,                        char  >  ,   recursive_remove_refwrap<
		tuple< reference_wrapper<int>,reference_wrapper<const char> >      >::type > ));
	BOOST_MPL_ASSERT(( is_same< 
		tuple<                   int       ,                  char ,                   float  >  ,   recursive_remove_refwrap<
		tuple< reference_wrapper<int const>,reference_wrapper<char>, reference_wrapper<float> >      >::type > ));
	BOOST_MPL_ASSERT(( is_same< 
		tuple<                   int ,                  char ,                   float       ,                   int  >  ,   recursive_remove_refwrap<
		tuple< reference_wrapper<int>,reference_wrapper<char>, reference_wrapper<float const>, reference_wrapper<int> >      >::type > ));

	// print_type< recursive_remove_refwrap< tuple< reference_wrapper<int>,reference_wrapper<char> > >::type >();

}



struct USED_MEMORY_CONST {
	size_t used_memory(int,float,char) const { return 0; }
};
struct USED_MEMORY {
	size_t used_memory(int,float,char) { return 0; }
};

template<typename T, class A, class B, class C>
struct HasUsedMemoryMethod
{
    template<typename U, size_t (U::*)(A,B,C) const> struct SFINAE {};
    template<typename U> static char Test(SFINAE<U, &U::used_memory>*);
    template<typename U> static int Test(...);
    static const bool Has = sizeof(Test<T>(0)) == sizeof(char);
};

struct CALL_OPER_CONST {
	template<class T>
	void operator()(T,float,char) const {}
};

SCHEME_HAS_MEMBER_FUNCTION_3(used_memory)
SCHEME_HAS_CONST_MEMBER_FUNCTION_3(used_memory)

TEST(META_UTIL,detect_function){

	BOOST_STATIC_ASSERT(( !HasUsedMemoryMethod<EMPTY            ,int,float,char>::Has ));
	BOOST_STATIC_ASSERT((  HasUsedMemoryMethod<USED_MEMORY_CONST,int,float,char>::Has ));

	BOOST_STATIC_ASSERT(( !has_const_member_fun_used_memory<EMPTY            ,size_t,int,float,char>::value ));
	BOOST_STATIC_ASSERT(( !has_const_member_fun_used_memory<USED_MEMORY      ,size_t,int,float,char>::value ));
	BOOST_STATIC_ASSERT((  has_const_member_fun_used_memory<USED_MEMORY_CONST,size_t,int,float,char>::value ));	
	BOOST_STATIC_ASSERT(( !has_member_fun_used_memory      <EMPTY            ,size_t,int,float,char>::value ));
	BOOST_STATIC_ASSERT((  has_member_fun_used_memory      <USED_MEMORY      ,size_t,int,float,char>::value ));	
	BOOST_STATIC_ASSERT(( !has_member_fun_used_memory      <USED_MEMORY_CONST,size_t,int,float,char>::value ));	

	BOOST_STATIC_ASSERT(( !has_const_call_oper_3<EMPTY          ,void,int,float,char>::value ));	
	BOOST_STATIC_ASSERT((  has_const_call_oper_3<CALL_OPER_CONST,void,int,float,char>::value ));	

}

}
}
}

