#include <gtest/gtest.h>

#include "riflib/objective/ObjectiveFunction.hh"
#include "riflib/objective/ObjectiveFunction_io.hh"

#include "riflib/util/meta/ref_wrap.hh"

#include <boost/mpl/vector.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/size.hpp>
#include <boost/static_assert.hpp>

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/io.hpp>

namespace scheme {
namespace objective {
namespace test {

namespace mpl = boost::mpl;
namespace bf = boost::fusion;
using std::cout;
using std::endl;


struct ScoreInt {
	typedef double Result;
	typedef int Interaction;
	double local_scale;
	ScoreInt():local_scale(1.0){}
	static std::string name(){ return "ScoreInt"; }
	template<class Config>
	Result operator()(Interaction a, Config const& c) const {
		return a*c.scale * local_scale;
	}
};
std::ostream & operator<<(std::ostream & out,ScoreInt const& si){ return out << si.name(); }

struct ScoreInt2 {
	typedef double Result;
	typedef int Interaction;
	static std::string name(){ return "ScoreInt2"; }
	template<class Config>
	Result operator()(Interaction const & a, Config const& c) const {
		return 2*a*c.scale;
	}
};
std::ostream & operator<<(std::ostream & out,ScoreInt2 const& si){ return out << si.name(); }
struct ScoreInt3 {
	typedef double Result;
	typedef int Interaction;
	static std::string name(){ return "ScoreInt3"; }
	template<class Config>
	Result operator()(Interaction const & a, Config const& c) const {
		return 3*a*c.scale;
	}
};
std::ostream & operator<<(std::ostream & out,ScoreInt3 const& si){ return out << si.name(); }

struct ScoreDouble {
	typedef double Result;
	typedef double Interaction;
	static std::string name(){ return "ScoreDouble"; }
	template<class Config>
	Result operator()(Interaction const & a, Config const& c) const {
		return a*c.scale;
	}
};
std::ostream & operator<<(std::ostream & out,ScoreDouble const& si){ return out << si.name(); }

struct ScoreIntDouble {
	typedef double Result;
	typedef std::pair<int,double> Interaction;
	static std::string name(){ return "ScoreIntDouble"; }
	template<class Config>
	Result operator()(int const & i,double const & a, Config const& c) const {
		return i*a*c.scale;
	}
	template<class Pair, class Config>
	Result operator()(Pair const & p, Config const& c) const {
		return this->operator()(p.first,p.second,c);
	}
};
std::ostream & operator<<(std::ostream & out,ScoreIntDouble const& si){ return out << si.name(); }

struct ConfigTest {
	double scale;
	ConfigTest():scale(1){}
};

TEST(ObjectiveFunction,detect_call_operator){
	BOOST_STATIC_ASSERT(( util::meta::has_const_call_oper_3<ScoreIntDouble,double,int const &,double const &,ConfigTest const &>::value ));
}

template<class Interactions >
struct SimpleInteractionSource {

	typedef util::meta::InstanceMap<Interactions,std::vector<mpl::_1> > MAP;

	MAP interactions_map_;

	typedef Interactions InteractionTypes;

	template<class Interaction>	std::vector<Interaction> &
	get_interactions() {
		return interactions_map_.template get<Interaction>();
	}
	template<class Interaction>	std::vector<Interaction> const &
	get_interactions() const {
		return interactions_map_.template get<Interaction>();
	}

	template<class Interaction>
	struct has_interaction : mpl::bool_<f::result_of::has_key<MAP,Interaction>::value> {};

	template<class Interaction>
	struct interaction_placeholder_type { typedef Interaction type; };

};
template<class I>
std::ostream & operator<<(std::ostream & out, SimpleInteractionSource<I> const & s){
	return out << s.interactions_map_;
}


TEST(ObjectiveFunction,basic_tests_local_and_global_config){
	typedef	ObjectiveFunction<
		mpl::list<
			ScoreInt,
			ScoreInt2,
			ScoreDouble,
			ScoreIntDouble
		>,
		ConfigTest
	> ObjFun;
	ObjFun score;

	typedef ObjFun::Results Results;
	typedef SimpleInteractionSource< mpl::vector<int,double,std::pair<int,double> > > InteractionSource;
	InteractionSource interaction_source;
	EXPECT_EQ( Results(0,0,0,0), score(interaction_source) );
	interaction_source.get_interactions<int>().push_back(1);
	EXPECT_EQ( Results(1,2,0,0), score(interaction_source) );
	interaction_source.get_interactions<double>().push_back(1.234);
	interaction_source.get_interactions<double>().push_back(2);
	interaction_source.get_interactions<double>().push_back(-2);
	EXPECT_EQ( Results(1,2,1.234,0), score(interaction_source) );
	interaction_source.get_interactions<std::pair<int,double> >().push_back(std::make_pair(1,1.5));
	EXPECT_EQ( Results(1,2,1.234,1.5), score(interaction_source) );

	score.default_config_.scale = 2.0;
	EXPECT_EQ( Results(2,4,2.468,3.0), score(interaction_source) );

	score.get_objective<ScoreInt>().local_scale = 2.0;
	EXPECT_EQ( Results(4,4,2.468,3.0), score(interaction_source) );

	std::ostringstream oss;
	oss << score << endl;
	// cout << oss.str();
}


TEST(ObjectiveFunction,test_results){
	typedef	ObjectiveFunction< mpl::vector< ScoreDouble, ScoreInt >, ConfigTest	> ObjFun;
	ObjFun score;
	// typedef util::meta::InstanceMap< mpl::vector<int,double>, std::vector<mpl::_1> > InteractionSource;
	typedef SimpleInteractionSource< mpl::vector<int,double> > InteractionSource;

	InteractionSource interaction_source;
	interaction_source.get_interactions<int>().push_back(1);
	interaction_source.get_interactions<int>().push_back(7);
	interaction_source.get_interactions<double>().push_back(1.2345);
	ObjFun::Results weights(2.0);
	ObjFun::Results results = score(interaction_source);
	float tot = results.sum();
	EXPECT_FLOAT_EQ( 9.2345f, tot );
	EXPECT_DOUBLE_EQ( 9.2345, results.sum() );
	EXPECT_DOUBLE_EQ( 18.469,  (results*weights).sum() );
	EXPECT_DOUBLE_EQ( 18.469,  (results+results).sum() );
	EXPECT_DOUBLE_EQ( 27.7035, (results*weights + results).sum() );
	EXPECT_DOUBLE_EQ( 9.2345,  (results*weights - results).sum() );
	EXPECT_DOUBLE_EQ( 9.2345,  (results*weights/weights).sum() );
}

template<class Interactions>
struct PlaceholderInteractionSource {
	typedef size_t Placeholder;
	typedef util::meta::InstanceMap<Interactions,std::vector<mpl::_1> > MAP;
	typedef util::meta::InstanceMap<Interactions,mpl::always<std::vector<Placeholder> > > MAP_int;

	MAP interactions_map_;
	MAP_int placeholder_map_;

	typedef Interactions InteractionTypes;

	template<class Interaction>
	void
	add_interaction(
		Interaction const & interaction
	){
		std::vector<Interaction> & vec = interactions_map_.template get<Interaction>();
		vec.push_back(interaction);
		placeholder_map_ .template get<Interaction>().push_back( vec.size()-1 );
		// cout << "add_interaction: " << interaction << " placeholder: " << vec.size()-1 << endl;
	}

	template<class Interaction>
	std::vector<Placeholder> const &
	get_interactions() const {
		return placeholder_map_.template get<Interaction>();
	}

	template<class Interaction>
	struct has_interaction : mpl::bool_<f::result_of::has_key<MAP,Interaction>::value> {};

	template<class Interaction>
	struct interaction_placeholder_type { typedef Placeholder type; };

	template<class Interaction>
	Interaction const &
	get_interaction_from_placeholder(
		Placeholder const & placeholder
	) const {
		// cout << "get_interaction_from_placeholder " << (interactions_map_.template get<Interaction>()[placeholder])
		//      << " placeholder " << placeholder << endl;
		return interactions_map_.template get<Interaction>()[placeholder];
	}

};
template<class I>
std::ostream operator<<(std::ostream & out, PlaceholderInteractionSource<I> const & s){
	return out << s.interactions_map_;
}


TEST(ObjectiveFunction,test_iteraction_placeholder){
	typedef	ObjectiveFunction<
		mpl::list<
			ScoreInt,
			ScoreInt2,
			ScoreInt3,
			ScoreDouble,
			ScoreIntDouble
		>,
		ConfigTest
	> ObjFun;
	ObjFun score;

	typedef ObjFun::Results Results;
	typedef PlaceholderInteractionSource< mpl::vector<int,double,std::pair<int,double> > > InteractionSource;
	InteractionSource interaction_source;
	EXPECT_EQ( Results(0,0,0,0,0), score(interaction_source) );
	interaction_source.add_interaction<int>(1);
	EXPECT_EQ( Results(1,2,3,0,0), score(interaction_source) );
	interaction_source.add_interaction<double>(1.234);
	interaction_source.add_interaction<double>(2);
	interaction_source.add_interaction<double>(-2);
	EXPECT_EQ( Results(1,2,3,1.234,0), score(interaction_source) );
	interaction_source.add_interaction<std::pair<int,double> >(std::make_pair(1,1.5));
	EXPECT_EQ( Results(1,2,3,1.234,1.5), score(interaction_source) );

	score.default_config_.scale = 2.0;
	EXPECT_EQ( Results(2,4,6,2.468,3.0), score(interaction_source) );

	score.get_objective<ScoreInt>().local_scale = 2.0;
	EXPECT_EQ( Results(4,4,6,2.468,3.0), score(interaction_source) );

}

TEST(ObjectiveFunction,ref_wrap_sanity){
	int i = 2;
	char c = 'b';
	std::pair<int,char> p(1,'a');
	std::pair<boost::reference_wrapper<int>,boost::reference_wrapper<char> > pr(boost::ref(i),boost::ref(c));
	EXPECT_NE( p.first, pr.first );	EXPECT_NE( p.second, pr.second );
	p = pr;
	EXPECT_EQ( p.first, pr.first );	EXPECT_EQ( p.second, pr.second );
	pr.first = boost::ref(i); pr.second = boost::ref(c);

	int & iref = pr.first;
	EXPECT_EQ(i,iref);
	iref = 4;
	EXPECT_EQ(i,4);

	// pr = p;
	// EXPECT_EQ( p.first, pr.first );

}


template<class IM_FROM, class IM_TO>
struct CopyIM {
	IM_FROM const & from_;
	IM_TO & to_;
	CopyIM(IM_FROM const & from,IM_TO & to) : from_(from),to_(to) {}
	template<class T>
	void operator()(util::meta::type2type<T>){
		typedef typename util::meta::recursive_remove_refwrap<T>::type Tnoref;
		// typedef std::pair<int,double> Tnoref;
		// cout << "FOO ";
		// util::meta::print_type<T>();
		// to_.template get<T>().resize(from_.template get<T>().size());
		// std::copy( from_.template get<T>().begin(), from_.template get<T>().end(), to_.template get<T>().begin() );
		// to_.template get<T>().clear();
		BOOST_FOREACH( Tnoref const & v, from_.template get<T>() ){
			to_.template get<T>().push_back( std::make_pair( boost::ref(v.first), boost::ref(v.second) ) );
		}
	}
};


struct vec_of_rem_refwrap {
	template<class T> struct apply { typedef
		std::vector<typename util::meta::recursive_remove_refwrap<T>::type> type; }; };

template<class Interactions >
struct RefInteractionSource {

	typedef util::meta::InstanceMap<Interactions,std::vector<mpl::_1> > MAPREF;
	typedef util::meta::InstanceMap<Interactions, vec_of_rem_refwrap > MAPVAL;

	MAPVAL imap_val_;
	MAPREF imap_ref_;

	// leave thes out here to test default behavoir
	// typedef Interactions InteractionTypes;

	void store_refs(){
		CopyIM<MAPVAL,MAPREF> cpim(imap_val_,imap_ref_);
		m::for_each<Interactions,util::meta::type2type<m::_1> >(cpim);
	}

	template<class Interaction>	std::vector<Interaction> const &
	get_interactions() const {
		return imap_ref_.template get<Interaction>();
	}

	template<class Interaction>	std::vector<typename util::meta::recursive_remove_refwrap<Interaction>::type > &
	get_interactions_vals() {
		return imap_val_.template get<Interaction>();
	}

	template<class Interaction>
	struct has_interaction : mpl::bool_<f::result_of::has_key<MAPREF,Interaction>::value> {};

	template<class Interaction>
	struct interaction_placeholder_type { typedef Interaction type; };

};
template<class I>
std::ostream & operator<<(std::ostream & out, RefInteractionSource<I> const & s){
	return out << s.imap_val_;
}

struct ScoreIntRefDoubleRef : ScoreIntDouble {
	typedef std::pair<boost::reference_wrapper<int const>, boost::reference_wrapper<double const> > Interaction;
	// template<class Config>
	// Result operator()(Interaction const & p, Config const& c) const {
	// 	int i = p.first;
	// 	double d = p.second;
	// 	this->template operator()<Config>(i,d,result,c);
	// }
};

TEST(ObjectiveFunction,tuple_of_ref_interactions){
	typedef	ObjectiveFunction<
		mpl::list<
			ScoreIntRefDoubleRef
		>,
		ConfigTest
	> ObjFun;
	ObjFun score;

	using boost::reference_wrapper;

	typedef ObjFun::Results Results;
	typedef std::pair<reference_wrapper<int const>,reference_wrapper<double const> > I1;
	typedef RefInteractionSource< mpl::vector<I1> > InteractionSource;
	InteractionSource interaction_source;

	EXPECT_EQ( Results(0), score(interaction_source) );
	interaction_source.get_interactions_vals<I1>().push_back(std::make_pair(1,1.5));
	interaction_source.store_refs();
	// std::vector< std::pair<int,double> > & tmp = interaction_source.get_interactions_vals<I1>();

	EXPECT_EQ( Results(1.5), score(interaction_source) );

	// std::ostringstream oss;
	// oss << score << endl;

}


struct ScoreIntWithPre {
	typedef double Result;
	typedef int Interaction;
	double local_scale;
	typedef m::true_ HasPre;
	ScoreIntWithPre():local_scale(1.0){}
	static std::string name(){ return "ScoreIntWithPre"; }
	template<class Config>
	Result operator()(Interaction a, Config const& c) const {
		return a*c.scale * local_scale;
	}
	template<class Scene, class Config>
	void pre( Scene const &, Result & r, Config const & c ) const {
		r += 1.0;
	}
};
std::ostream & operator<<(std::ostream & out,ScoreIntWithPre const& si){ return out << si.name(); }



TEST( ObjectiveFunction, test_pre )
{
	typedef	ObjectiveFunction<
		mpl::list<
			ScoreInt
		>,
		ConfigTest
	> ObjFunNoPre;
	ObjFunNoPre scorenopre;

	typedef	ObjectiveFunction<
		mpl::list<
			ScoreIntWithPre
		>,
		ConfigTest
	> ObjFunWithPre;
	ObjFunWithPre scorewithpre;

	typedef ObjFunNoPre::Results ResultsNoPre;
	typedef ObjFunWithPre::Results ResultsWithPre;
	typedef SimpleInteractionSource< mpl::vector<int,double,std::pair<int,double> > > InteractionSource;
	InteractionSource interaction_source;
	EXPECT_EQ( ResultsNoPre  (0), scorenopre  (interaction_source) );
	EXPECT_EQ( ResultsWithPre(1), scorewithpre(interaction_source) );
	interaction_source.get_interactions<int>().push_back(1);
	EXPECT_EQ( ResultsNoPre  (1), scorenopre  (interaction_source) );
	EXPECT_EQ( ResultsWithPre(2), scorewithpre(interaction_source) );

}

struct ScoreIntWithPost {
	typedef double Result;
	typedef int Interaction;
	double local_scale;
	typedef m::true_ HasPost;
	ScoreIntWithPost():local_scale(1.0){}
	static std::string name(){ return "ScoreIntWithPost"; }
	template<class Config>
	Result operator()(Interaction a, Config const& c) const {
		return a*c.scale * local_scale;
	}
	template<class Scene, class Config>
	void post( Scene const &, Result & r, Config const & c ) const {
		r += 1.0;
	}
};
std::ostream & operator<<(std::ostream & out,ScoreIntWithPost const& si){ return out << si.name(); }



TEST( ObjectiveFunction, test_post )
{
	typedef	ObjectiveFunction<
		mpl::list<
			ScoreInt
		>,
		ConfigTest
	> ObjFunNoPost;
	ObjFunNoPost scorenopost;

	typedef	ObjectiveFunction<
		mpl::list<
			ScoreIntWithPost
		>,
		ConfigTest
	> ObjFunWithPost;
	ObjFunWithPost scorewithpost;

	typedef ObjFunNoPost::Results ResultsNoPost;
	typedef ObjFunWithPost::Results ResultsWithPost;
	typedef SimpleInteractionSource< mpl::vector<int,double,std::pair<int,double> > > InteractionSource;
	InteractionSource interaction_source;
	EXPECT_EQ( ResultsNoPost  (0), scorenopost  (interaction_source) );
	EXPECT_EQ( ResultsWithPost(1), scorewithpost(interaction_source) );
	interaction_source.get_interactions<int>().push_back(1);
	EXPECT_EQ( ResultsNoPost  (1), scorenopost  (interaction_source) );
	EXPECT_EQ( ResultsWithPost(2), scorewithpost(interaction_source) );

}




struct ScoreIntWithScratch {
	typedef double Result;
	typedef int Interaction;
	double local_scale;
	typedef m::true_ HasPre;
	typedef m::true_ HasPost;
	typedef int Scratch;
	ScoreIntWithScratch():local_scale(1.0){}
	static std::string name(){ return "ScoreIntWithScratch"; }
	mutable Scratch * addr_of_scratch_should_stay_same;
	template<class Scene, class Config>
	void pre( Scene const &, Result & r, Scratch & s, Config const & c ) const {
		r += 1.0;
		addr_of_scratch_should_stay_same = &s;
	}
	template<class Config>
	Result operator()(Interaction a, Scratch & s, Config const& c ) const {
		if( addr_of_scratch_should_stay_same != &s ){
			std::cerr << "addr_of_scratch_should_stay_same" << std::endl;
			std::exit(-1);
		}
		return a*c.scale * local_scale;
	}
	template<class Scene, class Config>
	void post( Scene const &, Result & r, Scratch & s, Config const & c ) const {
		ASSERT_EQ( addr_of_scratch_should_stay_same, &s ); // verify that same scratch object is being passed
		r += 1.0;
	}
};
std::ostream & operator<<(std::ostream & out,ScoreIntWithScratch const& si){ return out << si.name(); }



TEST( ObjectiveFunction, test_scratch )
{
	typedef	ObjectiveFunction<
		mpl::list<
			ScoreInt
		>,
		ConfigTest
	> ObjFunNoScratch;
	ObjFunNoScratch scorenoscratch;

	typedef	ObjectiveFunction<
		mpl::list<
			ScoreIntWithScratch
		>,
		ConfigTest
	> ObjFunWithScratch;
	ObjFunWithScratch scorewithscratch;

	typedef ObjFunNoScratch::Results ResultsNoScratch;
	typedef ObjFunWithScratch::Results ResultsWithScratch;
	typedef SimpleInteractionSource< mpl::vector<int,double,std::pair<int,double> > > InteractionSource;
	InteractionSource interaction_source;
	EXPECT_EQ( ResultsNoScratch  (0), scorenoscratch  (interaction_source) );
	EXPECT_EQ( ResultsWithScratch(2), scorewithscratch(interaction_source) );
	interaction_source.get_interactions<int>().push_back(1);
	EXPECT_EQ( ResultsNoScratch  (1), scorenoscratch  (interaction_source) );
	EXPECT_EQ( ResultsWithScratch(3), scorewithscratch(interaction_source) );

}

}
}
}


