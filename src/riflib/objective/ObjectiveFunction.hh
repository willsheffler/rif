#ifndef INCLUDED_objective_ObjectiveFunction_HH
#define INCLUDED_objective_ObjectiveFunction_HH

#include "riflib/util/meta/util.hh"
#include "riflib/util/meta/intersect.hh"
#include "riflib/util/meta/InstanceMap.hh"

#include <boost/foreach.hpp>

#include <boost/mpl/eval_if.hpp>
// #include <boost/mpl/set.hpp>
// #include <boost/mpl/inserter.hpp>
// #include <boost/mpl/insert.hpp>
#include <boost/mpl/for_each.hpp>
// #include <boost/mpl/copy.hpp>
// #include <boost/mpl/copy_if.hpp>
#include <boost/mpl/transform.hpp>
// #include <boost/mpl/int.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/mpl.hpp>

// #define DEBUG_IO

namespace scheme {
namespace objective {

namespace m = boost::mpl;
namespace f = boost::fusion;

namespace traits {
	// struct result_type { template<class T> struct apply { typedef typename T::Result type; }; };
	struct interaction_type { template<class T> struct apply { typedef typename T::Interaction type; }; };
}

namespace impl {

	// TODO: is there any way these could return references?
	template< class Interaction, class PlaceHolder, class InteractionSource, bool>
	struct get_interaction_from_placeholder_impl {
		Interaction operator()(PlaceHolder const & placeholder, InteractionSource const &){
			return placeholder;
		}
	};
	template< class Interaction, class PlaceHolder, class InteractionSource>
	struct get_interaction_from_placeholder_impl<Interaction,PlaceHolder,InteractionSource,false> {
		Interaction operator()(PlaceHolder const & placeholder, InteractionSource const & source){
			return source.template get_interaction_from_placeholder<Interaction>(placeholder);
		}
	};
	template< class Interaction, class PlaceHolder, class InteractionSource >
	Interaction get_interaction_from_placeholder(PlaceHolder const & placeholder, InteractionSource const & source){
		return get_interaction_from_placeholder_impl<
			Interaction,PlaceHolder,InteractionSource,boost::is_same<PlaceHolder,Interaction>::value>()
			(placeholder,source);
	}

	///@brief helper class tests Objective::Interaction == Interaction
	template<class Interaction>
	struct objective_interactions_equal {
		template<class Objective>
		struct apply : m::bool_< boost::is_same<Interaction,typename Objective::Interaction>::value > {};
	};

	///@brief helper class builds mpl::vector of Objectives with Objective::Interaction == Interaction
	template<class Objectives>
	struct copy_objective_if_interaction_equal {
		template<class Interaction>
		struct apply : m::copy_if<Objectives,objective_interactions_equal<Interaction>, m::back_inserter<f::vector<> > > {};
	};

	typedef struct{} NullScratch;
	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(Scratch,NullScratch)


	template< class Result, class Objective, class Interaction, class Scratch, class Config >
	void call_objective( Result & result, Objective const & objective, Interaction const & interaction,
		                 Scratch & scratch, Config const & config ){
		result = objective( interaction, scratch, config );
	}
	template< class Result, class Objective, class Interaction, class Config >
	void call_objective( Result & result, Objective const & objective, Interaction const & interaction,
		                 NullScratch &, Config const & config ){
		result = objective( interaction, config );
	}
	template< class Result, class Objective, class Actor1, class Actor2, class Scratch, class Config >
	void call_objective( Result & result, Objective const & objective, Actor1 const & actor1, Actor2 const & actor2,
		                 Scratch & scratch, Config const & config ){
		result = objective( actor1, actor2, scratch, config );
	}
	template< class Result, class Objective, class Actor1, class Actor2, class Config >
	void call_objective( Result & result, Objective const & objective, Actor1 const & actor1, Actor2 const & actor2,
		                 NullScratch &, Config const & config ){
		result = objective( actor1, actor2, config );
	}

	///@brief helper functor to call objective for all Interaction instances in InteractionSource
	///@detail this one calls objective( interaction, config ) this is more general than the split pair one because it can support
	/// onebody or more than twobody interactions, but can be less efficient for pair interactions
	//
	template<
		class Interaction,
		class Results,
		class Scratches,
		class Config
	>
	struct EvalObjective {
		Interaction const & interaction;
		Results & results;
		Scratches & scratches;
		Config const & config;
		double weight;
		EvalObjective(
			Interaction const & i,
			Results & r,
			Scratches & s,
			Config const & c,
			double w
		) : interaction(i),results(r),scratches(s),config(c),weight(w) {}

		template<class Objective>
		void
		operator()(Objective const & objective) const {
			BOOST_STATIC_ASSERT(( f::result_of::has_key<typename Results::FusionType,Objective>::value ));
			#ifdef DEBUG_IO
			std::cout << "    EvalObjective:     Objective " << Objective::name() <<"( " << interaction  << " )" << std::endl;
			#endif
			typedef typename f::result_of::value_at_key<Results,Objective>::type Result;
			/// TODO: unpack args to allow for tuples of std::ref, would avoid extra copy
			Result result;// = objective( interaction, config );
			call_objective( result, objective, interaction, scratches.template get<Objective>(), config );
			results.template get<Objective>() += weight*result;
		}
	};

	///@brief helper functor to call objective for all Interaction instances in InteractionSource
	//@detail this one calls objective( actor1, actor2, config )
	template<
		class Interaction,
		class Results,
		class Scratches,
		class Config
	>
	struct EvalObjectiveSplitPair {
		typename Interaction::first_type const & actor_1;
		typename Interaction::second_type const & actor_2;
		Results & results;
		Scratches & scratches;
		Config const & config;
		double weight;
		EvalObjectiveSplitPair(
			typename Interaction::first_type const & i,
			typename Interaction::second_type const & j,
			Results & r,
			Scratches & s,
			Config const & c,
			double w
		) : actor_1(i),actor_2(j),results(r),scratches(s),config(c),weight(w) {}

		template<class Objective>
		void
		operator()(Objective const & objective) const {
			BOOST_STATIC_ASSERT(( f::result_of::has_key<typename Results::FusionType,Objective>::value ));
			#ifdef DEBUG_IO
			std::cout << "    EvalObjectiveSplitPair:     Objective " << Objective::name() <<"( " << interaction  << " )" << std::endl;
			#endif
			typedef typename f::result_of::value_at_key<Results,Objective>::type Result;
			/// TODO: unpack args to allow for tuples of std::ref, would avoid extra copy
			Result result;// = objective( actor_1, actor_2, config );
			call_objective( result, objective, actor_1, actor_2, scratches.template get<Objective>(), config );
			results.template get<Objective>() += weight*result;
		}
	};

	using m::false_;
	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(DefinesInteractionWeight,false_)
	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(UseVisitor,false_)
	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(HasPre,false_)
	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(HasPost,false_)

	template<class InteractionSource, class Placeholder>
	typename boost::disable_if<typename get_DefinesInteractionWeight_false_<InteractionSource>::type,double>::type
	get_interaction_weight(InteractionSource const & , Placeholder const & ){
		return 1.0;
	}
	template<class InteractionSource, class Placeholder>
	typename boost::enable_if<typename get_DefinesInteractionWeight_false_<InteractionSource>::type,double>::type
	get_interaction_weight(InteractionSource const & is, Placeholder const & ph){
		return is.get_weight_from_placeholder(ph);
	}


	// call pre( Result &, Scratch &, Config const & ) iff objective has HasPre && Scratch is not NullScratch
	template< class InteractionSource, class Objective, class Scratch, class Config >
	typename boost::enable_if<typename get_HasPre_false_<Objective>::type>::type
	eval_objective_pre(
		InteractionSource const & source,
		Objective const & objective,
		typename Objective::Result & result,
		Scratch & scratch,
		Config const & config
	){
		objective.pre( source, result, scratch, config );
	}
	// call pre( Result &, Config const & ) iff objective has HasPre and Scratch is NullScratch
	template< class InteractionSource, class Objective, class Config >
	typename boost::enable_if<typename get_HasPre_false_<Objective>::type>::type
	eval_objective_pre(
		InteractionSource const & source,
		Objective const & objective,
		typename Objective::Result & result,
		NullScratch & scratch,
		Config const & config
	){
		objective.pre( source, result, config );
	}

	// call post( Result &, Config const & ) iff objective has HasPost
	template< class InteractionSource, class Objective, class Scratch, class Config >
	typename boost::enable_if<typename get_HasPost_false_<Objective>::type>::type
	eval_objective_post(
		InteractionSource const & source,
		Objective const & objective,
		typename Objective::Result & result,
		Scratch & scratch,
		Config const & config
	){
		objective.post( source, result, scratch, config );
	}
	// call post( Result &, Config const & ) iff objective has HasPost
	template< class InteractionSource, class Objective, class Config >
	typename boost::enable_if<typename get_HasPost_false_<Objective>::type>::type
	eval_objective_post(
		InteractionSource const & source,
		Objective const & objective,
		typename Objective::Result & result,
		NullScratch & scratch,
		Config const & config
	){
		objective.post( source, result, config );
	}

	// no pre/post
	template< class InteractionSource, class Objective, class Scratch, class Config >
	typename boost::disable_if<typename get_HasPre_false_<Objective>::type>::type
	eval_objective_pre(
		InteractionSource const & source,
		Objective const &,
		typename Objective::Result &,
		Scratch &,
		Config const &
	){
		; // appropriate pre function no available, do nothing
	}
	template< class InteractionSource, class Objective, class Scratch, class Config >
	typename boost::disable_if<typename get_HasPost_false_<Objective>::type>::type
	eval_objective_post(
		InteractionSource const &,
		Objective const &,
		typename Objective::Result &,
		Scratch &,
		Config const &
	){
		; // appropriate post function no available, do nothing
	}


	template< class InteractionSource, class Interaction, class Results, class Scratches, class Config >
	struct EvalObjectivePre {
		InteractionSource const & source_;
		Results & results;
		Scratches & scratches;
		Config const & config;
		EvalObjectivePre( InteractionSource const & src, Results & r, Scratches & s, Config const & c ) : source_(src),results(r),scratches(s),config(c) {}
		template<class Objective>
		void operator()(Objective const & objective) const {
			BOOST_STATIC_ASSERT(( f::result_of::has_key<typename Results::FusionType,Objective>::value ));
			#ifdef DEBUG_IO
			std::cout << "    EvalObjectivePost:     Objective " << Objective::name() << std::endl;
			#endif
			eval_objective_pre( source_, objective, results.template get<Objective>(), scratches.template get<Objective>(), config );
		}
	};

	template< class InteractionSource, class Interaction, class Results, class Scratches, class Config >
	struct EvalObjectivePost {
		InteractionSource const & source_;
		Results & results;
		Scratches & scratches;
		Config const & config;
		EvalObjectivePost( InteractionSource const & src, Results & r, Scratches & s, Config const & c ) : source_(src),results(r),scratches(s),config(c) {}
		template<class Objective>
		void operator()(Objective const & objective) const {
			BOOST_STATIC_ASSERT(( f::result_of::has_key<typename Results::FusionType,Objective>::value ));
			#ifdef DEBUG_IO
			std::cout << "    EvalObjectivePost:     Objective " << Objective::name() << std::endl;
			#endif
			eval_objective_post( source_, objective, results.template get<Objective>(), scratches.template get<Objective>(), config );
		}
	};

	template< class InteractionSource, class ObjectiveMap, class Results, class Scratches, class Config >
	struct EvalObjectivesPre {
		InteractionSource const & source_;
		ObjectiveMap const & objective_map_;
		Results & results_;
		Scratches & scratches_;
		Config const & config_;
		EvalObjectivesPre( InteractionSource const & s, ObjectiveMap const & f, Results & r, Scratches & scratches, Config const & c )
		 : source_(s), objective_map_(f), results_(r), scratches_(scratches), config_(c) {}
		template<class Interaction>
		void operator()(util::meta::type2type<Interaction>) const {
			#ifdef DEBUG_IO
				std::cout << "    EvalObjectivesPre Interaction: ";
				util::meta::PrintType(std::cout).operator()(Interaction());
			#endif
			f::for_each(
				objective_map_.template get<Interaction>(),
				EvalObjectivePre< InteractionSource, Interaction, Results, Scratches, Config >( source_, results_, scratches_, config_ )
			);
		}

	};
	template< class InteractionSource, class ObjectiveMap, class Results, class Scratches, class Config >
	struct EvalObjectivesPost {
		InteractionSource const & source_;
		ObjectiveMap const & objective_map_;
		Results & results_;
		Scratches & scratches_;
		Config const & config_;
		EvalObjectivesPost( InteractionSource const & s, ObjectiveMap const & f, Results & r, Scratches & scratches, Config const & c )
		 : source_(s), objective_map_(f),results_(r),scratches_(scratches),config_(c) {}

		template<class Interaction>
		void operator()(util::meta::type2type<Interaction>) const {
			#ifdef DEBUG_IO
				std::cout << "    EvalObjectivesPost Interaction: ";
				util::meta::PrintType(std::cout).operator()(Interaction());
			#endif
			f::for_each(
				objective_map_.template get<Interaction>(),
				EvalObjectivePost< InteractionSource, Interaction, Results, Scratches, Config >( source_, results_, scratches_, config_ )
			);
		}

	};




	///@brief helper functor to call each objective for Interaction
	template<
		class InteractionSource,
		class ObjectiveMap,
		class Results,
		class Scratches,
		class Config,
		class Enable = void
	>
	struct EvalObjectives {};

	// this EvalObjectives is for using iteration, in the case that the InteractionSource dosen't support visitors
	template<
		class InteractionSource,
		class ObjectiveMap,
		class Results,
		class Scratches,
		class Config
	>
	struct EvalObjectives<
		InteractionSource,
		ObjectiveMap,
		Results,
		Scratches,
		Config,
		typename boost::disable_if<typename get_UseVisitor_false_<InteractionSource>::type>::type
	> {
		InteractionSource const & interaction_source_;
		ObjectiveMap const & objective_map_;
		Results & results_;
		Scratches & scratches_;
		Config const & config_;

		EvalObjectives(
			InteractionSource const & p,
			ObjectiveMap const & f,
			Results & r,
			Scratches & s,
			Config const & c
		) : interaction_source_(p),objective_map_(f),results_(r),scratches_(s),config_(c) {
			// std::cout << "USE ITERATION" << std::endl;
		}

		template<class Interaction>	void
		operator()(util::meta::type2type<Interaction>) const {
			BOOST_STATIC_ASSERT(( InteractionSource::template has_interaction<Interaction>::value ));
			#ifdef DEBUG_IO
				std::cout << "    EvalObjectives Interaction: ";
				util::meta::PrintType(std::cout).operator()(Interaction());
			#endif
			typedef typename InteractionSource::template interaction_placeholder_type<Interaction>::type Placeholder;
			BOOST_FOREACH(
				Placeholder const & interaction_placeholder,
				interaction_source_.template get_interactions<Interaction>()
			){
				double weight = get_interaction_weight(interaction_source_,interaction_placeholder);
				// std::cout << "        Interaction: " << interaction <<
									 // " Result: " << typeid(results_.template get<Objective>()).name() << std::endl;
				Interaction const & interaction =
					get_interaction_from_placeholder<
							Interaction,
							Placeholder,
							InteractionSource
						>( interaction_placeholder, interaction_source_ );
				f::for_each(
					objective_map_.template get<Interaction>(),
					EvalObjective< Interaction, Results , Scratches, Config >
					             ( interaction, results_, scratches_, config_, weight )
				);
			}
		}

	};

	template<
		class _Interaction,
		class Objectives,
		class Results,
		class Scratches,
		class Config
	>
	struct ObjectivesVisitor {

		typedef _Interaction Interaction;
		typedef m::false_ Symmetric;

		Objectives const & objectives_;
		Results & results_;
		Scratches & scratches_;
		Config const & config_;

		ObjectivesVisitor(
			Objectives const & o,
			Results & r,
			Scratches & s,
			Config const & c
		) : objectives_(o), results_(r), scratches_(s), config_(c){}

		void
		operator()( Interaction const & interaction, double weight=1.0) {
			f::for_each(
				objectives_,
				EvalObjective< Interaction, Results , Scratches, Config >
					         ( interaction, results_, scratches_, config_, weight )
			);
		}
		#ifdef CXX11
			template< class I = Interaction >
		#else
			template< class I >
		#endif
		typename boost::enable_if< util::meta::is_pair<I> , void >::type
		operator()(
			typename I::first_type const & a1,
			typename I::second_type const & a2,
			double weight=1.0
		) {
			f::for_each(
				objectives_,
				EvalObjectiveSplitPair< I, Results, Scratches, Config >
					         ( a1, a2, results_, scratches_, config_, weight )
			);
		}

	};


	// this EvalObjectives is for using the Visitor method instead of iteration, faster but requires support in the InteractionSource object
	template<
		class InteractionSource,
		class ObjectiveMap,
		class Results,
		class Scratches,
		class Config
	>
	struct EvalObjectives<
		InteractionSource,
		ObjectiveMap,
		Results,
		Scratches,
		Config,
		typename boost::enable_if<typename get_UseVisitor_false_<InteractionSource>::type>::type
	> {
		InteractionSource const & interaction_source_;
		ObjectiveMap const & objective_map_;
		Results & results_;
		Scratches & scratches_;
		Config const & config_;

		EvalObjectives(
			InteractionSource const & p,
			ObjectiveMap const & f,
			Results & r,
			Scratches & s,
			Config const & c
		) : interaction_source_(p),objective_map_(f),results_(r),scratches_(s),config_(c) {
			// std::cout << "USE VISITATION" << std::endl;
		}

		template<class Interaction>	void
		operator()(util::meta::type2type<Interaction>) const {
			BOOST_STATIC_ASSERT(( InteractionSource::template has_interaction<Interaction>::value ));
			#ifdef DEBUG_IO
				std::cout << "    EvalObjectives Interaction: ";
				util::meta::PrintType(std::cout).operator()(Interaction());
			#endif
			typedef typename f::result_of::value_at_key<ObjectiveMap,Interaction>::type Objectives;
			Objectives const & objectives = objective_map_.template get<Interaction>();
			ObjectivesVisitor<Interaction,Objectives,Results,Scratches,Config> visitor(objectives,results_,scratches_,config_);
			interaction_source_.visit(visitor);
		}

	};



	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(InteractionTypes,void)
	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(Result,double)

	// get_interaction_weight

}

///@brief Minimal Concept for Objective used in ObectiveFunction
struct ObjectiveConcept {
	///@typedef required type of result this functor generates,
	///@note must be convertable to double but could have extra data
	typedef double Result;
	///@typedef required type of input evaluatable data, could be a pair of tuple
	typedef int Interaction;
	///@brief return name of Objective
	static std::string name(){ return "ObjectiveConcept"; }
	///@brief evaluate the Objective
	///@param Interaction main body or interaction to evaluate
	///@param Result result should be stored here
	template<class Config>
	Result operator()( Interaction const& a, Config const& ) const { return a; }
};

///@brief a generic objective function for interacting bodies
///@tparam Objectives sequence of types modeling the Objective concept
///@tparam Config a global config object passed to each Objective
template<
	typename _Objectives,
	typename _Config
	>
struct ObjectiveFunction {

	typedef _Objectives Objectives;
	typedef _Config Config;

	///@typedef only for error checking, make sure Objectives are non-redundant
	typedef typename
		m::copy<
			Objectives,
			m::inserter<m::set<>, m::insert<m::_1,m::_2> >
		>::type
	UniqueObjectives;

	BOOST_STATIC_ASSERT(( m::size<UniqueObjectives>::value == m::size<Objectives>::value ));

	///@typedef unique InteractionTypes types
	typedef typename
		m::copy< typename // copy back into mpl::vector;
			m::copy< typename // copy to mpl::set so unique
				m::transform< // all InteractionTypes Type from Obectives
					Objectives,
					traits::interaction_type
					>::type,
				m::inserter<m::set<>,m::insert<m::_1,m::_2> >
				>::type,
			m::back_inserter<m::vector<> >
		>::type
	InteractionTypes;

	///@typedef mpl::vector of Objectives for each unique PetalsType
	typedef typename
		m::transform<
			InteractionTypes,
			impl::copy_objective_if_interaction_equal<Objectives>,
			m::back_inserter<m::vector<> >
		>::type
	InteractionObjectives;

	///@typedef ObjectiveMap
	typedef util::meta::InstanceMap<
		InteractionTypes,
		InteractionObjectives
	> ObjectiveMap;

	///@typedef Results
	typedef util::meta::NumericInstanceMap<
		Objectives,
		impl::get_Result_double<m::_1>
	> Results;

	///@typedef Scratches
	typedef util::meta::InstanceMap<
		Objectives,
		impl::get_Scratch_NullScratch<m::_1>
	> Scratches;

	///@typedef Weights
	typedef util::meta::NumericInstanceMap<
		Objectives,
		impl::get_Result_double<m::_1> // TODO: fix this to avoid weights storing extra stuff???
	> Weights;

	ObjectiveMap objective_map_;
	Config default_config_;
	Weights weights_;

	///@brief default c'tor, init weights_ to 1
	ObjectiveFunction() : weights_(1.0) {}

	///@brief accessor for Objectives, may be used to configure or initialize Objective instances
	template<class Objective>
	Objective &
	get_objective(){
		BOOST_STATIC_ASSERT(( true )); // shold check that objective is actuall in Objectives...
		return f::deref( f::find<Objective>(
			objective_map_.template get<typename Objective::Interaction>()
		));
	}

	///@brief evaluate a InteractionSource with specified Config
	///@param InteractionSource must implement the InteractionSource concept
	///@param Config override default Config
	///@return a NumericInstanceMap of ResultTypes defined by Objectives
	template<class InteractionSource>
	void
	operator()(
		InteractionSource const & source,
		Config const & config,
		Results & results
	) const {
		// make sure we only operate on interactions contained in source
		typedef typename impl::get_InteractionTypes_void<InteractionSource>::type SourceInteractionTypes;
		typedef typename m::eval_if<
				boost::is_same<void,SourceInteractionTypes>,
				InteractionTypes,
				util::meta::intersect<
					InteractionTypes,
					SourceInteractionTypes
				>
			>::type
			MutualInteractionTypes;
		BOOST_STATIC_ASSERT(( m::size<MutualInteractionTypes>::value ));

		Scratches scratches;

		#ifdef DEBUG_IO
			std::cout << "ObjectiveFunction pre" << std::endl;
		#endif
		m::for_each< MutualInteractionTypes, util::meta::type2type<m::_1> >(
			impl::EvalObjectivesPre<
				InteractionSource,
				ObjectiveMap,
				Results,
				Scratches,
				Config
			>( source, objective_map_, results, scratches, config )
		);

		#ifdef DEBUG_IO
			std::cout << "ObjectiveFunction operator()" << std::endl;
		#endif
		m::for_each< MutualInteractionTypes, util::meta::type2type<m::_1> >(
			impl::EvalObjectives<
				InteractionSource,
				ObjectiveMap,
				Results,
				Scratches,
				Config
			>( source, objective_map_, results, scratches, config )
		);

		#ifdef DEBUG_IO
			std::cout << "ObjectiveFunction post" << std::endl;
		#endif
		m::for_each< MutualInteractionTypes, util::meta::type2type<m::_1> >(
			impl::EvalObjectivesPost<
				InteractionSource,
				ObjectiveMap,
				Results,
				Scratches,
				Config
			>( source, objective_map_, results, scratches, config )
		);

	}

	///@brief evaluate a InteractionSource with specified Config
	///@param InteractionSource must implement the InteractionSource concept
	///@param Config override default Config
	///@return a NumericInstanceMap of ResultTypes defined by Objectives
	template<class InteractionSource>
	Results
	operator()(
		InteractionSource const & source,
		Config const & config
	) const {
		Results results;
		this->template operator()<InteractionSource>(source,config,results);
		return results * weights_;
	}

	///@brief evaluate a InteractionSource with default Config
	///@param InteractionSource must implement the InteractionSource concept
	///@return a NumericInstanceMap of ResultTypes defined by Objectives
	template<class InteractionSource>
	Results
	operator()(
		InteractionSource const & source
	) const {
		return this->operator()<InteractionSource>(source,default_config_);
	}

};


}
}

#endif
