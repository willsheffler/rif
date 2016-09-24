#ifndef INCLUDED_kinematics_Scene_HH
#define INCLUDED_kinematics_Scene_HH

#include "riflib/types.hh"
#include "riflib/util/assert.hh"
#include "riflib/kinematics/SceneBase.hh"
#include "riflib/kinematics/SceneIterator.hh"

#include "riflib/util/meta/InstanceMap.hh"
#include "riflib/util/container/ContainerInteractions.hh"

#include <boost/foreach.hpp>

#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/tuple/tuple.hpp>

#ifdef CEREAL
#include <cereal/access.hpp>
#include <cereal/types/memory.hpp>
#endif
// #include <boost/serialization/access.hpp>
// #include <boost/serialization/shared_ptr.hpp>

#include <vector>

namespace scheme {
namespace kinematics {

namespace m = boost::mpl;
namespace f = boost::fusion;

namespace impl {

	template< class ActorContainer >
	struct Conformation : util::meta::ContainerInstanceMap<ActorContainer> {
		typedef typename util::meta::ContainerInstanceMap<ActorContainer>::Keys Actors;
		Conformation(){}
		template<class Actor>
		void add_actor(Actor const & a){
			this->template get<Actor>().insert(this->template get<Actor>().end(),a);
		}
	};

	///@brief holds a Conformation pointer and an Xform
	template<
		class _Conformation,
		// class _Position,
		class Index
	>
	struct BodyTplt {
		typedef BodyTplt<_Conformation,Index> This;
		typedef _Conformation Conformation;
		typedef _Conformation ConformationConst; // TODO: fix this, requires figuring out Cereal load_and_construct
		// typedef _Position Position;

		BodyTplt() : conformation_() {}
		BodyTplt(shared_ptr<ConformationConst> c) : conformation_(c) {}

		BodyTplt( This const & other ){
			conformation_ = other.conformation_;
		}

		void deepcopy( This const & other ){
			conformation_ = make_shared<Conformation>( *other.conformation_ );
		}

		// BodyTplt( BodyTplt const & proto ){
		// 	conformation_ = make_shared<ConformationConst>( proto ); // NO!!! Confs are CONST
		// }

		// template<class Actor>
		// Actor
		// get_actor(
		// 	Index i
		// ) const {
		// 	return Actor(this->value()->template get<Actor>().at(i),position_);
		// }
		// template<class Actor>
		// Actor
		// get_actor_unsafe(Index i) const {
		// 	return Actor(this->value()->template get<Actor>()[i],position_);
		// }

		// Position const & position() const { return position_; }
		// void set_position(Position const & p){ position_ = p; }

		Conformation const & conformation() const { return *conformation_; }

		#ifdef CEREAL
		    // friend class boost::serialization::access;
		    friend class cereal::access;
	    	template<class Archive> void serialize(Archive & ar, const unsigned int ){
		    	// ar & position_;
		    	ar & conformation_;
		    }
	    #endif

	    bool operator==(This const & o) const { return *conformation_==*o.conformation_; }

	 private:
		// Position position_;
		shared_ptr<ConformationConst> conformation_;
	 };

	using m::true_;
	using m::false_;
	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(Symmetric,true_)
	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(RequireAbsolutePositioning,false_)

	template<class _Interaction>
	struct AccessVisitor {
		typedef _Interaction Interaction;

		Interaction i_;

		void operator()( Interaction const & i, double=0 ){ i_ = i; }

		#ifdef CXX11
			template< class I = Interaction >
		#else
			template< class I >
		#endif
		typename boost::enable_if< util::meta::is_pair<I> , void >::type
		operator()(
			typename I::first_type const & i,
			typename I::second_type const & j,
			double=0
		){ i_ = std::make_pair(i,j); }

		Interaction const & get() const { return i_; }
	};

	template<class Interaction>
	struct AccessVisitorAbsolute : AccessVisitor<Interaction> {
		typedef m::true_ RequireAbsolutePositioning;
	};

	template<class T,class I> struct get_placeholder_type {
		typedef typename std::pair<I,I> type; };
	template<class T1, class T2,class I> struct get_placeholder_type<std::pair<T1,T2>,I> {
		typedef typename std::pair<std::pair<I,I>,std::pair<I,I> > type; };


	template< class Scene >
	struct ActorAdder
	{
		Scene & scene_;
		boost::any const & actor_;
		typename Scene::Index ib_;
		bool & success_;

		ActorAdder(
			Scene & s,
			typename Scene::Index ib,
			boost::any const & a,
			bool & success
		) : scene_(s), ib_(ib), actor_(a), success_(success) {}

		template<class Actor>
		void operator()( Actor const & )
		{
		    try {
			    Actor tmp = boost::any_cast<Actor>( actor_ );
			    scene_.add_actor( ib_ , tmp );
		        success_ = true;
		    } catch(const boost::bad_any_cast &) {
		    	return;
		    }
		}
	};

	template< class Scene >
	struct ActorGetter
	{
		Scene const & scene_;
		boost::any const & actor_exemplar_;
		boost::any & actor_;
		typename Scene::Index ib_, ia_;
		bool & success_;
		ActorGetter(
			Scene const & s,
			typename Scene::Index ib,
			typename Scene::Index ia,
			boost::any const & aex,
			boost::any & a,
			bool & success
		) : scene_(s), ib_(ib), ia_(ia), actor_exemplar_(aex), actor_(a), success_(success) {}

		template<class Actor>
		void operator()( Actor const & )
		{
		    try {
			    boost::any_cast< util::meta::type2type<Actor> >( actor_exemplar_ );
			    actor_ = scene_.template get_actor<Actor>( ib_ , ia_ );
		        success_ = true;
		    } catch(const boost::bad_any_cast &) {
		    	return;
		    }
		}
	};

	template< class Scene >
	struct ActorGetterNonconst
	{
		Scene & scene_;
		boost::any const & actor_exemplar_;
		boost::any & actor_;
		typename Scene::Index ib_, ia_;
		bool & success_;
		ActorGetterNonconst(
			Scene & s,
			typename Scene::Index ib,
			typename Scene::Index ia,
			boost::any const & aex,
			boost::any & a,
			bool & success
		) : scene_(s), ib_(ib), ia_(ia), actor_exemplar_(aex), actor_(a), success_(success) {}

		template<class Actor>
		void operator()( Actor const & )
		{
		    try {
			    boost::any_cast< util::meta::type2type<Actor> >( actor_exemplar_ );
			    actor_ = &scene_.template get_nonconst_actor<Actor>( ib_ , ia_ );
		        success_ = true;
		    } catch(const boost::bad_any_cast &) {
		    	return;
		    }
		}
	};

	template< class Scene >
	struct ActorClearer
	{
		Scene & scene_;
		boost::any const & actor_exemplar_;
		typename Scene::Index ib_;
		bool & success_;
		ActorClearer(
			Scene & s,
			typename Scene::Index ib,
			boost::any const & aex,
			bool & success
		) : scene_(s), ib_(ib), actor_exemplar_(aex), success_(success) {}

		template<class Actor>
		void operator()( Actor const & )
		{
		    try {
			    boost::any_cast< Actor const & >( actor_exemplar_ );
			    scene_.template clear_actors<Actor>( ib_ );
		        success_ = true;
		    } catch(const boost::bad_any_cast &) {
		    	return;
		    }
		}
	};

	template< class Scene >
	struct ActorCounter
	{
		Scene const & scene_;
		boost::any const & actor_;
		typename Scene::Index ib_;
		int & count_;
		ActorCounter(
			Scene const & s,
			typename Scene::Index ib,
			boost::any const & a,
			int & count
		) : scene_(s), ib_(ib), actor_(a), count_(count) {}

		template<class Actor>
		void operator()( Actor const & )
		{
		    try {
			    boost::any_cast<util::meta::type2type<Actor> >( actor_ );
			    count_ = scene_.template num_actors<Actor>( ib_ );
		    } catch(const boost::bad_any_cast &) {
		    	return;
		    }
		}
	};


}


///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Scene / //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

	///@brief default implementation for inverse
	template<typename T> T inverse(T const & t){ return t.inverse(); }

	template<
		class _ActorContainers,
		class _Position,
		class _Index = uint64_t
	>
	struct Scene : public SceneBase<_Position,_Index> {

		typedef SceneBase<_Position,_Index> Base;
		typedef Scene<_ActorContainers,_Position,_Index> This;
		typedef _Position Position;
		typedef _Index Index;
		typedef impl::Conformation<_ActorContainers> Conformation;
		typedef impl::Conformation<_ActorContainers> ConformationConst;
		typedef typename Conformation::Actors Actors;
		typedef impl::BodyTplt<Conformation,Index> Body;
		typedef std::vector<Body> Bodies;
		typedef m::true_ DefinesInteractionWeight;
		typedef m::true_ UseVisitor;

		Bodies bodies_;

		Scene(Index nbodies=0) : SceneBase<Position,Index>() {
			for(Index i=0; i<nbodies; ++i) add_body();
			this->update_symmetry( (Index)bodies_.size() );
		}

		Scene( This const & o, bool deep = 1 ) : Base(o) {
			if( deep ){
				bodies_.resize(o.bodies_.size());
				for( int i = 0; i < bodies_.size(); ++i ){
					bodies_[i].deepcopy( o.bodies_[i] );
				}
			} else {
				bodies_ = o.bodies_; // is shallow copy
			}
		}
		virtual shared_ptr<Base> clone_deep() const {
			return make_shared<This>(*this,true);
		}
		virtual shared_ptr<Base> clone_shallow() const {
			return make_shared<This>(*this,false);
		}

		virtual ~Scene(){}


		// virtual Position position(Index i) const {
		// 	Index isym = this->sym_index_map(i);
		// 	// if( isym == 0 )	return                    bodies_.at(i).position();
		// 	/*else*/ return this->symframes_.at(isym) * bodies_.at(i).position();
		// }
		// virtual void set_position(Index i, Position const & newp){
		// 	assert( bodies_.size() == this->positions_.size() );
		// 	assert(i < bodies_.size());
		// 	this->set_position(newp);
		// }

		#ifdef CEREAL
		    // friend class boost::serialization::access;
		    friend class cereal::access;
	    	template<class Archive> void serialize(Archive & ar, const unsigned int ){
		    	ar & bodies_;
			    ar & this->symframes_;
		   		ar & this->n_sym_bodies_;
		   		ar & this->n_bodies_;
		   		ar & this->positions_;
	    	}
	    #endif
	    bool operator==(This const & o) const {
			assert( bodies_.size() == this->positions_.size() );
	    	return this->bodies_==o.bodies_ &&
	    	       this->symframes_==o.symframes_ &&
	    	       this->n_sym_bodies_==o.n_sym_bodies_ &&
	    	       this->positions_==o.positions_;
	    }


		/// mutators
		void add_body(){
			bodies_.push_back( make_shared<ConformationConst>() );
			this->positions_.push_back(Position::Identity());
			this->update_symmetry( (Index)bodies_.size() );
			assert( bodies_.size() == this->positions_.size() );
		}
		void add_body(shared_ptr<ConformationConst> const & c){
			bodies_.push_back(c);
			this->positions_.push_back(Position::Identity());
			this->update_symmetry( (Index)bodies_.size() );
			assert( bodies_.size() == this->positions_.size() );
		}
		// void add_body(Body const & b){
		// 	bodies_.push_back(b);
		// 	positions_.push_back(Position::Identity());
		// 	this->update_symmetry( (Index)bodies_.size() );
		// }
		// void set_body(Index i, shared_ptr<ConformationConst> const & c){ bodies_[i] = c; }
		// void set_body(Index i, Body const & b){ bodies_[i] = b; }
		Conformation & mutable_conformation_asym(Index i) { return const_cast<Conformation&>(bodies_.at(i).conformation()); }


		Conformation const & conformation(Index i) const { return bodies_.at(i%bodies_.size()).conformation(); }
		// Body const & body       (Index i) const { return bodies_.at(i); }

		/// mainly internal use
		Bodies const & __bodies_unsafe__() const { return bodies_; }
		std::vector<Position> & __symframes_unsafe__(){ return this->symframes_; }

		Position const & __position_unsafe_asym__(Index i) const { return this->positions_[i]; }
		Position       & __position_unsafe_asym__(Index i)       { return this->positions_[i]; }
		Position const __position_unsafe__(Index i) const {
			Index isym = this->sym_index_map(i);
			// if( isym == 0 )	return this->positions_[i];
			/*else*/ return this->symframes_[isym] * this->positions_[i];
		}
		Conformation const & __conformation_unsafe_asym__(Index i) const { return bodies_[i].conformation(); }
		Conformation const & __conformation_unsafe__(Index i) const { return bodies_[i%bodies_.size()].conformation(); }



		template<class Actor>
		typename boost::enable_if< impl::has_type_Position<Actor>, Actor >::type
		get_actor(Index ib, Index ia) const {
			return Actor( conformation(ib).template get<Actor>().at(ia), this->position(ib) );
		}
		template<class Actor>
		typename boost::disable_if< impl::has_type_Position<Actor>, Actor >::type
		get_actor(Index ib, Index ia) const {
			return conformation(ib).template get<Actor>().at(ia);
		}
		template<class Actor>
		void
		clear_actors(Index ib){
			mutable_conformation_asym(ib).template get<Actor>().clear();
		}


		template<class Actor>
		Actor &
		get_nonconst_actor( Index ib, Index ia ){
			return mutable_conformation_asym(ib).template get<Actor>().at(ia);
		}

		template<class Actor>
		void add_actor(Index ib, Actor const & a) {
			mutable_conformation_asym(ib).template add_actor<Actor>(a);
		}

		virtual bool add_actor( Index ib, boost::any const & a ){
			bool success = false;
			impl::ActorAdder<This> adder( *this, ib, a, success );
			boost::mpl::for_each<Actors>( adder );
			return success;
		}

		virtual bool get_actor( Index ib, Index ia, boost::any & a_in ) const {
			bool success = false;
			boost::any a;
			impl::ActorGetter<This> getter( *this, ib, ia, a_in, a, success );
			boost::mpl::for_each<Actors>( getter );
			if( success ) a_in = a;
			return success;
		}

		virtual bool get_nonconst_actor( Index ib, Index ia, boost::any & a_in ){
			bool success = false;
			boost::any a;
			impl::ActorGetterNonconst<This> getter( *this, ib, ia, a_in, a, success );
			boost::mpl::for_each<Actors>( getter );
			if( success ) a_in = a;
			return success;
		}

		virtual bool clear_actors( Index ib, boost::any const & a_in ) {
			bool success = false;
			impl::ActorClearer<This> clearer( *this, ib, a_in, success );
			boost::mpl::for_each<Actors>( clearer );
			return success;
		}

		template<class Actor>
		size_t num_actors(Index ib) const {
			return bodies_.at(ib).conformation().template get<Actor>().size();
		}

		virtual int num_actors( Index ib, boost::any const & a ) const {
			int num_actors = -1;
			impl::ActorCounter<This> counter( *this, ib, a, num_actors );
			boost::mpl::for_each<Actors>( counter );
			return num_actors;
		}

		//////////////////////////// actor iteration //////////////////////////////////

		template<class Actor>
		std::pair<
			SceneIter1B<This,Actor,Symmetric,ActorCopy>,
			SceneIter1B<This,Actor,Symmetric,ActorCopy>
		>
		get_actors() const {
			typedef SceneIter1B<This,Actor,Symmetric,ActorCopy> Iter;
			Iter beg = Iter::make_begin(*this);
			Iter end = Iter::make_end  (*this);
			return std::make_pair(beg,end);
		}
		template<class Actor>
		std::pair<
			SceneIter1B<This,Actor,NotSymmetric,ActorCopy>,
			SceneIter1B<This,Actor,NotSymmetric,ActorCopy>
		>
		get_actors_asym() const {
			typedef SceneIter1B<This,Actor,NotSymmetric,ActorCopy> Iter;
			Iter beg = Iter::make_begin(*this);
			Iter end = Iter::make_end  (*this);
			return std::make_pair(beg,end);
		}
		template<class Actor>
		std::pair<
			SceneIter1B<This,Actor,Symmetric,PlaceHolder>,
			SceneIter1B<This,Actor,Symmetric,PlaceHolder>
		>
		get_actors_placeholder() const {
			typedef SceneIter1B<This,Actor,Symmetric,PlaceHolder> Iter;
			Iter beg = Iter::make_begin(*this);
			Iter end = Iter::make_end  (*this);
			return std::make_pair(beg,end);
		}
		template<class Actor>
		std::pair<
			SceneIter1B<This,Actor,NotSymmetric,PlaceHolder>,
			SceneIter1B<This,Actor,NotSymmetric,PlaceHolder>
		>
		get_actors_placeholder_asym() const {
			typedef SceneIter1B<This,Actor,NotSymmetric,PlaceHolder> Iter;
			Iter beg = Iter::make_begin(*this);
			Iter end = Iter::make_end  (*this);
			return std::make_pair(beg,end);
		}

	///////////////////////////// Visitation, most efficient ///////////////////////////////////////

		template<class Visitor>
		typename boost::disable_if<util::meta::is_pair<typename Visitor::Interaction> >::type
		visit(Visitor & visitor) const {
			typedef typename Visitor::Interaction Actor;
			typedef typename f::result_of::value_at_key<Conformation,Actor>::type Container;
			Index const NBOD = (Index)bodies_.size();
			bool const symmetric = impl::get_Symmetric_true_<Visitor>::type::value;
			Index const NSYM = symmetric ? (Index)this->symframes_.size() : 1;
			for(Index i1 = 0; i1 < NBOD*NSYM; ++i1){
				Conformation const & c =   conformation(i1);
				Position     const & p = this->position(i1);
				Container const & container = c.template get<Actor>();
				BOOST_FOREACH( Actor const & a_0, container ){
					visit_1b_inner(visitor,a_0,p);
				}
			}

		}
		template<class Visitor, class Actor>
		typename boost::enable_if< impl::has_type_Position<Actor> >::type
		visit_1b_inner( Visitor & visitor, Actor const & a_0, Position const & p ) const {
			visitor( Actor(a_0,p) );
		}
		template<class Visitor, class Actor>
		typename boost::disable_if< impl::has_type_Position<Actor> >::type
		visit_1b_inner( Visitor & visitor, Actor const & a_0, Position const &   ) const {
			visitor( a_0 );
		}


		template<class Visitor>
		typename boost::enable_if<util::meta::is_pair<typename Visitor::Interaction> >::type
		visit(Visitor & visitor) const {
			typedef typename Visitor::Interaction::first_type Actor1;
			typedef typename Visitor::Interaction::second_type Actor2;
			typedef typename f::result_of::value_at_key<Conformation,Actor1>::type Container1;
			typedef typename f::result_of::value_at_key<Conformation,Actor2>::type Container2;
			typedef util::container::ContainerInteractions<typename Scene::Position,Container1,Container2,Index> ContInter;
			typedef typename ContInter::Range ContRange;
			typename util::container::get_citer<ContRange>::type iter,end;
			ContRange range;

			Index const NBOD = (Index)bodies_.size();
			Index const NSYM = (Index)this->symframes_.size();
			for(Index i1 = 0; i1 < NBOD*NSYM; ++i1){
				Conformation const & c1 =   conformation(i1);
				Position     const & p1 = this->position(i1);
				for(Index i2 = 0; i2 < NBOD*NSYM; ++i2){
					if( i1 >= NBOD && i2 >= NBOD ) continue;
					if( i1==i2 || (boost::is_same<Actor1,Actor2>::value && i2 <= i1) ) continue;
					Conformation const & c2 =   conformation(i2);
					Position     const & p2 = this->position(i2);

					// // simple version, ~10-15% faster in simple case
					// Index const NACT1 = c1.template get<Actor1>().size();
					// Index const NACT2 = c2.template get<Actor2>().size();
					// for(Index j1 = 0; j1 < NACT1; ++j1){
					// 	Actor1 a1( c1.template get<Actor1>()[j1], p1 );
					// 	for(Index j2 = 0; j2 < NACT2; ++j2){
					// 		Actor2 a2( c2.template get<Actor2>()[j2], p2 );
					// 		visitor( std::make_pair(a1,a2), i1<NBOD&&i2<NBOD?1.0:0.5 );
					// 	}
					// }

					Container1 const & container1 = c1.template get<Actor1>();
					Container2 const & container2 = c2.template get<Actor2>();
					Position const rel_pos = inverse(__position_unsafe__(i1))*( __position_unsafe__(i2) );
					ContInter::get_interaction_range( rel_pos, container1, container2, range );
					for( iter = get_cbegin(range),end  = get_cend(range); iter != end; ++iter){
						double const w = i1<NBOD&&i2<NBOD?1.0:0.5;
						Index j1,j2;
						boost::tie(j1,j2) = *iter;
						Actor1 const & a1_0( c1.template get<Actor1>()[j1] );
						Actor2 const & a2_0( c2.template get<Actor2>()[j2] );
						visit_2b_inner(visitor,a1_0,a2_0,p1,p2,rel_pos,w);
					}

				}
			}

		}

		///@brief visit_2b_inner specialization handles case where both actors are positionable (not fixed)
		template<class Visitor, class Actor1, class Actor2>
		typename boost::enable_if<
			m::and_<
				impl::get_RequireAbsolutePositioning_false_<Visitor>,
				impl::has_type_Position<Actor1> ,
				impl::has_type_Position<Actor2>
			> >::type
		visit_2b_inner(
			Visitor & visitor,
			Actor1 const & a1_0,
			Actor2 const & a2_0,
			Position const & p1,
			Position const & p2,
			Position const & ,
			double w
		) const {
			// cout << "ABSOLUTE POSITION" << endl;
			Actor1 a1( a1_0, p1 );
			Actor2 a2( a2_0, p2 );
			//TODO: wrap this call in a check for has_const_call_oper_3<Visitor,void,Actor1 const &,Actor2 const &,double>
			visitor.template operator()< std::pair<Actor1,Actor2> >( a1, a2, w );
		}

		///@brief visit_2b_inner specialization handles case where Actor1 is fixed
		template<class Visitor, class Actor1, class Actor2>
		typename boost::enable_if<
			m::and_<
				m::not_<impl::get_RequireAbsolutePositioning_false_<Visitor> >,
				impl::has_type_Position<Actor2>
			> >::type
		visit_2b_inner(
			Visitor & visitor,
			Actor1 const & a1_0,
			Actor2 const & a2_0,
			Position const & ,
			Position const & ,
			Position const & rel_pos,
			double w
		) const {
			// cout << "FIXED Actor1" << endl;
			Actor2 a2( a2_0, rel_pos );
			// TODO: remove requirement to form pair; will be more efficient if fixed a1 is passet through w/o copy
			// visitor( std::make_pair(a1_0,a2), w );
			visitor.template operator()< std::pair<Actor1,Actor2> >( a1_0, a2, w );
		}

		template<class Visitor, class Actor1, class Actor2>
		typename boost::enable_if<
			m::and_<
				        impl::has_type_Position<Actor1>   ,
				m::not_<impl::has_type_Position<Actor2> >
			> >::type
		visit_2b_inner(
			Visitor & visitor,
			Actor1 const & a1_0,
			Actor2 const & a2_0,
			Position const & ,
			Position const & ,
			Position const & rel_pos,
			double w
		) const {
			// cout << "FIXED Actor2" << endl;
			Actor1 a1( a1_0, inverse(rel_pos) );
			visitor.template operator()< std::pair<Actor1,Actor2> >( a1, a2_0, w );
		}

		///@brief visit_2b_inner specialization is an error and will not compile
		template<class Visitor, class Actor1, class Actor2>
		typename boost::enable_if<
			m::and_<
				m::not_<impl::has_type_Position<Actor1> > ,
				m::not_<impl::has_type_Position<Actor2> >
			> >::type
		visit_2b_inner(
			Visitor & ,
			Actor1 const & ,
			Actor2 const & ,
			Position const & ,
			Position const & ,
			double
		) const {
			cout << "This SHOULD NEVER HAPPEN" << endl;
			BOOST_STATIC_ASSERT((
				m::or_<
					impl::has_type_Position<Actor1> ,
					impl::has_type_Position<Actor2>
				>::value
			));
		}

	/////////////////////////// interaction iteration //////////////////////////////////

		// TODO: fix this!
		template<class Interaction> struct
		has_interaction : m::bool_<true> {}; // f::result_of::has_key<MAP,Interaction>::value> {};

		template<class Interaction>	struct
		interaction_placeholder_type {typedef typename impl::get_placeholder_type<Interaction,Index>::type type; };

		template<class Interaction> struct
		iter_type { typedef typename
			m::if_< util::meta::is_pair<Interaction>,
				typename m::if_< util::meta::is_homo_pair<Interaction>,
					SceneIter2B<This,Interaction,CountPairUpperTriangle>,
					SceneIter2B<This,Interaction,CountPairNoDuplicates>
				>::type,
				SceneIter1B<This,Interaction,NotSymmetric,PlaceHolder>
			>::type
			type; };

		template<class Interaction> struct
		interactions_type {
			typedef typename iter_type<Interaction>::type Iter;
			typedef std::pair<Iter,Iter> type;
		};
		template<class Interaction>
		typename interactions_type<Interaction>::type
		get_interactions() const {
			typedef typename iter_type<Interaction>::type Iter;
			Iter beg = Iter::make_begin(*this);
			Iter end = Iter::make_end  (*this);
			return std::make_pair(beg,end);
		}

		///@brief gets a concrete interaction (Actor or pair of Actors) from a SceneIter
		template<class Iter>
		typename Iter::Interaction
		get_interaction(
			Iter const & iter
		) const {
			return get_interaction<typename Iter::Interaction>(*iter);
		}

		///@brief gets a concrete interaction (Actor or pair of Actors) from a SceneIter
		template<class Iter>
		typename Iter::Interaction
		get_interaction_absolute(
			Iter const & iter
		) const {
			return get_interaction_absolute<typename Iter::Interaction>(*iter);
		}

		///@brief gets a concrete interaction (Actor or pair of Actors) from a placeholder
		template<class Interaction>
		typename boost::disable_if< util::meta::is_pair<Interaction>, Interaction >::type
		get_interaction(
			typename impl::get_placeholder_type<Interaction,Index>::type const & ph
		) const {
			Conformation const & c =   conformation(ph.first);
			Position     const & p = this->position(ph.first);
			impl::AccessVisitor<Interaction> visitor;
			visit_1b_inner( visitor, c.template get<Interaction>()[ph.second], p );
			return visitor.get();
		}

		///@brief gets a concrete interaction (Actor or pair of Actors) from a placeholder
		template<class Interaction>
		typename boost::disable_if< util::meta::is_pair<Interaction>, Interaction >::type
		get_interaction_absolute(
			typename impl::get_placeholder_type<Interaction,Index>::type const & ph
		) const {
			return get_interaction(ph);
		}

		template<class Interaction,class Visitor>
		typename boost::enable_if< util::meta::is_pair<Interaction>, Interaction >::type
		get_interaction_from_placeholder_2b(
			typename impl::get_placeholder_type<Interaction,Index>::type const & ph
		) const {
			typedef typename Interaction::first_type Actor1;
			typedef typename Interaction::second_type Actor2;
			Index i1 = ph.first .first, i2 = ph.first .second;
			Index j1 = ph.second.first, j2 = ph.second.second;
			Actor1 const & a1_0 = conformation(i1).template get<Actor1>()[j1];
			Actor2 const & a2_0 = conformation(i2).template get<Actor2>()[j2];
			Position const & p1 = this->position(i1);
			Position const & p2 = this->position(i2);
			Position const rel_pos = inverse(__position_unsafe__(i1))*( __position_unsafe__(i2) );
			Visitor visitor;
			double dummy=0;
			visit_2b_inner(visitor,a1_0,a2_0,p1,p2,rel_pos,dummy);
			return visitor.get();
		}

		///@brief gets a concrete interaction (Actor or pair of Actors) from a placeholder
		template<class Interaction>
		typename boost::enable_if< util::meta::is_pair<Interaction>, Interaction >::type
		get_interaction(
			typename impl::get_placeholder_type<Interaction,Index>::type const & ph
		) const {
			return get_interaction_from_placeholder_2b<Interaction,impl::AccessVisitor<Interaction> >(ph);
		}

		///@brief gets a concrete interaction (Actor or pair of Actors) from a placeholder
		template<class Interaction>
		typename boost::enable_if< util::meta::is_pair<Interaction>, Interaction >::type
		get_interaction_absolute(
			typename impl::get_placeholder_type<Interaction,Index>::type const & ph
		) const {
			return get_interaction_from_placeholder_2b<Interaction,impl::AccessVisitorAbsolute<Interaction> >(ph);
		}

		///@brief default weight is 1.0
		template<class PlaceHolder>
		double get_weight_from_placeholder(PlaceHolder const &) const {
			return 1.0;
		}
		///@brief for symmetry, if either body not in asym unit, weight should be 0.5
		double get_weight_from_placeholder(std::pair<std::pair<Index,Index>,std::pair<Index,Index> > const & ph ) const {
			return ( ph.first.first < this->nbodies_asym() && ph.first.second < this->nbodies_asym() ) ? 1.0 : 0.5;
		}


	};

	template< class Scene, class MetaData >
	struct ActorDumpPDB {
		Scene const & scene;
		std::ostream & out;
		MetaData const & meta;
		ActorDumpPDB(
			Scene const & s,
			std::ostream & o,
			MetaData const & m
		) : scene(s), out(o), meta(m) {}
		template<class Actor>
		void operator()( Actor const & ){
			typedef typename Scene::Index Index;
			for( Index ib = 0; ib < scene.n_bodies_; ++ib ){
				for( Index ia = 0; ia < scene.template num_actors<Actor>(ib); ++ia ){
					Actor const & a( scene.template get_actor<Actor>(ib,ia) );
					write_pdb( out, a, meta );
				}
			}
		}
	};


	template< class Scene, class MetaData >
	void write_pdb( std::ostream & out, Scene const & scene, MetaData const & meta ){
		ActorDumpPDB<Scene,MetaData> dumper( scene, out, meta );
		boost::mpl::for_each< typename Scene::Actors >( dumper );
	}


}
}

#endif
