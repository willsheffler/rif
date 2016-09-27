#ifndef INCLUDED_kinematics_SceneIterator_HH
#define INCLUDED_kinematics_SceneIterator_HH

#include "riflib/types.hpp"
#include "riflib/util/meta/util.hpp"
#include "riflib/util/meta/util.hpp"
#include "riflib/util/meta/print_type.hpp"
#include "riflib/util/meta/InstanceMap.hpp"
#include "riflib/util/container/ContainerInteractions.hpp"

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/mpl.hpp>

#include <boost/iterator/iterator_facade.hpp>

namespace scheme {
namespace kinematics {

namespace m = boost::mpl;
namespace f = boost::fusion;
using std::cout;
using std::endl;

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  1 body /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

	namespace impl {
		SCHEME_HAS_MEMBER_TYPE(Position)
	}

	// SymPolicies
	template<class Scene, class Actor>
	struct Symmetric {
		typedef typename Scene::Index Index;
		typedef typename Scene::Conformation Conf;
		typedef typename Scene::Position Pos;		
		static Index nbodies(Scene const & s) { return s.nbodies(); }
		static Conf const & conformation(Scene const & s, Index i){ return s.__conformation_unsafe__(i); }
		static Pos  const & position(Scene const & s, Index i){ return s.__position_unsafe__(i); }		
	};
	template<class Scene, class Actor>
	struct NotSymmetric {
		typedef typename Scene::Index Index;
		typedef typename Scene::Conformation Conf;
		typedef typename Scene::Position Pos;		
		static Index nbodies(Scene const & s) { return s.nbodies_asym(); }
		static Conf const & conformation(Scene const & s, Index i){ return s.__conformation_unsafe_asym__(i); }
		static Pos  const & position(Scene const & s, Index i){ return s.__position_unsafe_asym__(i); }		
	};

	// Deref Policies
	template<class Scene, class Actor>
	struct PlaceHolder {
		typedef typename Scene::Index Index;
		typedef std::pair<Index,Index> Value;
		static Value deref(Scene const &, Index ib, Index ia){ return Value(ib,ia); }
	};	
	template<class Scene, class Actor>
	struct ActorCopy {
		typedef typename Scene::Index Index;
		typedef Actor Value;
		typedef typename Scene::Position Pos;
		static Value deref(Scene const & s, Index ib, Index ia){ 
			return s.template get_interaction<Actor>(std::make_pair(ib,ia));
		}
	};	

	template<
		class Scene,
		class Interaction,
		template<class,class> class SymPolicy,
		template<class,class> class DerefPolicy
	>
	struct SceneIter1B : boost::iterator_facade<
							SceneIter1B<Scene,Interaction,SymPolicy,DerefPolicy>,
							typename DerefPolicy<Scene,Interaction>::Value const,
							boost::forward_traversal_tag,
							typename DerefPolicy<Scene,Interaction>::Value const
						>					
	{
		BOOST_STATIC_ASSERT(( !util::meta::is_pair<Interaction>::value ));
		typedef SceneIter1B<Scene,Interaction,SymPolicy,DerefPolicy> THIS;
		typedef Interaction Actor;
		typedef SymPolicy<Scene,Interaction> Sym;
		typedef DerefPolicy<Scene,Interaction> Deref;
		typedef typename Scene::Index Index;
		SceneIter1B(){}
		SceneIter1B(Scene const & s, Index ib, Index ia) : scene_(&s),ibody_(ib),iactor_(ia){}
		static THIS make_begin(Scene const & s){ return THIS(s,next_nonempty(s,(Index)0),(Index)0); }
		static THIS make_end  (Scene const & s){ return THIS(s,  Sym::nbodies(s)       ,(Index)0); }
		static Index next_nonempty(Scene const & s, Index i) {
			// cout << "next_nonempty " << i;
			while( i < Sym::nbodies(s) && 0==Sym::conformation(s,i).template get<Actor>().size() ) ++i;
			// cout << " -> " << i << endl;
			return i;
		}
	private:
	    friend class boost::iterator_core_access;
		void
		increment(){
            assert( ibody_ < Sym::nbodies(*scene_) );
			++iactor_;
			if( iactor_ == Sym::conformation(*scene_,ibody_).template get<Actor>().size() ){
				iactor_ = 0;
				ibody_ = next_nonempty(*scene_,++ibody_);
			}
		}
		typename Deref::Value const 
		dereference() const { 
			return Deref::deref(*scene_,ibody_,iactor_);
		}
		bool
		equal(THIS const & o) const {
			return scene_==o.scene_ && ibody_==o.ibody_ && iactor_==o.iactor_;
		}
		Scene const * scene_;
		Index ibody_,iactor_;
	};




///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// 2 body //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
	// TODO: fix end generation for speed and correctness

	template<class I> struct  make_pairpair { typedef std::pair<std::pair<I,I>,std::pair<I,I> > type; };

	struct CountPairNoDuplicates {
		template<class Index>
		static
		void
		next(Index & i, Index & j, size_t s0, size_t s1){
			++j;
			j += i==j ? 1 : 0; // skip iff j- == i
			if( j >= ((i>=s0)?s0:s1) ){
				j = 0;
				++i;
			}
			// cout << "NEXT " << i << " " << j << endl;
		}
		// static size_t get_i_end(size_t n){ return n; }
		// static size_t get_j_end(size_t  ){ return 0; }
		// template<class Index>
		// static void get_end(size_t , size_t s1, Index & i, Index & j){
		// 	i = s1;
		// 	j = 0;
		// }
		static size_t max_body1(size_t , size_t s1){ return s1; }
	};
	struct CountPairUpperTriangle {
		template<class Index>
		static
		void
		next(Index & i, Index & j, size_t , size_t s1){
			++j;
			if( j >= s1 ){
				++i;
				j = i+1;
			}			
		}
		// // static size_t get_i_end(size_t n){ if(n==1) return 1; return n==0?0:n-1; }
		// // static size_t get_j_end(size_t n){ if(n==1) return 2; return n==0?0:n; }
		// template<class Index>
		// static void get_end(size_t s0, size_t s1, Index & i, Index & j){
		// 	if     ( s0==0 ) { i = j = 0; }
		// 	else if( s0==1 ) { i = 1; j = 2; }
		// 	else {
		// 		i = s0-1;
		// 		j = s0;
		// 		i += s0==s1?0:1;
		// 		j += s0==s1?0:1;
		// 	}
		// }
		static size_t max_body1(size_t s0, size_t ){ return s0; }
	};

	template<class _Scene,class _Interaction, class _CountPair>
	struct SceneIter2B : boost::iterator_facade<
							SceneIter2B<_Scene,std::pair<
								typename _Interaction::first_type,
								typename _Interaction::second_type>,
								_CountPair >,
							typename make_pairpair<typename _Scene::Index>::type const,
							boost::forward_traversal_tag,
							typename make_pairpair<typename _Scene::Index>::type const
						>
	{
		typedef _Scene Scene;
		typedef _Interaction Interaction;
		typedef _CountPair CountPair;
		BOOST_STATIC_ASSERT(( util::meta::is_pair<Interaction>::value ));
		typedef typename Interaction::first_type Actor1;
		typedef typename Interaction::second_type Actor2;
		typedef SceneIter2B<Scene,std::pair<Actor1,Actor2>,CountPair> THIS;
		typedef typename Scene::Index Index;
		typedef typename Scene::Conformation Conf;
		typedef typename f::result_of::value_at_key<Conf,Actor1>::type Container1;
		typedef typename f::result_of::value_at_key<Conf,Actor2>::type Container2;
		typedef util::container::ContainerInteractions<typename Scene::Position,Container1,Container2,Index> ContInter;
		typedef typename ContInter::Range ContRange;

		SceneIter2B(){}
		SceneIter2B(Scene const & s, Index ib1, Index ib2) : scene_(&s),ibody1_(ib1),ibody2_(ib2){}
		void
		update_range(){ // X * A = B, X = B * ~A
			assert( ibody1_ < scene_->nbodies_asym() || ibody2_ < scene_->nbodies_asym() );
			ContInter::get_interaction_range(
				scene_->__position_unsafe__(ibody2_) * inverse( scene_->__position_unsafe__(ibody1_) ) ,
				scene_->__conformation_unsafe__(ibody1_).template get<Actor1>()  ,
				scene_->__conformation_unsafe__(ibody2_).template get<Actor2>()  ,
				range_  );
			iter_ = get_cbegin(range_);
			end_  = get_cend(range_);
			// cout << "UPDATE_RANGE " << ibody1_ << " " << ibody2_ << " beg " << *iter_ << " end " << *end_ << endl;
		}
		static
		THIS 
		make_begin(
			Scene const & s
		){
			Index i=0,j=s.nbodies()>1?1:0;
			next_while_empty(s,i,j); 
			THIS beg(s,i,j); 
			// cout << "MAKE BEGIN " << i << " " << j << " " << *beg << endl;
			if( i < (Index)s.nbodies     () && j < (Index)s.nbodies      () && 
			   (i < (Index)s.nbodies_asym() || j < (Index)s.nbodies_asym()) ) beg.update_range();
			return beg;
		}
		static
		THIS
		make_end(
			Scene const & s
		){
			// Index i,j;
			// CountPair::get_end(s.nbodies_asym(),s.nbodies(),i,j);
			return THIS(s,std::numeric_limits<Index>::max(),std::numeric_limits<Index>::max());
		}
		bool is_universal_end() const {
			return ibody1_ == std::numeric_limits<Index>::max() && ibody2_ == std::numeric_limits<Index>::max();
		}
		static
		void
		next_while_empty(
			Scene const & s,
			Index & i,
			Index & j
		) {
			// cout << "next_while_empty " << i<<"-"<<j << endl;
 			while( i < CountPair::max_body1(s.nbodies_asym(),s.nbodies()) && j < s.nbodies() && ( 
 				0==s.__conformation_unsafe__(i).template get<Actor1>().size() ||
			    0==s.__conformation_unsafe__(j).template get<Actor2>().size() )
			){
				// cout << "  NEXT " << i << " " << j;
 				CountPair::next(i,j,s.nbodies_asym(),s.nbodies());
 				// cout << " -> " << i << " " << j << endl;
 			}
			// cout << " -> " << i << "-" << j << endl;
		}
	private:
	    friend class boost::iterator_core_access;
		void
		increment(){
			// cout << "INCREMENT iter_ " << *iter_;
			++iter_;
			// cout << " -> " << *iter_ << " end " << *end_ << endl;
			if( iter_==end_ ){ // end_ of body/body inter, must update range_
				CountPair::next(ibody1_,ibody2_,scene_->nbodies_asym(),scene_->nbodies());
				// cout << "NEXT_NONEMPTY " << ibody1_ << " " << ibody2_;
				next_while_empty(*scene_,ibody1_,ibody2_);
				// cout << " -> " << ibody1_ << " " << ibody2_ << endl;
				if( ibody1_ < scene_->nbodies     () && ibody2_ < scene_->nbodies     () && 
				   (ibody1_ < scene_->nbodies_asym() || ibody2_ < scene_->nbodies_asym() )) update_range();
				// std::cout << "redo range_ " << ibody1_ << " " << ibody2_ << std::endl;
				// std::cout << "SceneIter2B NEW BODIES: " << ibody1_ << " " << ibody2_ << std::endl;
			}
			// not valid when at end
			// assert( ibody1_ < scene_->nbodies_asym() );
			// assert( ibody2_ < scene_->nbodies() );
			// assert( iter_->first  < scene_->body(ibody1_).conformation().template get<Actor1>().size() );
			// assert( iter_->second < scene_->body(ibody2_).conformation().template get<Actor2>().size() );
		}
		typename make_pairpair<Index>::type const
		dereference() const { 
			return std::make_pair(std::make_pair(ibody1_,ibody2_),*iter_);
		}
		bool
		equal(THIS const & o) const {
			assert(scene_==o.scene_);
			// hacky but works for checking if at end
			if( o.is_universal_end() ){
				return ( ( ibody1_ >= scene_->nbodies     () || ibody2_ >= scene_->nbodies     () ) ||
				         ( ibody1_ >= scene_->nbodies_asym() && ibody2_ >= scene_->nbodies_asym() ) );
			} else if( is_universal_end() ) {
				return ( ( o.ibody1_ >= o.scene_->nbodies     () || o.ibody2_ >= o.scene_->nbodies     () ) ||
				         ( o.ibody1_ >= o.scene_->nbodies_asym() && o.ibody2_ >= o.scene_->nbodies_asym() ) );				
			} else {
				std::cerr << "SceneIter2B can only compare to end!" << std::endl;
				std::exit(-1);
			}
			// // std::cout << "EQUAL " << ibody1_ << " " << ibody2_ << " / " << o.ibody1_ << " " << o.ibody2_ << std::endl;
			// return scene_==o.scene_ && ibody1_==o.ibody1_ && ibody2_==o.ibody2_;
		}
		Scene const * scene_;
		Index ibody1_,ibody2_;
		ContRange range_;
		typename util::container::get_citer<ContRange>::type iter_,end_;
	};



}
}

#endif
