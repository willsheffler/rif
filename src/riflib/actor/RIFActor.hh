#ifndef INCLUDED_actor_RIFActor_HH
#define INCLUDED_actor_RIFActor_HH

#include "riflib/objective/hash/XformMap.hh"

namespace scheme {
namespace actor {




	template<class _Position, class _XformMap>
	struct RIFActor {

		///@brief Position type, leave out to make actor "Fixed"
		typedef _Position Position;
		typedef _XformMap XformMap;
		typedef RIFActor<Position,XformMap> THIS;

		Position position_;
		scheme::shared_ptr<XformMap> xform_map_;


		RIFActor() : position_(), xform_map_() {}

		RIFActor(Position const & p, shared_ptr<XformMap> xmap) : position_(p), xform_map_(xmap) {}

		RIFActor(
			RIFActor const & actor0,
			Position const & moveby
		){
			xform_map_ = actor0.xform_map_;
			set_position(moveby*actor0.position());
		}

		void 
		set_position(
			Position const & pos
		){ position_ = pos; }

		Position const &
		position() const { return position_; }

		bool operator==(THIS const & o) const { return o.position_==position_ && o.xform_map_==xform_map_; }

		///@brief necessary for testing only
		bool operator<(THIS const & o) const { return std::make_pair(position_,xform_map_) < std::make_pair(o.position_,o.xform_map_); }

		///@brief necessary for testing only
	  	template<class Archive> void serialize(Archive & ar, const unsigned int ){
	  		ar & position_;
	  		ar & xform_map_;
	  	}

	};


template<class X, class M>
std::ostream & operator<<(std::ostream & out,RIFActor<X,M> const& si){ return out << "RIFActor"; }

}
}

#endif
