#ifndef INCLUDED_actor_ActorConcept_io_HH
#define INCLUDED_actor_ActorConcept_io_HH

#include "riflib/actor/ActorConcept.hh"
#include <iostream>

namespace scheme {
namespace actor {



	template< class P, class D >
	std::ostream & operator<<(std::ostream & out, ActorConcept<P,D> const & a){
		return out << "Actor( " << a.position_ << ", " << a.data_ << " )";
	}

	// template< class P >
	// std::ostream & operator<<(std::ostream & out, ActorConcept<P,int> const & a){
	// 	return out << "ADI( " << a.position_ << ", " << a.data_ << " )";
	// }

	// template< class P >
	// std::ostream & operator<<(std::ostream & out, ActorConcept<P,char> const & a){
	// 	return out << "ADC( " << a.position_ << ", " << a.data_ << " )";
	// }



}
}

#endif
