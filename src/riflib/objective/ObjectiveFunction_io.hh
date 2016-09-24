#ifndef INCLUDED_objective_ObjectiveFunction_io_HH
#define INCLUDED_objective_ObjectiveFunction_io_HH

#include "riflib/util/meta/print_type.hh"

// #define DEBUG_IO 

namespace scheme {
namespace objective {

namespace m = boost::mpl;
namespace f = boost::fusion;

template<class O,class C> struct ObjectiveFunction;

template<typename O,typename C>
std::ostream & operator<<(std::ostream & out, ObjectiveFunction<O,C> const & obj){
	typedef ObjectiveFunction<O,C> OBJ;
	out << "ObjectiveFunction" << std::endl;
	out << "    Interactions:" << std::endl;
		m::for_each<typename OBJ::InteractionTypes>(util::meta::PrintType(out,"        "));
	out << "    Objectives:" << std::endl;		
	f::for_each( (typename OBJ::ObjectiveMap::FusionType&)obj.objective_map_, util::meta::PrintBFMapofVec(out,"        ") );
	out << "    RAW:" << std::endl;
		m::for_each<typename OBJ::Objectives>(util::meta::PrintType(out,"        "));
	return out;
}

}
}

#endif
