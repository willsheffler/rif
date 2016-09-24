#ifndef INCLUDED_scheme_nest_maps_util_get_quaternion_component_HH
#define INCLUDED_scheme_nest_maps_util_get_quaternion_component_HH

#include "riflib/util/meta/util.hh"
#include <boost/utility/enable_if.hpp>

namespace scheme {
namespace nest {
namespace pmap {


/// all this get_quaternion_component mess is to handle eigen quaternions
template<class Q, class Index>
typename boost::enable_if<
	util::meta::has_const_subscript_oper<Q,typename Q::value_type const &,Index> ,
	typename Q::value_type const &    >::type
get_quaternion_component(Q const & q, Index i){ return q[i]; }

template<class Q, class Index>
typename boost::enable_if<
	util::meta::has_subscript_oper<Q,typename Q::value_type &,Index> ,
	typename Q::value_type &    >::type
get_quaternion_component(Q & q, Index i){ return q[i]; }

template<class Q, class Index>
typename boost::disable_if<
	util::meta::has_const_subscript_oper<Q,typename Q::Scalar const &,Index> ,
	typename Q::Scalar const &    >::type
get_quaternion_component(Q const & q, Index i){ return q.coeffs()[i]; }

template<class Q, class Index>
typename boost::disable_if<
	util::meta::has_subscript_oper<Q,typename Q::Scalar &,Index> ,
	typename Q::Scalar &    >::type
get_quaternion_component(Q & q, Index i){ return q.coeffs()[i]; }

}
}
}

#endif
