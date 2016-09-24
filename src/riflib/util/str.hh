#ifndef INCLUDED_scheme_util_str_HH
#define INCLUDED_scheme_util_str_HH


#include <string>
#include <sstream>

namespace scheme {

template<class T>
std::string str(T const & t){
	std::ostringstream oss;
	oss << t;
	return oss.str();
}

}

#endif
