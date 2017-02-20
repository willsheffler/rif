#ifndef INCLUDED_scheme_util_str_HH
#define INCLUDED_scheme_util_str_HH

#include <sstream>
#include <string>

namespace scheme {

template <class T>
std::string str(T const& t) {
  std::ostringstream oss;
  oss << t;
  return oss.str();
}
}

#endif
