#include <stdint.h>

#ifndef CXX11
#include <boost/static_assert.hpp>
BOOST_STATIC_ASSERT_MSG(false, "c++11 required for rif stuff");
#endif
// #ifndef CXX14
//  #include <boost/static_assert.hpp>
//  BOOST_STATIC_ASSERT_MSG( false, "c++14 required for rif stuff");
// #endif

#include <memory>

namespace rif {

using std::shared_ptr;
using std::unique_ptr;
using std::weak_ptr;
using std::make_shared;
using std::enable_shared_from_this;
}
