#include <gtest/gtest.h>

#include "dump_pdb_atom.hpp"

namespace scheme { namespace io { namespace test {

using std::cout;
using std::endl;

TEST(IO,dum_pdb_atom){
	AtomData a;
	std::ostringstream oss1,oss2;
	dump_pdb_atom(oss1,0,0,0,a);
	oss2 << "ATOM      0 ATOM RES A   0       0.000   0.000   0.000  1.00  0.00        ELEM" << endl;
	ASSERT_EQ( oss1.str() , oss2.str() );
}

}}}
