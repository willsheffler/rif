#ifndef INCLUDED_rosetta_atype_map_HH
#define INCLUDED_rosetta_atype_map_HH

#include <string>

namespace scheme {
namespace rosetta {

int rosetta_atom_type(std::string const& resname, std::string const& atomname);
}
}

#endif
