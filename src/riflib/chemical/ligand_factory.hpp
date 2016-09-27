#ifndef INCLUDED_chemical_ligand_factory_HH
#define INCLUDED_chemical_ligand_factory_HH

#include "riflib/rosetta/atype_map.hpp"

#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

namespace scheme { namespace chemical {

std::vector<std::string>
get_pdb_lines(
	std::string resn,
	bool hydrogen
);

struct RosettaAtypeMap {
	static int get_atom_type(
		std::string const & aname,
		std::string const & rname
	){
		return rosetta::rosetta_atom_type(aname,rname);
	}
};


template<class Atom, class AtypeMap = RosettaAtypeMap>
struct LigandFactory {

	Atom make_atom_pdbline(std::string const & _line, int atype = -1){
		/*
ATOM      1  N1  BTN X   1       0.696 -12.422   3.375  1.00 20.00           N 	*/
		std::string line = boost::trim_copy(_line);
		int const N = (int)line.size();
		BOOST_VERIFY( N > 53 );
		double x = boost::lexical_cast<double>( boost::trim_copy( line.substr(28,10) ) );
		double y = boost::lexical_cast<double>( boost::trim_copy( line.substr(38,8) ) );
		double z = boost::lexical_cast<double>( boost::trim_copy( line.substr(46,8) ) );
		std::string atomname = boost::trim_copy( line.substr(12,4) );
		std::string resname  = boost::trim_copy( line.substr(16,4) );
		std::string elem = ""; if(N > 68) elem = boost::trim_copy( line.substr(68) );
		float occ  = 1.0; if(N > 58) occ  = boost::lexical_cast<float>( boost::trim_copy( line.substr(54,6) ) );
		float bfac = 0.0; if(N > 64) bfac = boost::lexical_cast<float>( boost::trim_copy( line.substr(60,6) ) );
		if(atype==-1){
			atype = AtypeMap::get_atom_type(resname,atomname);
			// BOOST_VERIFY(atype >= 0);
			// std::cout << resname << " " << atomname << " " << atype << std::endl;
		}
		return Atom(
			typename Atom::Position(x,y,z),
			atype,
			atomname,
			resname,
			line[21],
			boost::lexical_cast<int>( boost::trim_copy( line.substr(22,6) ) ),
			boost::lexical_cast<int>( boost::trim_copy( line.substr( 6,6) ) ),
			elem,
			line.substr(0,6)=="HETATM",
			occ,
			bfac
		);

	}

	template<class Iter>
	void make_biotin_minimal(Iter i){
		*i++ = make_atom_pdbline("ATOM      1  N1  BTN X   1       0.696 -12.422   3.375  1.00 20.00           N", 7  );
		*i++ = make_atom_pdbline("ATOM      2  S1  BTN X   1       0.576  -9.666   5.336  1.00 20.00           S", 17 );
		*i++ = make_atom_pdbline("ATOM      3  C1  BTN X   1      -0.523 -10.824   6.189  1.00 20.00           C", 3  );
		*i++ = make_atom_pdbline("ATOM      4  N2  BTN X   1      -1.324 -12.123   4.201  1.00 20.00           N", 7  );
		*i++ = make_atom_pdbline("ATOM      5  C2  BTN X   1      -0.608 -12.327   3.072  1.00 20.00           C", 3  );
		*i++ = make_atom_pdbline("ATOM      6  O1  BTN X   1      -1.125 -12.422   1.933  1.00 20.00           O", 13 );
		*i++ = make_atom_pdbline("ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C", 3  );
		*i++ = make_atom_pdbline("ATOM      8  C4  BTN X   1       0.953 -12.267   4.780  1.00 20.00           C", 6  );
		*i++ = make_atom_pdbline("ATOM      9  C5  BTN X   1       1.765 -11.040   5.134  1.00 20.00           C", 4  );
		*i++ = make_atom_pdbline("ATOM     10  C6  BTN X   1      -1.836 -10.395   6.850  1.00 20.00           C", 4  );
	}

	template<class Iter>
	void make_atoms( Iter outiter, std::string resn, bool hydrogen = true ){
		std::vector<std::string> lines = get_pdb_lines(resn,hydrogen);
		BOOST_FOREACH(std::string const & line,lines){
			*outiter++ = make_atom_pdbline(line);
		}
	}
};

}}

#endif
