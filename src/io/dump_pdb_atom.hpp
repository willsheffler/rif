#ifndef INCLUDED_io_dump_pdb_atom_HH
#define INCLUDED_io_dump_pdb_atom_HH

#include "chemical/AtomData.hpp"

#include <boost/assert.hpp>
#include <cassert>
#include <iostream>
#include <sstream>

namespace scheme {
namespace io {

using chemical::AtomData;

inline void dump_pdb_atom(std::ostream &out, double x, double y, double z,
                          AtomData const &a = AtomData()) {
  // std::string atomname,resname,elem;
  // int atomnum,resnum;
  // char chain;
  // bool ishet;
  // float occ,bfac;
  BOOST_VERIFY(a.atomname.size() < 5);
  BOOST_VERIFY(a.resname.size() < 4);
  BOOST_VERIFY(a.elem.size() < 11);
  BOOST_VERIFY(x < 10000 && x > -1000);
  BOOST_VERIFY(y < 10000 && y > -1000);
  BOOST_VERIFY(z < 10000 && z > -1000);
  // cout << "ATOM   1604  C   GLU A 220       5.010  12.933   1.553  1.00 41.10
  // C" << endl;
  char buf[128];
  std::string aname = a.atomname;
  if (aname.size() == 1) aname = aname + "  ";
  if (aname.size() == 2) aname = aname + " ";
  snprintf(buf, 128, "%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
           a.ishet ? "HETATM" : "ATOM  ", a.atomnum, aname.c_str(),
           a.resname.c_str(), a.chain, a.resnum, x, y, z, a.occ, a.bfac,
           a.elem.c_str());
  out << buf;
}

template <class XYZ>
inline void dump_pdb_atom_resname_atomname(std::ostream &out,
                                           std::string const &resname,
                                           std::string const &atomname,
                                           XYZ const &xyz) {
  AtomData a;
  a.resname = resname;
  a.atomname = atomname;
  dump_pdb_atom(out, xyz[0], xyz[1], xyz[2], a);
}

template <class XYZ>
inline void dump_pdb_atom(std::ostream &out, std::string const &elem,
                          XYZ const &xyz) {
  AtomData a;
  a.elem = elem;
  dump_pdb_atom(out, xyz[0], xyz[1], xyz[2], a);
}

template <class XYZ>
inline void dump_pdb_atom(std::ostream &out, std::string elem, int resi,
                          XYZ const &xyz) {
  AtomData a;
  a.elem = elem;
  a.resnum = resi;
  dump_pdb_atom(out, xyz[0], xyz[1], xyz[2], a);
}

template <class XYZ>
inline void dump_pdb_atom(std::ostream &out, int index, XYZ const &xyz) {
  AtomData a;
  a.atomnum = index;
  dump_pdb_atom(out, xyz[0], xyz[1], xyz[2], a);
}

template <class XYZ>
inline void dump_pdb_atom(std::ostream &out, XYZ const &xyz,
                          AtomData const &a) {
  dump_pdb_atom(out, xyz[0], xyz[1], xyz[2], a);
}

template <class XYZ>
inline void dump_pdb_atom(std::ostream &out, XYZ const &xyz, AtomData const &a,
                          int resi) {
  AtomData d(a);
  d.resnum = resi;
  dump_pdb_atom(out, xyz[0], xyz[1], xyz[2], d);
}

template <class Atom>
inline void dump_pdb_atom(std::ostream &out, Atom const &atom, int iatom = -1,
                          int ires = -1, char chain = 0) {
  AtomData d = atom.data();
  if (iatom > 0) d.atomnum = iatom;
  if (ires > 0) d.resnum = ires;
  if (chain != 0) d.chain = chain;
  dump_pdb_atom(out, atom.position()[0], atom.position()[1], atom.position()[2],
                d);
}

template <class Atom>
inline std::string dump_pdb_atom(Atom const &atom) {
  std::ostringstream o;
  dump_pdb_atom(o, atom.position()[0], atom.position()[1], atom.position()[2],
                atom.data());
  return o.str().substr(0, o.str().size() - 1);
}
}
}

#endif
