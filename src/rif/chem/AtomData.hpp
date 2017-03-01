#ifndef INCLUDED_chemical_AtomData_HH
#define INCLUDED_chemical_AtomData_HH

#include <string>

#include <iostream>

namespace rif {
namespace chemical {

struct AtomData {
  // ATOM   1608  CD  GLU A 220       3.348   8.506  -0.163  1.00 47.39 C
  std::string atomname, resname, elem;
  int atomnum, resnum;
  char chain;
  bool ishet;
  float occ, bfac;
  static std::string default_atomname() { return "ATOM"; }
  static std::string default_resname() { return "RES"; }
  static char default_chain() { return 'A'; }
  static int default_resnum() { return 0; }
  static int default_atomnum() { return 0; }
  static std::string default_elem() { return "ELEM"; }
  static bool default_ishet() { return false; }
  static float default_occ() { return 1.0; }
  static float default_bfac() { return 0.0; }
  AtomData(std::string const &_atomname = AtomData::default_atomname(),
           std::string const &_resname = AtomData::default_resname(),
           char _chain = AtomData::default_chain(),
           int _resnum = AtomData::default_resnum(),
           int _atomnum = AtomData::default_atomnum(),
           std::string const &_elem = AtomData::default_elem(),
           bool _ishet = AtomData::default_ishet(),
           float _occ = AtomData::default_occ(),
           float _bfac = AtomData::default_bfac())
      : atomname(_atomname),
        resname(_resname),
        elem(_elem),
        atomnum(_atomnum),
        resnum(_resnum),
        chain(_chain),
        ishet(_ishet),
        occ(_occ),
        bfac(_bfac) {}

  bool operator==(AtomData const &o) const {
    return atomname == o.atomname && resname == o.resname && elem == o.elem &&
           atomnum == o.atomnum && resnum == o.resnum && chain == o.chain &&
           ishet == o.ishet && occ == o.occ && bfac == o.bfac;
  }
};

inline std::ostream &operator<<(std::ostream &out, AtomData const &a) {
  if (a.atomname != AtomData::default_atomname()) out << a.atomname << ", ";
  if (a.resname != AtomData::default_resname()) out << a.resname << ", ";
  if (a.chain != AtomData::default_chain()) out << a.chain << ", ";
  if (a.resnum != AtomData::default_resnum()) out << a.resnum << ", ";
  if (a.atomnum != AtomData::default_atomnum()) out << a.atomnum << ", ";
  if (a.elem != AtomData::default_elem()) out << a.elem << ", ";
  if (a.ishet != AtomData::default_ishet()) out << a.ishet << ", ";
  if (a.occ != AtomData::default_occ()) out << a.occ << ", ";
  if (a.bfac != AtomData::default_bfac()) out << a.bfac << ", ";
  return out;
}
}
}

#endif
