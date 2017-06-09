#ifndef INCLUDED_actor_Atom_HH
#define INCLUDED_actor_Atom_HH

#include <iostream>
#include "chem/AtomData.hpp"
#include "io/dump_pdb_atom.hpp"
#include "numeric/util.hpp"
#include "types.hpp"

namespace rif {
namespace actor {

template <class _Position>
struct Atom {
  using Position = _Position;
  using THIS = Atom<Position>;
  using Scalar = typename Position::Scalar;

  Atom() : pos(0, 0, 0), atype(-1), rtype(-1), anum(-1) {}

  template <class P>
  Atom(P p, int16_t type = -1, int8_t restype = -1, int8_t atomnum = -1)
      : pos(p[0], p[1], p[2]), atype(type), rtype(restype), anum(atomnum) {}

  template <class Xform>
  Atom(Atom a, Xform moveby) {
    pos = moveby * a.position();
    atype = a.atype;
    rtype = a.rtype;
    anum = a.anum;
  }

  Atom(Position p, Atom a) {
    pos = p;
    atype = a.atype;
    rtype = a.rtype;
    anum = a.anum;
  }

  int type() const { return atype; }
  int restype() const { return rtype; }
  int atomnum() const { return anum; }
  void set_type(int16_t i) { atype = i; }
  void set_restype(int8_t i) { rtype = i; }
  void set_atomnum(int8_t i) { anum = i; }

  Position position() const { return pos; }

  template <class P>
  void set_position(P const &pos) {
    pos[0] = pos[0];
    pos[1] = pos[1];
    pos[2] = pos[2];
  }

  bool operator==(Atom<Position> o) const {
    return numeric::approx_eq(o.pos, pos) && o.atype == atype &&
           o.rtype == rtype && o.anum == anum;
  }

  auto &operator[](size_t i) { return pos[i]; }
  auto const &operator[](size_t i) const { return pos[i]; }

  Position pos;
  int8_t atype;
  int8_t anum;
  int16_t rtype;
};

template <class P>
std::ostream &operator<<(std::ostream &out, Atom<P> const &x) {
  return out << "Atom( " << x.position() << ", " << x.type() << " )";
}

template <class P, class MetaData>
void write_pdb(std::ostream &out, Atom<P> const &a, MetaData const &meta) {
  io::dump_pdb_atom(out, a.position(),
                    meta.atom_data(a.restype(), a.atomnum()));
}

template <class V>
Atom<V> operator+(Atom<V> a, V v) {
  Atom<V> b(a);
  b.pos += v;
  return b;
}
template <class V>
Atom<V> operator+(V v, Atom<V> a) {
  Atom<V> b(a);
  b.pos += v;
  return b;
}
template <class V>
Atom<V> operator-(Atom<V> a, V v) {
  Atom<V> b(a);
  b.pos -= v;
  return b;
}
template <class V>
Atom<V> operator-(V v, Atom<V> a) {
  Atom<V> b(a);
  b.pos = v - a.pos;
  return b;
}
template <class V>
Atom<V> operator*(Atom<V> a, float f) {
  Atom<V> b(a);
  b.pos = f * a.pos;
  return b;
}
template <class V>
Atom<V> operator/(Atom<V> a, float f) {
  Atom<V> b(a);
  b.pos = a.pos / f;
  return b;
}
template <class M, class V>
Atom<V> operator*(M m, Atom<V> a) {
  Atom<V> b(a);
  b.pos = m * a.pos;
  return b;
}

// using chemical::AtomData;

// template <class _Position>
// struct AtomWithData {
//   typedef _Position Position;

//   AtomWithData()
//       : pos(0, 0, 0), atype(0), data_(rif::make_shared<AtomData>()) {}
//   // data_(new AtomData) {}

//   template <class P>
//   AtomWithData(P const &p, int type = 0,
//                std::string const &atomname = AtomData::default_atomname(),
//                std::string const &resname = AtomData::default_resname(),
//                char chain = AtomData::default_chain(),
//                int resnum = AtomData::default_resnum(),
//                int atomnum = AtomData::default_atomnum(),
//                std::string const &elem = AtomData::default_elem(),
//                bool ishet = AtomData::default_ishet(),
//                float occ = AtomData::default_occ(),
//                float bfac = AtomData::default_bfac())
//       : pos(p[0], p[1], p[2]),
//         atype(type),
//         data_(rif::make_shared<AtomData>(  // TODO: why do I need rif:: here?
//             // data_(new AtomData(
//             atomname, resname, chain, resnum, atomnum, elem, ishet, occ,
//             bfac)) {}

//   // // TODO: compare raw ptr for speed
//   // ~AtomWithData(){ delete data_; }

//   template <class Xform>
//   AtomWithData(AtomWithData const &a, Xform const &moveby) {
//     pos = moveby * a.position();
//     atype = a.atype;
//     data_ = a.data_;
//   }

//   int type() const { return atype; }
//   void set_type(int i) { atype = i; }
//   AtomData const &data() const { return *data_; }
//   AtomData &nonconst_data() { return *data_; }

//   Position const &position() const { return pos; }

//   template <class P>
//   void set_position(P const &pos) {
//     pos[0] = pos[0];
//     pos[1] = pos[1];
//     pos[2] = pos[2];
//   }

//   bool operator==(AtomWithData<Position> const &o) const {
//     return numeric::approx_eq(o.pos, pos) && o.atype == atype &&
//            o.data_ == data_;
//   }

//  private:
//   int atype;
//   shared_ptr<AtomData> data_;
//   // AtomData *data_;
//   Position pos;
// };

// template <class P>
// std::ostream &operator<<(std::ostream &out, AtomWithData<P> const &x) {
//   return out << "AtomWithData( " << x.position() << ", " << x.type() << ", "
//              << x.data() << " )";
// }

// template <class P, class MetaData>
// void write_pdb(std::ostream &out, AtomWithData<P> const &a, MetaData const &)
// {
//   io::dump_pdb_atom(out, a.position(), a.data());
// }
}
}

namespace std {
template <class P>
struct is_pod<rif::actor::Atom<P>> : public std::integral_constant<bool, true> {
};
}

#endif
