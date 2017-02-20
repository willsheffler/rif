#ifndef INCLUDED_actor_VoxelActor_HH
#define INCLUDED_actor_VoxelActor_HH

#include "objective/voxel/VoxelArray.hpp"

namespace scheme {
namespace actor {

template <class _Position, class _Float>
struct VoxelActor {
  typedef VoxelActor<_Position, _Float> THIS;

  ///@brief Position type, leave out to make actor "Fixed"
  // typedef _Position Position;
  typedef _Float Float;
  typedef objective::voxel::VoxelArray<3, Float> VoxelArray;
  // typedef std::vector<std::vector<shared_ptr< VoxelArray > > > Voxels;
  typedef std::vector<std::vector<VoxelArray*> > Voxels;

  // Position position_;
  Voxels voxels_;

  VoxelActor() {}

  VoxelActor(Voxels const& v) : voxels_(v) {}

  // VoxelActor( Position const & p, Voxels const * v ) :  voxels_(v) {}

  // VoxelActor(
  // 	THIS const & actor0,
  // 	Position const & moveby
  // ){
  // 	voxels_ = actor0.voxels_; // better to use ptr? 3 words vs 1?
  // 	set_position(moveby*actor0.position());
  // }

  // void
  // set_position(
  // 	Position const & pos
  // ){ position_ = pos; }

  // Position const &
  // position() const { return position_; }

  Voxels const& voxels() const { return voxels_; }

  // bool operator==(THIS const & o) const { return o.position_==position_ &&
  // o.voxels_==voxels_; }

  // ///@brief necessary for testing only
  // bool operator<(THIS const & o) const { return
  // std::make_pair(position_,voxels_) < std::make_pair(o.position_,o.voxels_);
  // }

  ///@brief necessary for testing only
  template <class Archive>
  void serialize(Archive& ar, const unsigned int) {
    // ar & position_;
    // ar & voxels_;
    std::exit(-1);
  }
};

template <class X, class F>
std::ostream& operator<<(std::ostream& out, VoxelActor<X, F> const& va) {
  return out << "VoxelActor";
}

template <class X, class F, class MetaData>
void write_pdb(std::ostream&, VoxelActor<X, F> const&, MetaData const&) {
  ;
}

template <class VoxelActor, class Atom, bool REPL_ONLY = false>
struct Score_Voxel_vs_Atom {
  typedef float Result;
  typedef std::pair<VoxelActor, Atom> Interaction;
  static std::string name() { return "Score_Voxel_vs_Atom"; }
  template <class Config>
  Result operator()(VoxelActor const& v, Atom const& a, Config const& c) const {
    // std::cout << "score voxel vs atom " << a.data().atomname << std::endl;
    // std::cout << "   resl " << c << std::endl;
    // std::cout << "   type " << a.type() << std::endl;
    // std::cout << "   pos  " << a.position().transpose() << std::endl;
    // std::cout << "     LB " << v.voxels()[c][a.type()]->lb_ << std::endl;
    // std::cout << "     UB " << v.voxels()[c][a.type()]->ub_ << std::endl;
    float score = v.voxels()[c][a.type()]->at(a.position()[0], a.position()[1],
                                              a.position()[2]);
    // std::cout << "  score " << score << std::endl;
    if (REPL_ONLY)
      return std::max(0.0f, score);
    else
      return a.type() > 17 ? std::max(0.0f, score) : score;
  }
  template <class Pair, class Config>
  Result operator()(Pair const& p, Config const& c) const {
    return this->operator()(p.first, p.second, c);
  }
};
template <class A, class B>
std::ostream& operator<<(std::ostream& out,
                         Score_Voxel_vs_Atom<A, B> const& si) {
  return out << si.name();
}
}
}

#endif
