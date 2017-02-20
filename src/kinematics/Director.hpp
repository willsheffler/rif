#ifndef INCLUDED_kinematics_Director_HH
#define INCLUDED_kinematics_Director_HH

#include "kinematics/Scene.hpp"
#include "nest/MultiNest.hpp"
#include "types.hpp"

#include <boost/any.hpp>

#include <vector>

/*
        NOTES

        the boost any stuff does have a performance impact...

        AS IS: set_scene rate: 1.26338e+06 / sec
        without new/del       ~1.75000e+06 / sec

        directly setting scene from NEST without boost::any or virtual func
   calls or multinest:
        set_scene rate: 1.86553e+07 / sec [       OK ] Director.basic_test (37
   ms)


*/

namespace scheme {
namespace kinematics {

template <class _Position, class _BigIndex = uint64_t, class _Index = uint64_t>
struct Director {
  typedef _Position Position;
  typedef _BigIndex BigIndex;
  typedef _Index Index;
  typedef SceneBase<Position, Index> Scene;

  virtual bool set_scene(BigIndex const& i, int resl, Scene& scene) const = 0;

  virtual BigIndex size(int resl) const = 0;
};

template <class _Nest>
struct NestDirector
    : public Director<typename _Nest::Value, typename _Nest::Index,
                      typename _Nest::Index> {
  typedef _Nest Nest;
  typedef typename Nest::Value Position;
  typedef typename Nest::Index BigIndex;
  typedef typename Nest::Index Index;
  typedef SceneBase<Position, Index> Scene;

  Index ibody_;
  Nest nest_;

  NestDirector() : ibody_(0), nest_() {}
  NestDirector(Index ibody) : ibody_(ibody), nest_() {}
  template <class A>
  NestDirector(A const& a, Index ibody) : ibody_(ibody), nest_(a) {}
  template <class A, class B>
  NestDirector(A const& a, B const& b, Index ibody)
      : ibody_(ibody), nest_(a, b) {}
  template <class A, class B, class C>
  NestDirector(A const& a, B const& b, C const& c, Index ibody)
      : ibody_(ibody), nest_(a, b, c) {}
  template <class A, class B, class C, class D>
  NestDirector(A const& a, B const& b, C const& c, D const& d, Index ibody)
      : ibody_(ibody), nest_(a, b, c, d) {}

  Nest const& nest() const { return nest_; }

  virtual bool set_scene(BigIndex const& i, int resl, Scene& scene) const {
    Position p;
    bool success = nest_.get_state(i, resl, p);
    if (!success) return false;
    scene.set_position(ibody_, p);
    return true;
  }

  virtual BigIndex size(int resl) const { return nest_.size(resl); }
};

template <class Nest>
std::ostream& operator<<(std::ostream& out, NestDirector<Nest> const& d) {
  out << "NestDirector " << d.nest();
  return out;
}

class SymElem {};

template <class _Position, class _Index = uint64_t>
struct SceneTree {
  typedef SceneTree<_Position, _Index> THIS;
  typedef _Position Position;
  typedef _Index Index;
  typedef SceneBase<Position, Index> Scene;
  typedef nest::NestBase<Index> Nest;
  typedef shared_ptr<Nest> NestP;
  typedef shared_ptr<SceneTree> SceneTreeP;

  Position edge_xform_;
  std::vector<NestP> position_nests_;
  std::vector<NestP> sym_elem_nests_;
  std::vector<NestP> conformation_nests_;

  std::vector<Index> bodies_;
  std::vector<SymElem> sym_elems_;

  std::vector<SceneTreeP> children_;

  // mutable Position tmp_;

  SceneTree(Position const& edge_xform = Position::Identity())
      : edge_xform_(edge_xform) {}

  void add_body(Index ibody) { bodies_.push_back(ibody); }
  void add_position_nest(NestP nestp) { position_nests_.push_back(nestp); }
  void add_child(SceneTreeP child) { children_.push_back(child); }

  void get_nests(std::vector<NestP>& nests) const {
    nests.insert(nests.end(), position_nests_.begin(), position_nests_.end());
    nests.insert(nests.end(), sym_elem_nests_.begin(), sym_elem_nests_.end());
    nests.insert(nests.end(), conformation_nests_.begin(),
                 conformation_nests_.end());
    BOOST_FOREACH (SceneTreeP child, children_) { child->get_nests(nests); }
  }

  void get_dofs(std::vector<boost::any>& dofs) const {
    for (int ipos = 0; ipos < position_nests_.size(); ++ipos) {
      // std::cerr << "create temporary Position" << std::endl;
      dofs.push_back(new Position);
      // Position *tmp = &tmp_;
      // dofs.push_back( tmp );
    }
    if (sym_elem_nests_.size())
      throw std::logic_error("sym_elem_nests not implemented");
    if (conformation_nests_.size())
      throw std::logic_error("conformation nests not implemented");
    BOOST_FOREACH (SceneTreeP child, children_) { child->get_dofs(dofs); }
  }
  void clear_empty_dofs(std::vector<boost::any>::iterator idofs) const {
    for (int ipos = 0; ipos < position_nests_.size(); ++ipos) {
      // std::cerr << "delete temporary Position" << std::endl;
      delete boost::any_cast<Position*>(*idofs);
      ++idofs;
    }
    if (sym_elem_nests_.size())
      throw std::logic_error("sym_elem_nests not implemented");
    if (conformation_nests_.size())
      throw std::logic_error("conformation nests not implemented");
    BOOST_FOREACH (SceneTreeP child, children_) {
      child->clear_empty_dofs(idofs);
    }
  }

  void set_scene(Position const& parent_position,
                 std::vector<boost::any>::const_iterator idof,
                 Scene& scene) const {
    Position position(edge_xform_ * parent_position);
    for (int ipos = 0; ipos < position_nests_.size(); ++ipos) {
      // std::cerr << "set_scene get nest xform" << std::endl;
      Position& nest_xform(*boost::any_cast<Position*>(*idof));
      // std::cerr << "set_scene get nest xform DONE" << std::endl;
      position = nest_xform * position;
      ++idof;
    }
    if (sym_elem_nests_.size())
      throw std::logic_error("sym_elem_nests not implemented");
    if (conformation_nests_.size())
      throw std::logic_error("conformation nests not implemented");
    BOOST_FOREACH (Index ibody, bodies_) {
      scene.set_position(ibody, position);
    }
    BOOST_FOREACH (SceneTreeP child, children_) {
      child->set_scene(parent_position, idof, scene);
    }
  }
};

// MultiNest(Nests const & nests)
// void add_nest( NestP nestp )

template <class _Position, class _BigIndex = uint64_t, class _Index = uint64_t>
struct TreeDirector : public Director<_Position, _BigIndex, _Index> {
  typedef _BigIndex BigIndex;
  typedef _Index Index;
  typedef _Position Position;
  // typedef SceneTree<Position,Index> SceneTree;
  typedef shared_ptr<SceneTree<Position, Index> > SceneTreeP;
  typedef SceneBase<Position, Index> Scene;
  typedef nest::NestBase<Index> Nest;
  typedef shared_ptr<Nest> NestP;
  typedef nest::MultiNest<Index, BigIndex> MultiNest;

  SceneTreeP root_;
  MultiNest multinest_;

  TreeDirector(SceneTreeP root) : root_(root) { init(); }

  void init() {
    std::vector<NestP> nests;
    root_->get_nests(nests);
    multinest_.init(nests);
  }

  virtual bool set_scene(BigIndex const& i, int resl, Scene& scene) const {
    std::vector<boost::any> tmp_dofs;
    // std::cerr << "get_dofs" << std::endl;
    root_->get_dofs(tmp_dofs);
    assert(tmp_dofs.size() == multinest_.size());
    // std::cerr << "get_states" << std::endl;
    bool success = multinest_.get_states(i, resl, tmp_dofs);
    if (!success) return false;
    // std::cerr << *boost::any_cast<numeric::X1dim*>(tmp_dofs.front()) <<
    // std::endl;
    assert(root_ != 0);
    // std::cerr << "set_scene" << std::endl;
    root_->set_scene(Position::Identity(), tmp_dofs.begin(), scene);
    // std::cerr << "clear_empty_dofs" << std::endl;
    root_->clear_empty_dofs(tmp_dofs.begin());
    return true;
  }

  virtual BigIndex size(int resl) const { return multinest_.size(resl); }

 private:
  TreeDirector() {}
};
}
}

#endif
