#ifndef INCLUDED_kinematics_SceneBase_HH
#define INCLUDED_kinematics_SceneBase_HH

#include "types.hpp"
#include "util/meta/util.hpp"

#include <boost/any.hpp>
#include <vector>

namespace scheme {
namespace kinematics {

template <class _Position, class _Index = uint64_t>
struct SceneBase {
  typedef _Index Index;
  typedef _Position Position;
  typedef SceneBase<Position, Index> THIS;

  // members
  std::vector<Position> symframes_;  // include identity at first position
  Index n_sym_bodies_, n_bodies_;
  std::vector<Position> positions_;

  SceneBase()
      : n_bodies_(0), n_sym_bodies_(0), symframes_(1, Position::Identity()) {}

  virtual ~SceneBase() {}

  virtual shared_ptr<THIS> clone_shallow() const {
    return clone_deep();
  }  // default to deep copy
  virtual shared_ptr<THIS> clone_deep() const = 0;

  Position position(Index i) const {
    Index isym = this->sym_index_map(i);
    return this->symframes_.at(isym) * positions_.at(i);
  }
  void set_position(Index i, Position const &newp) { positions_.at(i) = newp; }
  template <class Actor>
  bool add_actor(Index ib, Actor const &actor) {
    boost::any any_actor = actor;
    return add_actor(ib, any_actor);
  }
  virtual bool add_actor(Index ib, boost::any const &a) { return false; }

  template <class Actor>
  Actor get_actor(Index ib, Index ia) const {
    boost::any any_actor = util::meta::type2type<Actor>();
    get_actor(ib, ia, any_actor);
    return boost::any_cast<Actor>(any_actor);
  }
  virtual bool get_actor(Index ib, Index ia, boost::any &a) const {
    return false;
  }

  template <class Actor>
  Actor &get_nonconst_actor(Index ib, Index ia) {
    boost::any any_actor = util::meta::type2type<Actor>();
    get_nonconst_actor(ib, ia, any_actor);
    return *boost::any_cast<Actor *>(any_actor);
  }
  virtual bool get_nonconst_actor(Index ib, Index ia, boost::any &a) {
    return false;
  }

  template <class Actor>
  int num_actors(Index ib) const {
    boost::any any_actor = util::meta::type2type<Actor>();
    return num_actors(ib, any_actor);
  }
  virtual int num_actors(Index ib, boost::any const &a) const { return -1; }

  template <class Actor>
  bool clear_actors(Index ib) {
    Actor a;
    boost::any any_actor = a;
    return clear_actors(ib, any_actor);
  }
  virtual bool clear_actors(Index ib, boost::any const &a) { return false; }

  // symmetry stuff, doesn't need to be Conformation-sepcific
  void update_symmetry(Index nbodies) {
    n_bodies_ = nbodies;
    n_sym_bodies_ = n_bodies_ * ((Index)symframes_.size());
  }
  Index sym_index_map(Index &i) const {
    Index isym = i / n_bodies_;
    i = i % n_bodies_;
    return isym;
  }

  void set_symmetry(std::vector<Position> const &sym) {
    symframes_ = sym;
    update_symmetry(n_bodies_);
  }
  void add_symframe(Position const &symframe) {
    symframes_.push_back(symframe);
    update_symmetry(n_bodies_);
  }
  std::vector<Position> const &symframes() const { return symframes_; }
  Index nbodies() const { return n_sym_bodies_; }
  Index nbodies_asym() const { return n_bodies_; }
};
}
}

#endif
