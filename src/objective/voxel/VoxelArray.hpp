#ifndef INCLUDED_objective_voxel_VoxelArray_HH
#define INCLUDED_objective_voxel_VoxelArray_HH

#include <boost/assert.hpp>
#include <boost/multi_array.hpp>
#include <boost/type_traits.hpp>
#include "util/SimpleArray.hpp"
#include "util/assert.hpp"

#ifdef CEREAL
#include <cereal/access.hpp>
#endif
//#include <boost/serialization/access.hpp>
//#include <boost/serialization/split_member.hpp>

namespace scheme {
namespace objective {
namespace voxel {

template <size_t _DIM, class _Float = float>
struct VoxelArray : boost::multi_array<_Float, _DIM> {
  BOOST_STATIC_ASSERT((_DIM > 0));
  typedef boost::multi_array<_Float, _DIM> BASE;
  typedef VoxelArray<_DIM, _Float> THIS;
  static size_t const DIM = _DIM;
  typedef _Float Float;
  typedef util::SimpleArray<DIM, typename BASE::size_type> Indices;
  typedef util::SimpleArray<DIM, Float> Bounds;
  Bounds lb_, ub_, cs_;

  VoxelArray() {}

  template <class F1, class F2, class F3>
  VoxelArray(F1 const &lb, F2 const &ub, F3 const &cs)
      : lb_(lb), ub_(ub), cs_(cs) {
    Indices extents = floats_to_index(ub_);
    // std::cout << extents << std::endl;
    this->resize(extents + Indices(1));  // pad by one
  }

  template <class Floats>
  Indices floats_to_index(Floats const &f) const {
    Indices ind;
    for (int i = 0; i < DIM; ++i) {
      Float tmp = ((f[i] - lb_[i]) / cs_[i]);
      // assert(tmp >= 0.0);
      ind[i] = tmp;
    }
    // std::cout << "floats_to_index " << f << " " << ind << std::endl;
    return ind;
  }

  Bounds indices_to_center(Indices const &idx) const {
    Bounds c;
    for (int i = 0; i < DIM; ++i) {
      c[i] = (idx[i] + 0.5) * cs_[i] + lb_[i];
    }
    return c;
  }

  template <class Floats>
  typename boost::disable_if<boost::is_arithmetic<Floats>, Float const &>::type
  operator[](Floats const &floats) const {
    return this->operator()(floats_to_index(floats));
  }

  template <class Floats>
  typename boost::disable_if<boost::is_arithmetic<Floats>, Float &>::type
  operator[](Floats const &floats) {
    return this->operator()(floats_to_index(floats));
  }

  Float at(Float f, Float g, Float h) const {
    Indices idx = floats_to_index(Bounds(f, g, h));
    if (idx[0] < this->shape()[0] && idx[1] < this->shape()[1] &&
        idx[2] < this->shape()[2])
      return this->operator()(idx);
    else
      return 0.0;
  }

  template <class V>
  Float at(V const &v) const {
    Indices idx = floats_to_index(Bounds(v[0], v[1], v[2]));
    if (idx[0] < this->shape()[0] && idx[1] < this->shape()[1] &&
        idx[2] < this->shape()[2])
      return this->operator()(idx);
    else
      return 0.0;
  }

  // void write(std::ostream & out) const {
  // 	out.write( (char const*)&lb_, sizeof(Bounds) );
  // 	out.write( (char const*)&ub_, sizeof(Bounds) );
  // 	out.write( (char const*)&cs_, sizeof(Bounds) );
  // 	for(size_t i = 0; i < DIM; ++i) out.write( (char
  // const*)&(this->shape()[i]), sizeof() );
  // 	out.write( (char const*)this->data(), this->num_elements()*sizeof(Float)
  // );
  // }
  // void read(std::istream & in){
  // 	in.read( (char*)&lb_, sizeof(Bounds) );
  // 	in.read( (char*)&ub_, sizeof(Bounds) );
  // 	in.read( (char*)&cs_, sizeof(Bounds) );
  // 	in.read( (char*)this->data(), this->num_elements()*sizeof(Float) );
  // }
  bool operator==(THIS const &o) const {
    return lb_ == o.lb_ && ub_ == o.ub_ && cs_ == o.cs_ &&
           (BASE const &)o == (BASE const &)*this;
  }

#ifdef CEREAL
  // friend class boost::serialization::access;
  friend class cereal::access;  // befriend the cereal version of access
#endif

  template <class Archive>
  void save(Archive &ar, const unsigned int) const {
    BOOST_VERIFY(boost::is_pod<Float>::type::value);
    ar &lb_;
    ar &ub_;
    ar &cs_;
    for (size_t i = 0; i < DIM; ++i) ar & this->shape()[i];
    for (size_t i = 0; i < this->num_elements(); ++i) ar & this->data()[i];
  }
  template <class Archive>
  void load(Archive &ar, const unsigned int) {
    BOOST_VERIFY(boost::is_pod<Float>::type::value);
    ar &lb_;
    ar &ub_;
    ar &cs_;
    Indices extents;
    for (size_t i = 0; i < DIM; ++i) ar &extents[i];
    this->resize(extents);
    for (size_t i = 0; i < this->num_elements(); ++i) ar & this->data()[i];
  }
  void save(std::ostream &out) const {
    BOOST_VERIFY(boost::is_pod<Float>::type::value);
    out.write((char *)&lb_, sizeof(Bounds));
    out.write((char *)&ub_, sizeof(Bounds));
    out.write((char *)&cs_, sizeof(Bounds));
    for (size_t i = 0; i < DIM; ++i) {
      out.write((char *)(&(this->shape()[i])),
                sizeof(typename BASE::size_type));
    }
    for (size_t i = 0; i < this->num_elements(); ++i)
      out.write((char *)(&(this->data()[i])), sizeof(Float));
  }
  void load(std::istream &in) {
    BOOST_VERIFY(boost::is_pod<Float>::type::value);
    ALWAYS_ASSERT(in.good());
    in.read((char *)&lb_, sizeof(Bounds));
    ALWAYS_ASSERT(in.good());
    in.read((char *)&ub_, sizeof(Bounds));
    ALWAYS_ASSERT(in.good());
    in.read((char *)&cs_, sizeof(Bounds));
    ALWAYS_ASSERT(in.good());
    // todo: should I check these against the c'tor values? if not, should add
    // default ctor?
    Indices extents;
    for (size_t i = 0; i < DIM; ++i) {
      in.read((char *)(&(extents[i])), sizeof(typename BASE::size_type));
    }
    ALWAYS_ASSERT(in.good());
    this->resize(extents);
    for (size_t i = 0; i < this->num_elements(); ++i)
      in.read((char *)(&(this->data()[i])), sizeof(Float));
    ALWAYS_ASSERT(in.good());
  }
  // BOOST_SERIALIZATION_SPLIT_MEMBER()
};

template <size_t D, class F>
std::ostream &operator<<(std::ostream &out, VoxelArray<D, F> const &v) {
  out << "VoxelArray( lb: " << v.lb_ << " ub: " << v.ub_ << " cs: " << v.cs_
      << " nelem: " << v.num_elements() << " sizeof_val: " << sizeof(F) << " )";
  return out;
}
}
}
}

#endif
