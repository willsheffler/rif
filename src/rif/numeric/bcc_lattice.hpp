#pragma once

#include "util/SimpleArray.hpp"
#include "util/template_loop.hpp"

#include <iostream>

namespace rif {
namespace numeric {

/**
 * @brief      n-dimensional BCC lattice
 *
 * @tparam     DIM     { description }
 * @tparam     _Float  { description }
 * @detail     will hit lower bound but not upper... for example
 *      if range is 0..1 and 2 cells, x values could be 0, 0.25, 0.5 and 0.75
 * ((not 1.0!!))
 */
template <int DIM, class _Float = double, class _Index = uint64_t>
struct BCC {
  using Float = _Float;
  using Index = _Index;
  using Indices = util::SimpleArray<DIM, Index>;
  using Floats = util::SimpleArray<DIM, Float>;
  BOOST_STATIC_ASSERT((DIM > 2));

  Indices nside_, nside_prefsum_;
  Floats lower_, width_, lower_cen_, half_width_;

  bool operator!=(BCC<DIM, Float, Index> o) const { return !(*this == o); }
  bool operator==(BCC<DIM, Float, Index> o) const {
    return nside_ == o.nside_ && lower_ == o.lower_ && width_ == o.width_;
  }

  BCC() {}

  template <class Sizes>
  BCC(Sizes sizes, Floats lower = Floats(0), Floats upper = Floats(1)) {
    init(sizes, lower, upper);
  }

  template <class Sizes>
  void init(Sizes sizes, Floats lower, Floats upper) {
    nside_ = sizes;
    lower_ = lower;
    for (size_t i = 0; i < DIM; ++i) nside_prefsum_[i] = nside_.prod(i);
    width_ = (upper - lower_) / nside_.template cast<Float>();
    half_width_ = width_ / 2.0;
    lower_cen_ = lower_ + half_width_;
    uint64_t totsize = 2;
    for (size_t i = 0; i < DIM; ++i) {
      if (totsize * nside_[i] < totsize)
        throw std::invalid_argument("Index Type is too narrow");
      totsize *= nside_[i];
    }
  }

  Index size() const noexcept { return nside_.prod() * 2; }

  Floats operator[](Index index) const noexcept {
    bool odd = index & 1;
    Indices indices = ((index >> 1) / nside_prefsum_) % nside_;
    return this->get_center(indices, odd);
  }

  Floats get_center(Indices indices, bool odd) const noexcept {
    return lower_cen_ + width_ * indices.template cast<Float>() +
           (odd ? half_width_ : 0);
  }

  Indices get_indices(Floats value, bool &odd) const noexcept {
    value = (value - lower_) / width_;
    Indices const indices = value.template cast<Index>();
    value = value - indices.template cast<Float>() - 0.5;
    Indices const corner_indices = indices - (value < 0).template cast<Index>();
    odd = (0.25 * DIM) < fabs((value.sign() * value).sum());
    return odd ? corner_indices : indices;
  }

  Index operator[](Floats value) const noexcept {
    bool odd;
    Indices indices = get_indices(value, odd);
    Index index = (nside_prefsum_ * indices).sum();
    return (index << 1) + odd;
  }

  template <class Iiter>
  void neighbors(Index index, Iiter iter, bool edges = false,
                 bool edges2 = false) const noexcept {
    *iter++ = index;
    bool odd = index & 1;
    Indices indices = ((index >> 1) / nside_prefsum_) % nside_;
    // std::cout << indices << std::endl;
    for (Index i = 0; i < DIM; ++i) {
      indices[i] += 1;
      // std::cout << indices << " " << i1 << std::endl;
      if ((indices < nside_).sum() == DIM)
        *iter++ = (nside_prefsum_ * indices).sum() << 1 | odd;
      indices[i] -= 2;
      // std::cout << indices << " " << i2 << std::endl;
      if ((indices < nside_).sum() == DIM)
        *iter++ = (nside_prefsum_ * indices).sum() << 1 | odd;
      indices[i] += 1;  // reset
    }
    odd = !odd;
    Index sodd = odd ? -1 : 1;
    for (Index i = 0; i < (1 << DIM); ++i) {
      Indices corner(indices);
      for (int d = 0; d < DIM; ++d) corner[d] += ((i >> d) & 1) ? sodd : 0;
      // std::cout << corner << std::endl;
      if ((corner < nside_).sum() == DIM)
        *iter++ = (nside_prefsum_ * corner).sum() << 1 | odd;
    }
    if (edges) {
      odd = !odd;
      for (Index i = 0; i < DIM - 1; ++i) {
        for (Index j = i + 1; j < DIM; ++j) {
          indices[i] += 1;
          indices[j] += 1;  // +1,+1
          // std::cout << indices << " " << i1 << std::endl;
          if ((indices < nside_).sum() == DIM)
            *iter++ = (nside_prefsum_ * indices).sum() << 1 | odd;
          indices[i] -= 2;  // -1,+1
          // std::cout << indices << " " << i2 << std::endl;
          if ((indices < nside_).sum() == DIM)
            *iter++ = (nside_prefsum_ * indices).sum() << 1 | odd;
          indices[j] -= 2;  // -1,-1
          // std::cout << indices << " " << i2 << std::endl;
          if ((indices < nside_).sum() == DIM)
            *iter++ = (nside_prefsum_ * indices).sum() << 1 | odd;
          indices[i] += 2;  // +1,-1
          // std::cout << indices << " " << i2 << std::endl;
          if ((indices < nside_).sum() == DIM)
            *iter++ = (nside_prefsum_ * indices).sum() << 1 | odd;
          // reset
          indices[i] -= 1;
          indices[j] += 1;
        }
      }
    }
  }
};

template <int DIM, class Float, class Index>
std::ostream &operator<<(std::ostream &out, BCC<DIM, Float, Index> bcc) {
  return out << "lb " << bcc.lower_ << " w " << bcc.width_;
}

template <int DIM, class Float, class Index = uint64_t>
struct Cubic {
  typedef util::SimpleArray<DIM, Index> Indices;
  typedef util::SimpleArray<DIM, Float> Floats;
  BOOST_STATIC_ASSERT((DIM > 2));

  Indices nside_, nside_prefsum_;
  Floats lower_, width_, lower_cen_, half_width_;

  Cubic() {}

  template <class Sizes>
  Cubic(Sizes sizes, Floats lower = Floats(0), Floats upper = Floats(1)) {
    init(sizes, lower, upper);
  }

  template <class Sizes>
  void init(Sizes sizes, Floats lower = Floats(0), Floats upper = Floats(1)) {
    nside_ = sizes;
    lower_ = lower;
    for (size_t i = 0; i < DIM; ++i) nside_prefsum_[i] = nside_.prod(i);
    width_ = (upper - lower_) / nside_.template cast<Float>();
    half_width_ = width_ / 2.0;
    lower_cen_ = lower_ + half_width_;
  }

  Index size() const noexcept { return nside_.prod(); }

  Floats operator[](Index index) const noexcept {
    Indices indices = (index / nside_prefsum_) % nside_;
    return get_center(indices);
  }

  Floats get_center(Indices indices) const noexcept {
    return lower_cen_ + width_ * indices.template cast<Float>();
  }

  Indices get_indices(Floats value) const noexcept {
    value = (value - lower_) / width_;
    return value.template cast<Index>();
  }

  Index operator[](Floats value) const noexcept {
    Indices indices = get_indices(value);
    return (nside_prefsum_ * indices).sum();
  }

  template <class Iiter>
  void neighbors(Index index, Iiter iter, bool = false) const noexcept {
    Indices idx0 = (index / nside_prefsum_) % nside_;
    Indices threes(1);
    for (int d = 1; d < DIM; ++d) threes[d] = 3 * threes[d - 1];
    for (int i = 0; i < threes[DIM - 1] * 3; ++i) {
      Indices idx(idx0);
      for (int d = 0; d < DIM; ++d) idx[d] += ((i / threes[d]) % 3) - 1;
      // std::cout << i << " " << (idx-idx0).template cast<int>()+1 <<
      // std::endl;
      if ((idx < nside_).sum() == DIM) *iter++ = (nside_prefsum_ * idx).sum();
    }
  }
};
}
}
