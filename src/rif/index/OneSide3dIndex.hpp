// gross but efficient code from c/gpu stuff many years ago
// todo: refactor wtih visitor to remove extra code

#pragma once

#include "rif/eigen_types.hpp"
#include "rif/util/assert.hpp"

namespace rif {
namespace index {

template <class Pt>
Pt &getPt(Pt &p) {
  return p;
}

template <class Pt, class Payload>
Pt &getPt(std::pair<Pt, Payload> &p) {
  return p.firt;
}

template <class Pt>
void copyPayload(Pt, Pt) {}

template <class Pt, class Payload>
void copyPayload(std::pair<Pt, Payload> &a, std::pair<Pt, Payload> &b) {
  a.second = b.second;
}

/// iff Payload is void, type stored will be Pt, else std::pair<Pt, Payload>
template <class _Pt, class _Payload = void, class _ShortIdx = uint16_t>
struct OneSide3dIndex {
  using Pt = _Pt;
  using Payload = _Payload;
  using ShortIdx = _ShortIdx;
  using F = typename Pt::Scalar;
  using Value = typename std::conditional<std::is_same<Payload, void>::value,
                                          Pt, std::pair<Pt, Payload>>::type;
  using ShortIdxPair = std::pair<ShortIdx, ShortIdx>;

  F width_, width2_;  // cell width
  size_t Npts;
  Value const *values_ = nullptr;
  ShortIdxPair const *pindex_ = nullptr;
  int xdim_, ydim_, zdim_;
  float xmx_, ymx_, zmx_;
  V3f translation_;

  OneSide3dIndex() {}

  template <class C>
  OneSide3dIndex(C const &pts, F width) {
    Value const *ptr = &pts[0];
    init(width, ptr, pts.size());
  }

  OneSide3dIndex(F width, Value const *ptr, size_t Npts) {
    init(width, ptr, Npts);
  }

  /// basic visitor to count contacts, good one to take as an example
  struct CountVisitor {
    int result = 0;
    bool visit(Pt, Value, float) {
      ++result;
      return false;  // no early termination, want full count
    }
  };

  struct ContactVisitor {
    bool result = false;
    bool visit(Pt, Pt, float) {
      result = true;
      return true;  // return true to terminate early
    }
  };

  template <class Point>
  int nbcount(Point query) const {
    CountVisitor myvisitor;
    visit(query, myvisitor);
    return myvisitor.result;
  }

  template <class Point>
  bool contact(Point query) const {
    ContactVisitor myvisitor;
    visit(query, myvisitor);
    return myvisitor.result;
  }

  template <class Point>
  int brute_nbcount(Point query) const {
    CountVisitor myvisitor;
    brute_visit(query, myvisitor);
    return myvisitor.result;
  }

  template <class Point>
  int brute_contact(Point query) const {
    ContactVisitor myvisitor;
    brute_visit(query, myvisitor);
    return myvisitor.result;
  }

  void init(F width, Value const *ptr, size_t Npts);

  virtual ~OneSide3dIndex() {
    if (values_) delete values_;
    if (pindex_) delete pindex_;
  }

  bool sanity_check() const {
    for (int ix = 0; ix < xdim_; ++ix) {
      for (int iy = 0; iy < ydim_; ++iy) {
        for (int iz = 0; iz < zdim_; ++iz) {
          // std::cout << ix << " " << iy << " " << iz << endl;
          ushort const ig = ix + xdim_ * iy + ydim_ * xdim_ * iz;
          ushort const igl = pindex_[ig].first;
          ushort const igu = pindex_[ig].second;
          for (int i = igl; i < igu; ++i) {
            // float const & x(values_[i].first);
            float const &y(getPt(values_[i])[1]);
            float const &z(getPt(values_[i])[2]);
            // if(i==igl) std::cout << endl;
            // bool xc = width_*(float)ix <= x && x <=
            // width_*(float)(ix+1);
            bool yc = width_ * (float)iy <= y && y <= width_ * (float)(iy + 1);
            bool zc = width_ * (float)iz <= z && z <= width_ * (float)(iz + 1);
            if (/*!xc||*/ !yc || !zc) {
              ALWAYS_ASSERT_MSG(false,
                                "insanity in OneSide3dIndex::sanity_check");
            }
          }
        }
        return true;
      }
    }
    return true;
  }

  template <class OrigVisitor>
  struct DistSqCheckVisitor {
    float mywidth2;
    OrigVisitor &visitor;
    DistSqCheckVisitor(float w, OrigVisitor &v) : mywidth2(w), visitor(v) {}
    bool visit(Pt v, Pt c) {
      float d2 = (v[0] - c[0]) * (v[0] - c[0]) + (v[1] - c[1]) * (v[1] - c[1]) +
                 (v[2] - c[2]) * (v[2] - c[2]);
      if (d2 <= mywidth2) return visitor.visit(v, c, d2);
      return false;
    }
  };

  template <typename Visitor, typename Point>
  void visit(Point const &query, Visitor &visitor) const {
    DistSqCheckVisitor<Visitor> d2_checking_visitor(width2_, visitor);
    visit_lax(query, d2_checking_visitor);
  }

  template <typename Visitor, typename Point>
  void visit_lax(Point query_in, Visitor &visitor) const {
    Pt query(query_in);
    query[0] += translation_[0];
    query[1] += translation_[1];
    query[2] += translation_[2];
    float x = query[0];
    float y = query[1];
    float z = query[2];
    if (x < -width_ || y < -width_ || z < -width_) return;  // worth it iff
    if (x > xmx_ || y > ymx_ || z > zmx_) return;           // worth it iff
    int const ix = (x < 0) ? 0 : std::min(xdim_ - 1, (int)(x / width_));
    int const iy0 = (y < 0) ? 0 : (int)(y / width_);
    int const iz0 = (z < 0) ? 0 : (int)(z / width_);
    int const iyl = std::max(0, iy0 - 1);
    int const izl = std::max(0, iz0 - 1);
    int const iyu = std::min((int)(ydim_), iy0 + 2);
    int const izu = std::min((int)(zdim_), (int)(iz0 + 2));
    for (int iy = iyl; iy < iyu; ++iy) {
      for (int iz = izl; iz < izu; ++iz) {
        int const ig = ix + xdim_ * iy + xdim_ * ydim_ * iz;
        assert(ig < xdim_ * ydim_ * zdim_);
        assert(ix < xdim_);
        assert(iy < ydim_);
        assert(iz < zdim_);
        int const &igl = pindex_[ig].first;
        int const &igu = pindex_[ig].second;
        for (int i = igl; i < igu; ++i) {
          if (visitor.visit(query, values_[i])) return;
        }
      }
    }
  }

  template <typename Visitor, typename Point>
  void brute_visit(Point const &query, Visitor &visitor) const {
    DistSqCheckVisitor<Visitor> d2_checking_visitor(width2_, visitor);
    brute_visit_lax(query, d2_checking_visitor);
  }

  /// testing utility function
  template <typename Visitor, typename Point>
  void brute_visit_lax(Point query_in, Visitor &visitor) const {
    Pt query(query_in);
    query[0] += translation_[0];
    query[1] += translation_[1];
    query[2] += translation_[2];
    for (size_t i = 0; i < Npts; ++i) {
      visitor.visit(query, values_[i]);
    }
  }
};

template <class Pt, class Payload, class ShortIdx>
void OneSide3dIndex<Pt, Payload, ShortIdx>::init(F width, Value const *ptr,
                                                 size_t _Npts) {
  Npts = _Npts;
  assert(Npts);

  width_ = width;
  width2_ = width * width;

  F xmn = 9e9, ymn = 9e9, zmn = 9e9;
  F xmx = -9e9, ymx = -9e9, zmx = -9e9;
  for (int i = 0; i < Npts; ++i) {
    xmn = std::min(xmn, getPt(ptr[i])[0]);
    ymn = std::min(ymn, getPt(ptr[i])[1]);
    zmn = std::min(zmn, getPt(ptr[i])[2]);
    xmx = std::max(xmx, getPt(ptr[i])[0]);
    ymx = std::max(ymx, getPt(ptr[i])[1]);
    zmx = std::max(zmx, getPt(ptr[i])[2]);
  }

  xdim_ = (int)((xmx - xmn + 0.0001) / width_ + 0.999999);
  ydim_ = (int)((ymx - ymn + 0.0001) / width_ + 0.999999);
  zdim_ = (int)((zmx - zmn + 0.0001) / width_ + 0.999999);
  assert(xdim_ < 9999);
  assert(ydim_ < 9999);
  assert(zdim_ < 9999);
  int const gsize = xdim_ * ydim_ * zdim_;
  ShortIdxPair *tmp_pindex = new ShortIdxPair[gsize];
  ShortIdxPair *tmp_pindex2 = new ShortIdxPair[gsize];

  for (int i = 0; i < gsize; ++i) {
    tmp_pindex2[i].first = 0;
    tmp_pindex2[i].second = 0;
  }
  // TR<<"atom "<<Npts<<" grid1 "<<xdim_*ydim_*zdim_<<" "<<xdim_<<"
  // "<<ydim_<<" "<<zdim_<<std::endl;

  for (int i = 0; i < Npts; ++i) {
    int ix = (int)((getPt(ptr[i])[0] - xmn /*+FUDGE*/) / width_);
    int iy = (int)((getPt(ptr[i])[1] - ymn /*+FUDGE*/) / width_);
    int iz = (int)((getPt(ptr[i])[2] - zmn /*+FUDGE*/) / width_);
    assert(ix >= 0);
    assert(iy >= 0);
    assert(iz >= 0);
    assert(ix < xdim_);
    assert(iy < ydim_);
    assert(iz < zdim_);
    int ig = ix + xdim_ * iy + xdim_ * ydim_ * iz;
    assert(ig >= 0);
    assert(ig < 9999999);
    ++(tmp_pindex2[ig].second);
  }
  for (int i = 1; i < gsize; ++i)
    tmp_pindex2[i].first = tmp_pindex2[i - 1].first + tmp_pindex2[i - 1].second;
  for (int i = 1; i < gsize; ++i)
    tmp_pindex2[i].second = tmp_pindex2[i].first + tmp_pindex2[i].second;
  for (int iz = 0; iz < zdim_; ++iz)
    for (int iy = 0; iy < ydim_; ++iy)
      for (int ix = 0; ix < xdim_; ++ix) {
        int const ixl = (int)std::max(0, (int)ix - 1);
        int const ixu = std::min(xdim_ - 1u, ix + 1u);
        int const ig0 = xdim_ * iy + xdim_ * ydim_ * iz;
        tmp_pindex[ix + ig0].first = tmp_pindex2[ixl + ig0].first;
        tmp_pindex[ix + ig0].second = tmp_pindex2[ixu + ig0].second;
      }
  // for(int iz = 0; iz < zdim_; ++iz) for(int iy = 0; iy < ydim_; ++iy)
  // for(int ix = 0; ix < xdim_; ++ix) {
  //       int i = ix+xdim_*iy+xdim_*ydim_*iz;
  //       TR<<ix<<" "<<iy<<" "<<iz<<" "<<I(3,tmp_pindex2[i].first)<<"
  //       "<<I(3,tmp_pindex2[i].second) <<" "<<I(3,tmp_pindex[i].first)<<"
  //       "<<I(3,tmp_pindex[i].second)<<std::endl;
  //     }
  pindex_ = tmp_pindex;
  Value *values = new Value[Npts + 4];  // space for 4 overflow ptr
  for (int i = 0; i < 4; ++i) {
    getPt(values[Npts + i])[0] = 9e9;
    getPt(values[Npts + i])[1] = 9e9;
    getPt(values[Npts + i])[2] = 9e9;
  }
  ushort *gridc = new ushort[gsize];
  for (int i = 0; i < gsize; ++i) gridc[i] = 0;
  for (int i = 0; i < Npts; ++i) {
    int const ix = (int)((getPt(ptr[i])[0] - xmn /*+FUDGE*/) / width_);
    int const iy = (int)((getPt(ptr[i])[1] - ymn /*+FUDGE*/) / width_);
    int const iz = (int)((getPt(ptr[i])[2] - zmn /*+FUDGE*/) / width_);
    int const ig = ix + xdim_ * iy + xdim_ * ydim_ * iz;
    int const idx = tmp_pindex2[ig].first + gridc[ig];
    getPt(values[idx])[0] = getPt(ptr[i])[0] - xmn /*+FUDGE*/;
    getPt(values[idx])[1] = getPt(ptr[i])[1] - ymn /*+FUDGE*/;
    getPt(values[idx])[2] = getPt(ptr[i])[2] - zmn /*+FUDGE*/;
    copyPayload(values[idx], ptr[i]);
    ++(gridc[ig]);
  }
  values_ = values;
  translation_[0] = /* FUDGE*/ -xmn;
  translation_[1] = /* FUDGE*/ -ymn;
  translation_[2] = /* FUDGE*/ -zmn;
  xmx_ = xmx - xmn /*+FUDGE*/ + width_;
  ymx_ = ymx - ymn /*+FUDGE*/ + width_;
  zmx_ = zmx - zmn /*+FUDGE*/ + width_;

  delete gridc;
  delete tmp_pindex2;
}
}  // end index ns
}  // end rif ns
