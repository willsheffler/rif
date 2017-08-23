#pragma once

#include "rif/eigen_types.hpp"
#include "rif/util/assert.hpp"
#include "rif/util/maybe_pair.hpp"

namespace rif {
namespace index {

using namespace rif::util;

/// iff Payload is void, type stored will be Pt, else std::pair<Pt, Payload>
template <class _Pt, class _Payload = void, class _ShortIdx = uint16_t>
struct stripe_index_3d {
  using This = stripe_index_3d<_Pt, _Payload, _ShortIdx>;
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

  stripe_index_3d() {}
  template <class C>
  stripe_index_3d(F width, C const &pts) {
    Value const *ptsp = &pts[0];
    init(width, ptsp, pts.size());
  }
  template <class C, class D>
  stripe_index_3d(F width, C const &pts, D const &vals) {
    Pt const *ptsp = &pts[0];
    Payload const *valsp = &vals[0];
    init(width, ptsp, pts.size(), valsp);
  }
  virtual ~stripe_index_3d() {
    if (values_) delete values_;
    if (pindex_) delete pindex_;
  }
  void init(F width, Pt const *ptr, size_t Npts, Payload const *ptr2 = nullptr);

  size_t size() const { return this->Npts; }

  template <class Point>
  Pt prepare_query_point(Point q_in) const {
    Pt q;
    for (int i = 0; i < 3; ++i) q[i] = q_in[i] + translation_[i];
    return q;
  }
  Pt prepare_result_point(Pt r_in) const {
    Pt r;
    for (int i = 0; i < 3; ++i) r[i] = r_in[i] - translation_[i];
    return r;
  }

  /// basic visitor to count contacts, good one to take as an example
  struct CountVisitor {
    int result = 0;
    bool visit(Pt, Value, float) {
      ++result;
      return false;  // no early termination, want full count
    }
  };
  struct NBExistsVisitor {
    bool result = false;
    bool visit(Pt, Value, float) {
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
  bool nbexists(Point query) const {
    NBExistsVisitor myvisitor;
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
  int brute_nbexists(Point query) const {
    NBExistsVisitor myvisitor;
    brute_visit(query, myvisitor);
    return myvisitor.result;
  }

  //////////// neighbor values lookup ///////////

  template <class OutputIt>
  struct NeighboringValuesVisitor {
    OutputIt outp;
    This const *parent;
    NeighboringValuesVisitor(This const *parent, OutputIt outp)
        : parent(parent), outp(outp) {}
    bool visit(Pt, Value value, float) {
      Value nbr = value;
      get_first_if_pair(nbr) =
          parent->prepare_result_point(get_first_if_pair(nbr));
      *outp++ = nbr;
      return false;
    }
  };
  template <class OutputIt>
  NeighboringValuesVisitor<OutputIt> neighboring_values_visitor(
      OutputIt out) const {
    return NeighboringValuesVisitor<OutputIt>(this, out);
  }
  template <class Point>
  std::vector<Value> neighboring_values(Point query) const {
    std::vector<Value> result;
    auto myvisitor = neighboring_values_visitor(std::back_inserter(result));
    visit(query, myvisitor);
    return result;
  }

  //////////// neighbor points lookup ///////////

  template <class OutputIt>
  struct NeighboringPointsVisitor {
    OutputIt outp;
    This const *parent;
    NeighboringPointsVisitor(This const *parent, OutputIt outp)
        : parent(parent), outp(outp) {}
    bool visit(Pt, Value value, float) {
      Pt nbr = get_first_if_pair(value);
      *outp++ = parent->prepare_result_point(nbr);
      return false;
    }
  };
  template <class OutputIt>
  NeighboringPointsVisitor<OutputIt> neighboring_points_visitor(
      OutputIt out) const {
    return NeighboringPointsVisitor<OutputIt>(this, out);
  }
  template <class Point>
  std::vector<Pt> neighboring_points(Point query) const {
    std::vector<Point> result;
    auto myvisitor = neighboring_points_visitor(std::back_inserter(result));
    visit(query, myvisitor);
    return result;
  }
  template <class Point>
  std::vector<Pt> neighboring_points_brute(Point query) const {
    std::vector<Point> result;
    auto myvisitor = neighboring_points_visitor(std::back_inserter(result));
    brute_visit(query, myvisitor);
    return result;
  }

  //////////// neighbor payloads lookup ///////////

  template <class OutputIt>
  struct NeighboringPayloadsVisitor {
    OutputIt outp;
    NeighboringPayloadsVisitor(OutputIt outp) : outp(outp) {}
    bool visit(Pt, Value value, float) {
      *outp++ = get_second_if_pair(value);
      return false;
    }
  };
  template <class OutputIt>
  NeighboringPayloadsVisitor<OutputIt> neighboring_payloads_visitor(
      OutputIt out) const {
    return NeighboringPayloadsVisitor<OutputIt>(out);
  }
  template <class Point>
  std::vector<Payload> neighboring_payloads(Point query) const {
    std::vector<Payload> result;
    auto myvisitor = neighboring_payloads_visitor(std::back_inserter(result));
    visit(query, myvisitor);
    return result;
  }
  template <class Point>
  std::vector<Payload> neighboring_payloads_brute(Point query) const {
    std::vector<Payload> result;
    auto myvisitor = neighboring_payloads_visitor(std::back_inserter(result));
    brute_visit(query, myvisitor);
    return result;
  }

  template <class OrigVisitor>
  struct DistSqCheckVisitor {
    float mywidth2;
    OrigVisitor &visitor;
    DistSqCheckVisitor(float w, OrigVisitor &v) : mywidth2(w), visitor(v) {}
    bool visit(Pt q, Value c_in) {
      Pt c = get_first_if_pair(c_in);
      float d2 = (q[0] - c[0]) * (q[0] - c[0]) + (q[1] - c[1]) * (q[1] - c[1]) +
                 (q[2] - c[2]) * (q[2] - c[2]);
      if (d2 <= mywidth2) return visitor.visit(q, c_in, d2);
      return false;
    }
  };
  template <class OrigVisitor>
  DistSqCheckVisitor<OrigVisitor> distsq_check_visitor(OrigVisitor &v) const {
    return DistSqCheckVisitor<OrigVisitor>(width2_, v);
  }

  template <typename Visitor, typename Point>
  void visit(Point const &query, Visitor &visitor) const {
    auto d2_checking_visitor = distsq_check_visitor(visitor);
    visit_lax(query, d2_checking_visitor);
  }
  template <typename Visitor, typename Point>
  void brute_visit(Point const &query, Visitor &visitor) const {
    DistSqCheckVisitor<Visitor> d2_checking_visitor(width2_, visitor);
    brute_visit_lax(query, d2_checking_visitor);
  }

  // finally, the actual guts of the thing... this and init
  template <typename Visitor, typename Point>
  void visit_lax(Point query_in, Visitor &visitor) const {
    Pt query = prepare_query_point(query_in);
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
  void brute_visit_lax(Point query_in, Visitor &visitor) const {
    Pt query = prepare_query_point(query_in);
    for (size_t i = 0; i < Npts; ++i) {
      visitor.visit(query, values_[i]);
    }
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
            Pt pt = get_first_if_pair(values_[i]);
            float const &y(pt[1]);
            float const &z(pt[2]);
            // if(i==igl) std::cout << endl;
            // bool xc = width_*(float)ix <= x && x <=
            // width_*(float)(ix+1);
            bool yc = width_ * (float)iy <= y && y <= width_ * (float)(iy + 1);
            bool zc = width_ * (float)iz <= z && z <= width_ * (float)(iz + 1);
            if (/*!xc||*/ !yc || !zc) {
              ALWAYS_ASSERT_MSG(false,
                                "insanity in stripe_index_3d::sanity_check");
            }
          }
        }
        return true;
      }
    }
    return true;
  }
};

template <class Pt>
void copyPayload(Pt, void const *, int) {}

template <class Pt, class Payload>
void copyPayload(std::pair<Pt, Payload> &a, Payload const *b, int i) {
  a.second = b[i];
}

/// this is basically sorting the points along three different axes
/// so that we can loop over 9 stripes instead of 27 cells
template <class Pt, class Payload, class ShortIdx>
void stripe_index_3d<Pt, Payload, ShortIdx>::init(F width, Pt const *ptr,
                                                  size_t _Npts,
                                                  Payload const *ptr2) {
  Npts = _Npts;
  assert(Npts);

  width_ = width;
  width2_ = width * width;

  F xmn = 9e9, ymn = 9e9, zmn = 9e9;
  F xmx = -9e9, ymx = -9e9, zmx = -9e9;
  for (int i = 0; i < Npts; ++i) {
    xmn = std::min(xmn, ptr[i][0]);
    ymn = std::min(ymn, ptr[i][1]);
    zmn = std::min(zmn, ptr[i][2]);
    xmx = std::max(xmx, ptr[i][0]);
    ymx = std::max(ymx, ptr[i][1]);
    zmx = std::max(zmx, ptr[i][2]);
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
    int ix = (int)((ptr[i][0] - xmn /*+FUDGE*/) / width_);
    int iy = (int)((ptr[i][1] - ymn /*+FUDGE*/) / width_);
    int iz = (int)((ptr[i][2] - zmn /*+FUDGE*/) / width_);
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
  Value *tmpvalues = new Value[Npts + 4];  // space for 4 overflow ptr
  for (int i = 0; i < 4; ++i) {
    get_first_if_pair(tmpvalues[Npts + i])[0] = 9e9;
    get_first_if_pair(tmpvalues[Npts + i])[1] = 9e9;
    get_first_if_pair(tmpvalues[Npts + i])[2] = 9e9;
  }
  ushort *gridc = new ushort[gsize];
  for (int i = 0; i < gsize; ++i) gridc[i] = 0;
  for (int i = 0; i < Npts; ++i) {
    int const ix = (int)((ptr[i][0] - xmn /*+FUDGE*/) / width_);
    int const iy = (int)((ptr[i][1] - ymn /*+FUDGE*/) / width_);
    int const iz = (int)((ptr[i][2] - zmn /*+FUDGE*/) / width_);
    int const ig = ix + xdim_ * iy + xdim_ * ydim_ * iz;
    int const idx = tmp_pindex2[ig].first + gridc[ig];
    get_first_if_pair(tmpvalues[idx])[0] = ptr[i][0] - xmn /*+FUDGE*/;
    get_first_if_pair(tmpvalues[idx])[1] = ptr[i][1] - ymn /*+FUDGE*/;
    get_first_if_pair(tmpvalues[idx])[2] = ptr[i][2] - zmn /*+FUDGE*/;
    copyPayload(tmpvalues[idx], ptr2, i);
    ++(gridc[ig]);
  }
  values_ = tmpvalues;
  translation_[0] = /* FUDGE*/ -xmn;
  translation_[1] = /* FUDGE*/ -ymn;
  translation_[2] = /* FUDGE*/ -zmn;

  xmx_ = xmx - xmn /*+FUDGE*/ + width_;
  ymx_ = ymx - ymn /*+FUDGE*/ + width_;
  zmx_ = zmx - zmn /*+FUDGE*/ + width_;

  delete gridc;
  delete tmp_pindex2;

  {
    F xmn = 9e9, ymn = 9e9, zmn = 9e9;
    F xmx = -9e9, ymx = -9e9, zmx = -9e9;
    for (int i = 0; i < Npts; ++i) {
      xmn = std::min(xmn, get_first_if_pair(values_[i])[0]);
      ymn = std::min(ymn, get_first_if_pair(values_[i])[1]);
      zmn = std::min(zmn, get_first_if_pair(values_[i])[2]);
      xmx = std::max(xmx, get_first_if_pair(values_[i])[0]);
      ymx = std::max(ymx, get_first_if_pair(values_[i])[1]);
      zmx = std::max(zmx, get_first_if_pair(values_[i])[2]);
    }
    // std::cout << "delta " << translation_.transpose() << std::endl;
    // std::cout << "mnval " << xmn << " " << ymn << " " << zmn << std::endl;
    // std::cout << "mxval " << xmx << " " << ymx << " " << zmx << std::endl;
    ALWAYS_ASSERT(fabs(xmn) < 0.001);
    ALWAYS_ASSERT(fabs(ymn) < 0.001);
    ALWAYS_ASSERT(fabs(zmn) < 0.001);
  }
}
}  // end index ns
}  // end rif ns
