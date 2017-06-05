// gross but efficient code from c/gpu stuff many years ago
// todo: refactor wtih visitor to remove extra code

#pragma once

#include "rif/eigen_types.hpp"
#include "rif/util/assert.hpp"

namespace rif {
namespace index {

template <class Pt, class Idx = uint16_t>
struct OneSide3dIndex {
  using F = typename Pt::Scalar;
  typedef struct { Idx x, y; } ushort2;
  typedef V3<F> Vec;

  F width_, width2_;  // cell width
  size_t Npts;
  Pt const *pts_ = nullptr;
  ushort2 const *pindex_ = nullptr;
  int xdim_, ydim_, zdim_;
  float xmx_, ymx_, zmx_;
  Vec translation_;

  OneSide3dIndex() {}

  template <class C>
  OneSide3dIndex(C const &pts, F width) {
    Pt const *ptr = &pts[0];
    init(width, ptr, pts.size());
  }

  OneSide3dIndex(F width, Pt const *ptr, size_t Npts) {
    init(width, ptr, Npts);
  }

  void init(F width, Pt const *ptr, size_t Npts);

  virtual ~OneSide3dIndex() {
    if (pts_) delete pts_;
    if (pindex_) delete pindex_;
  }

  bool sanity_check() const {
    for (int ix = 0; ix < xdim_; ++ix) {
      for (int iy = 0; iy < ydim_; ++iy) {
        for (int iz = 0; iz < zdim_; ++iz) {
          // std::cout << ix << " " << iy << " " << iz << endl;
          ushort const ig = ix + xdim_ * iy + ydim_ * xdim_ * iz;
          ushort const igl = pindex_[ig].x;
          ushort const igu = pindex_[ig].y;
          for (int i = igl; i < igu; ++i) {
            // float const & x(pts_[i].x);
            float const &y(pts_[i].y());
            float const &z(pts_[i].z());
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

  int brute_nbcount(Vec const &v_in) {
    Vec const v = v_in + translation_;
    int count = 0;
    for (size_t i = 0; i < Npts; ++i) {
      Pt const &a2 = pts_[i];
      float const d2 = (v[0] - a2[0]) * (v[0] - a2[0]) +
                       (v[1] - a2[1]) * (v[1] - a2[1]) +
                       (v[2] - a2[2]) * (v[2] - a2[2]);
      if (d2 <= width2_) {
        ++count;
      }
    }
    return count;
  }

  // todo: express with visitor instead iff pref same
  int nbcount(Vec const &v_in) const {
    Vec const v = v_in + translation_;
    float x = v.x();
    float y = v.y();
    float z = v.z();
    if (x < -width_ || y < -width_ || z < -width_) return 0;  // worth it iff
    if (x > xmx_ || y > ymx_ || z > zmx_) return 0;           // worth it iff
    int count = 0;
    int const ix =
        (x < 0) ? 0 : std::min(xdim_ - 1, static_cast<int>(x / width_));
    int const iy0 = (y < 0) ? 0 : static_cast<int>(y / width_);
    int const iz0 = (z < 0) ? 0 : static_cast<int>(z / width_);
    int const iyl = std::max(0, iy0 - 1);
    int const izl = std::max(0, iz0 - 1);
    int const iyu = std::min(static_cast<int>(ydim_), iy0 + 2);
    int const izu =
        std::min(static_cast<int>(zdim_), static_cast<int>(iz0) + 2);
    for (int iy = iyl; iy < iyu; ++iy) {
      for (int iz = izl; iz < izu; ++iz) {
        int const ig = ix + xdim_ * iy + xdim_ * ydim_ * iz;
        assert(ig < xdim_ * ydim_ * zdim_);
        assert(ix < xdim_);
        assert(iy < ydim_);
        assert(iz < zdim_);
        int const &igl = pindex_[ig].x;
        int const &igu = pindex_[ig].y;
        for (int i = igl; i < igu; ++i) {
          Pt const a2 = pts_[i];
          float const d2 = (x - a2.x()) * (x - a2.x()) +
                           (y - a2.y()) * (y - a2.y()) +
                           (z - a2.z()) * (z - a2.z());
          if (d2 <= width2_) {
            ++count;
          }
        }
      }
    }
    return count;
  }

  bool contact(Vec const &v_in) const {
    Vec const v = v_in + translation_;
    float x = v.x();
    float y = v.y();
    float z = v.z();
    if (x < -width_ || y < -width_ || z < -width_)
      return false;                                      // worth it iff
    if (x > xmx_ || y > ymx_ || z > zmx_) return false;  // worth it iff
    int const ix =
        (x < 0) ? 0 : std::min(xdim_ - 1, static_cast<int>(x / width_));
    int const iy0 = (y < 0) ? 0 : static_cast<int>(y / width_);
    int const iz0 = (z < 0) ? 0 : static_cast<int>(z / width_);
    int const iyl = std::max(0, iy0 - 1);
    int const izl = std::max(0, iz0 - 1);
    int const iyu = std::min(static_cast<int>(ydim_), iy0 + 2);
    int const izu =
        std::min(static_cast<int>(zdim_), static_cast<int>(iz0) + 2);
    for (int iy = iyl; iy < iyu; ++iy) {
      for (int iz = izl; iz < izu; ++iz) {
        int const ig = ix + xdim_ * iy + xdim_ * ydim_ * iz;
        assert(ig < xdim_ * ydim_ * zdim_);
        assert(ix < xdim_);
        assert(iy < ydim_);
        assert(iz < zdim_);
        int const &igl = pindex_[ig].x;
        int const &igu = pindex_[ig].y;
        for (int i = igl; i < igu; ++i) {
          Pt const a2 = pts_[i];
          float const d2 = (x - a2.x()) * (x - a2.x()) +
                           (y - a2.y()) * (y - a2.y()) +
                           (z - a2.z()) * (z - a2.z());
          if (d2 <= width2_) {
            return true;
          }
        }
      }
    }
    return false;
  }

  template <typename Visitor, typename P>
  void visit(P const &v_in, Visitor &visitor) const {
    P v(v_in);
    v.x() += translation_.x();
    v.y() += translation_.y();
    v.z() += translation_.z();
    float x = v.x();
    float y = v.y();
    float z = v.z();
    if (x < -width_ || y < -width_ || z < -width_) return;  // worth it iff
    if (x > xmx_ || y > ymx_ || z > zmx_) return;           // worth it iff
    int const ix = (x < 0) ? 0 : std::min(xdim_ - 1, (int)(x / width_));
    int const iy0 = (y < 0) ? 0 : y / width_;
    int const iz0 = (z < 0) ? 0 : z / width_;
    int const iyl = std::max(0, iy0 - 1);
    int const izl = std::max(0, iz0 - 1);
    int const iyu = std::min((int)ydim_, iy0 + 2);
    int const izu = std::min((int)zdim_, (int)iz0 + 2);
    for (int iy = iyl; iy < iyu; ++iy) {
      for (int iz = izl; iz < izu; ++iz) {
        int const ig = ix + xdim_ * iy + xdim_ * ydim_ * iz;
        assert(ig < xdim_ * ydim_ * zdim_);
        assert(ix < xdim_);
        assert(iy < ydim_);
        assert(iz < zdim_);
        int const &igl = pindex_[ig].x;
        int const &igu = pindex_[ig].y;
        for (int i = igl; i < igu; ++i) {
          Pt const &c = *((P *)pts_ + i);  // DANGER! assume type P has
                                           // same in-memory prefix as
                                           // Pt (works for
                                           // xyzVector<float>)
          float const d2 = (x - c.x()) * (x - c.x()) +
                           (y - c.y()) * (y - c.y()) +
                           (z - c.z()) * (z - c.z());
          if (d2 <= width2_) {
            visitor.visit(v, c, d2);
          }
        }
      }
    }
  }

  template <typename Visitor, typename P>
  void visit_lax(P const &v_in, Visitor &visitor) const {
    P v(v_in);
    v.x() += translation_.x();
    v.y() += translation_.y();
    v.z() += translation_.z();
    float x = v.x();
    float y = v.y();
    float z = v.z();
    if (x < -width_ || y < -width_ || z < -width_) return;  // worth it iff
    if (x > xmx_ || y > ymx_ || z > zmx_) return;           // worth it iff
    int const ix =
        (x < 0) ? 0 : std::min(xdim_ - 1, static_cast<int>(x / width_));
    int const iy0 = (y < 0) ? 0 : static_cast<int>(y / width_);
    int const iz0 = (z < 0) ? 0 : static_cast<int>(z / width_);
    int const iyl = std::max(0, iy0 - 1);
    int const izl = std::max(0, iz0 - 1);
    int const iyu = std::min(static_cast<int>(ydim_), iy0 + 2);
    int const izu =
        std::min(static_cast<int>(zdim_), static_cast<int>(iz0 + 2));
    for (int iy = iyl; iy < iyu; ++iy) {
      for (int iz = izl; iz < izu; ++iz) {
        int const ig = ix + xdim_ * iy + xdim_ * ydim_ * iz;
        assert(ig < xdim_ * ydim_ * zdim_);
        assert(ix < xdim_);
        assert(iy < ydim_);
        assert(iz < zdim_);
        int const &igl = pindex_[ig].x;
        int const &igu = pindex_[ig].y;
        for (int i = igl; i < igu; ++i) {
          P const &c = *((P *)(pts_ + i));
          visitor.visit(v, c);
        }
      }
    }
  }
};

template <class Pt, class Idx>
void OneSide3dIndex<Pt, Idx>::init(F width, Pt const *ptr, size_t _Npts) {
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

  xdim_ = static_cast<int>((xmx - xmn + 0.0001) / width_ + 0.999999);
  ydim_ = static_cast<int>((ymx - ymn + 0.0001) / width_ + 0.999999);
  zdim_ = static_cast<int>((zmx - zmn + 0.0001) / width_ + 0.999999);
  assert(xdim_ < 9999);
  assert(ydim_ < 9999);
  assert(zdim_ < 9999);
  int const gsize = xdim_ * ydim_ * zdim_;
  ushort2 *tmp_pindex = new ushort2[gsize];
  ushort2 *tmp_pindex2 = new ushort2[gsize];

  for (int i = 0; i < gsize; ++i) {
    tmp_pindex2[i].y = 0;
    tmp_pindex2[i].x = 0;
  }
  // TR<<"atom "<<Npts<<" grid1 "<<xdim_*ydim_*zdim_<<" "<<xdim_<<"
  // "<<ydim_<<" "<<zdim_<<std::endl;

  for (int i = 0; i < Npts; ++i) {
    int ix = static_cast<int>((ptr[i][0] - xmn /*+FUDGE*/) / width_);
    int iy = static_cast<int>((ptr[i][1] - ymn /*+FUDGE*/) / width_);
    int iz = static_cast<int>((ptr[i][2] - zmn /*+FUDGE*/) / width_);
    assert(ix >= 0);
    assert(iy >= 0);
    assert(iz >= 0);
    assert(ix < xdim_);
    assert(iy < ydim_);
    assert(iz < zdim_);
    int ig = ix + xdim_ * iy + xdim_ * ydim_ * iz;
    assert(ig >= 0);
    assert(ig < 9999999);
    ++(tmp_pindex2[ig].y);
  }
  for (int i = 1; i < gsize; ++i)
    tmp_pindex2[i].x = tmp_pindex2[i - 1].x + tmp_pindex2[i - 1].y;
  for (int i = 1; i < gsize; ++i)
    tmp_pindex2[i].y = tmp_pindex2[i].x + tmp_pindex2[i].y;
  for (int iz = 0; iz < zdim_; ++iz)
    for (int iy = 0; iy < ydim_; ++iy)
      for (int ix = 0; ix < xdim_; ++ix) {
        int const ixl = (int)std::max(0, (int)ix - 1);
        int const ixu = std::min(xdim_ - 1u, ix + 1u);
        int const ig0 = xdim_ * iy + xdim_ * ydim_ * iz;
        tmp_pindex[ix + ig0].x = tmp_pindex2[ixl + ig0].x;
        tmp_pindex[ix + ig0].y = tmp_pindex2[ixu + ig0].y;
      }
  // for(int iz = 0; iz < zdim_; ++iz) for(int iy = 0; iy < ydim_; ++iy)
  // for(int ix = 0; ix < xdim_; ++ix) {
  //       int i = ix+xdim_*iy+xdim_*ydim_*iz;
  //       TR<<ix<<" "<<iy<<" "<<iz<<" "<<I(3,tmp_pindex2[i].x)<<"
  //       "<<I(3,tmp_pindex2[i].y) <<" "<<I(3,tmp_pindex[i].x)<<"
  //       "<<I(3,tmp_pindex[i].y)<<std::endl;
  //     }
  pindex_ = tmp_pindex;
  Pt *gatom = new Pt[Npts + 4];  // space for 4 overflow ptr
  for (int i = 0; i < 4; ++i) {
    gatom[Npts + i][0] = 9e9;
    gatom[Npts + i][1] = 9e9;
    gatom[Npts + i][2] = 9e9;
  }
  ushort *gridc = new ushort[gsize];
  for (int i = 0; i < gsize; ++i) gridc[i] = 0;
  for (int i = 0; i < Npts; ++i) {
    int const ix = static_cast<int>((ptr[i][0] - xmn /*+FUDGE*/) / width_);
    int const iy = static_cast<int>((ptr[i][1] - ymn /*+FUDGE*/) / width_);
    int const iz = static_cast<int>((ptr[i][2] - zmn /*+FUDGE*/) / width_);
    int const ig = ix + xdim_ * iy + xdim_ * ydim_ * iz;
    int const idx = tmp_pindex2[ig].x + gridc[ig];
    gatom[idx][0] = ptr[i][0] - xmn /*+FUDGE*/;
    gatom[idx][1] = ptr[i][1] - ymn /*+FUDGE*/;
    gatom[idx][2] = ptr[i][2] - zmn /*+FUDGE*/;
    ++(gridc[ig]);
  }
  pts_ = gatom;
  translation_[0] = /* FUDGE*/ -xmn;
  translation_[1] = /* FUDGE*/ -ymn;
  translation_[2] = /* FUDGE*/ -zmn;
  xmx_ = xmx - xmn /*+FUDGE*/ + width_;
  ymx_ = ymx - ymn /*+FUDGE*/ + width_;
  zmx_ = zmx - zmn /*+FUDGE*/ + width_;
  // for(int iz = 0; iz < zdim(); ++iz) for(int iy = 0; iy < ydim(); ++iy)
  // for(int ix = 0; ix < xdim(); ++ix) {
  //       int i = ix+xdim_*iy+xdim_*ydim_*iz;
  //       TR<<"GRID CELL "<<ix<<" "<<iy<<" "<<iz<<std::endl;
  //       for(int ig = tmp_pindex2[i].x; ig < tmp_pindex2[i].y; ++ig) {
  //       TR<<F(7,3,gatom[ig].x)<<" "<<F(7,3,gatom[ig].y)<<"
  //       "<<F(7,3,gatom[ig].z)<<std::endl;
  //     }
  //   }
  delete gridc;
  delete tmp_pindex2;
}
}  // end index ns
}  // end rif ns
