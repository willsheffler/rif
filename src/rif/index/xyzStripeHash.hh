
template <class T, class Idx = uint16_t>
struct xyzStripeHash {
  typedef struct { Idx x, y; } IdxPair;
  typedef V3f Vec;

  // iterators:
  template <class C>
  struct iter_base : public std::iterator<std::input_iterator_tag, float> {
    iter_base(Ball const *p) : p_(p) {}
    C &operator=(C const &r) {
      p_ = r.p_;
      return *this;
    }
    C &operator++() {
      ++p_;
      return static_cast<C &>(*this);
    }
    bool operator!=(C const &r) const { return (p_ != r.p_); }
    bool operator==(C const &r) const { return (p_ == r.p_); }

   protected:
    Ball const *p_;
  };
  struct const_iterator : public iter_base<const_iterator> {
    const_iterator(Ball const *p) : iter_base<const_iterator>(p) {}
    Vec const &operator*() { return *((Vec const *)(this->p_)); }
    Vec const *operator->() { return ((Vec const *)(this->p_)); }
    float radius() { return this->p_->lj_radius(); }
  };

 public:
  xyzStripeHash(float grid_size = 0.0,
                utility::vector1<Ball> const &balls = utility::vector1<Ball>());

  void init(utility::vector1<Ball> const &balls);

  virtual ~xyzStripeHash() {
    if (grid_balls_) delete grid_balls_;
    if (grid_stripe_) delete grid_stripe_;
  }

  const_iterator begin() const { return const_iterator(grid_balls_); }
  const_iterator end() const { return const_iterator(grid_balls_ + nballs_); }

  bool sanity_check() const;

  std::string debug_pdb(Xform const &x = Xform::identity()) const;

  int nbcount(Vec const &v_in) const;
  int nbcount_raw(Vec const &v) const;
  bool clash(Vec const &v_in) const;
  bool clash_not_resid(Vec const &v_in, int const &resid,
                       int const &resid2 = 0) const;
  bool clash_raw(Vec const &v) const;
  float clash_amount(Vec const &v_in) const;

  // @brief Check if dist(b_in, hb) < (b_in.radius + hb.radius) for any hb in
  // hash.
  //
  // Input ball provided in global coordinate frame, clashes are evaluated in
  // global frame.
  //
  // Returns 1-based ball index if clash is found, 0 otherwise.
  int clash_check_ball(Ball const &b) const;

  // @brief Generate residue mapping (r_t, r) where:
  // 	any(dist(b_t, b) < (b_t.radius + b.radius)) for
  // 		{b_t in target_balls | b_t.resi = r_t}, {b in this | b.resi ==
  // r}
  //
  // Populated residue_pairs with the first identified clash per r_t, meaning
  // that
  // all clashing resi in test_balls are included in mapping, however not all
  // clashing
  // resi in this are guarenteed to be included.
  //
  // If residue_pairs already contains an entry for resi r_t the residue will
  // not
  // be checked, and the mapping will not be updated.
  //
  // @returns true if any clashing pair was identified, false otherwise.
  bool clash_check_residue_pairs(utility::vector1<Ball> const &test_balls,
                                 std::map<Size, Size> &residue_pairs) const;

  void fill_pairs(xyzVector_float const &v, int const &ir,
                  utility::vector1<std::pair<int, int> > &pairs,
                  float maxd2 = 0.0) const;

  template <typename Visitor, typename P>
  void visit(P const &v_in, Visitor &visitor) const {
    P v(v_in);
    v.x() += translation_.x();
    v.y() += translation_.y();
    v.z() += translation_.z();
    float x = v.x();
    float y = v.y();
    float z = v.z();
    if (x < -grid_size_ || y < -grid_size_ || z < -grid_size_)
      return;                                      // worth it iff
    if (x > xmx_ || y > ymx_ || z > zmx_) return;  // worth it iff
    int const ix = (x < 0) ? 0 : numeric::min(xdim_ - 1, (int)(x / grid_size_));
    int const iy0 = (y < 0) ? 0 : y / grid_size_;
    int const iz0 = (z < 0) ? 0 : z / grid_size_;
    int const iyl = numeric::max(0, iy0 - 1);
    int const izl = numeric::max(0, iz0 - 1);
    int const iyu = numeric::min((int)ydim_, iy0 + 2);
    int const izu = numeric::min((int)zdim_, (int)iz0 + 2);
    for (int iy = iyl; iy < iyu; ++iy) {
      for (int iz = izl; iz < izu; ++iz) {
        int const ig = ix + xdim_ * iy + xdim_ * ydim_ * iz;
        assert(ig < xdim_ * ydim_ * zdim_);
        assert(ix < xdim_);
        assert(iy < ydim_);
        assert(iz < zdim_);
        int const &igl = grid_stripe_[ig].x;
        int const &igu = grid_stripe_[ig].y;
        for (int i = igl; i < igu; ++i) {
          Ball const &c = *((P *)grid_balls_ + i);  // DANGER! assume type P has
                                                    // same in-memory prefix as
                                                    // Ball (works for
                                                    // xyzVector<float>)
          float const d2 = (x - c.x()) * (x - c.x()) +
                           (y - c.y()) * (y - c.y()) +
                           (z - c.z()) * (z - c.z());
          if (d2 <= grid_size2_) {
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
    if (x < -grid_size_ || y < -grid_size_ || z < -grid_size_)
      return;                                      // worth it iff
    if (x > xmx_ || y > ymx_ || z > zmx_) return;  // worth it iff
    int const ix =
        (x < 0) ? 0 : numeric::min(xdim_ - 1, static_cast<int>(x / grid_size_));
    int const iy0 = (y < 0) ? 0 : static_cast<int>(y / grid_size_);
    int const iz0 = (z < 0) ? 0 : static_cast<int>(z / grid_size_);
    int const iyl = numeric::max(0, iy0 - 1);
    int const izl = numeric::max(0, iz0 - 1);
    int const iyu = numeric::min(static_cast<int>(ydim_), iy0 + 2);
    int const izu =
        numeric::min(static_cast<int>(zdim_), static_cast<int>(iz0 + 2));
    for (int iy = iyl; iy < iyu; ++iy) {
      for (int iz = izl; iz < izu; ++iz) {
        int const ig = ix + xdim_ * iy + xdim_ * ydim_ * iz;
        assert(ig < xdim_ * ydim_ * zdim_);
        assert(ix < xdim_);
        assert(iy < ydim_);
        assert(iz < zdim_);
        int const &igl = grid_stripe_[ig].x;
        int const &igu = grid_stripe_[ig].y;
        for (int i = igl; i < igu; ++i) {
          P const &c = *((P *)(grid_balls_ + i));
          visitor.visit(v, c);
        }
      }
    }
  }

  Ball const *grid_atoms() const { return grid_balls_; }
  Size size() const { return nballs_; }
  int natom() const { return nballs_; }
  int xdim() const { return xdim_; }
  int ydim() const { return ydim_; }
  int zdim() const { return zdim_; }
  float grid_size() const { return grid_size_; }
  float grid_size2() const { return grid_size2_; }
  xyzVector_float const &translation() const { return translation_; }
  xyzVector<Real> translation_real() const {
    return xyzVector<Real>(translation_.x(), translation_.y(),
                           translation_.z());
  }

  IdxPair const *grid_stripe() const { return grid_stripe_; }

  Ball const &ball(Size const &ib) const {
    assert(1ul <= ib && ib <= (Size)nballs_);
    return grid_balls_[ib - 1];
  }
  xyzVector_float xyz(Size const &ib) const {
    assert(1ul <= ib && ib <= (Size)nballs_);
    return grid_balls_[ib - 1].xyz() - translation_;
  }
  Size resi(Size const &ib) const {
    assert(1ul <= ib && ib <= (Size)nballs_);
    return grid_balls_[ib - 1].resi();
  }

 private:
  float grid_size_, grid_size2_;
  int nballs_;
  Ball const *grid_balls_;
  IdxPair const *grid_stripe_;
  int xdim_, ydim_, zdim_;
  float xmx_, ymx_, zmx_;
  // numeric::xyzMatrix<Real> rotation_;
  numeric::xyzVector<float> translation_;
  // neighbor_iterator neighbor_end_;
};

}  // namespace hashing
}  // namespace geometry
}  // namespace numeric

#endif
