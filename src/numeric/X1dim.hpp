#ifndef INCLUDED_numeric_X1dim_HH
#define INCLUDED_numeric_X1dim_HH

#include <cmath>
#include <iostream>

namespace scheme {
namespace numeric {

///@brief Simple linear "Xform" class for testing
struct X1dim {
  double val_;
  X1dim() : val_(0) {}
  X1dim(double d) : val_(d) {}
  bool operator==(X1dim const &o) const { return o.val_ == val_; }
  bool operator<(X1dim const &o) const { return val_ < o.val_; }
  template <class Archive>
  void serialize(Archive &ar, const unsigned int) {
    ar &val_;
  }
  static X1dim Identity() { return X1dim(0.0); }
  double &operator[](int i) {
    assert(i == 0);
    return val_;
  }
  double operator[](int i) const {
    assert(i == 0);
    return val_;
  }
};
inline X1dim operator*(X1dim a, X1dim b) { return X1dim(a.val_ + b.val_); }
inline std::ostream &operator<<(std::ostream &out, X1dim const &x) {
  return out << x.val_;
}
inline X1dim inverse(X1dim x) { return X1dim(-x.val_); }
inline double distance(X1dim a, X1dim b) { return std::abs(a.val_ - b.val_); }
}
}

#endif
