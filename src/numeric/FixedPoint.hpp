#ifndef INCLUDED_numeric_FixedPoint_HH
#define INCLUDED_numeric_FixedPoint_HH

#include <limits>

namespace rif {
namespace numeric {

template <int Divisor, class Data = unsigned char>
struct FixedPoint {
  typedef std::numeric_limits<Data> NL;
  Data data_;
  FixedPoint() { data_ = 0; }
  // FixedPoint(double d){ *this = d; }
  FixedPoint(float d) { *this = d; }
  void operator=(float d) {
    data_ = std::max((int)NL::min(),
                     std::min((int)NL::max(), (int)(d * (float)Divisor)));
  }
  // void operator =(float  d){ (Data)(data_ = d*(float )Divisor); }
  // operator float(){ return (double)data_/(double)Divisor; }
  operator float() const { return (float)data_ / (float)Divisor; }
  // operator int(){ return (int)data_; }
};
}
}

#endif
