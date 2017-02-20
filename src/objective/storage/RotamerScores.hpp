#ifndef INCLUDED_objective_storage_RotamerScores_HH
#define INCLUDED_objective_storage_RotamerScores_HH

#include <boost/lexical_cast.hpp>
#include "util/SimpleArray.hpp"
#include "util/assert.hpp"

#include <vector>

namespace scheme {
namespace objective {
namespace storage {

struct Empty {};

template <class _Data = uint16_t, int _RotamerBits = 9, int _Divisor = -13>
struct RotamerScore {
  BOOST_STATIC_ASSERT(_Divisor < 0);
  typedef RotamerScore<_Data, _RotamerBits, _Divisor> THIS;
  typedef _Data Data;
  static const int RotamerBits = _RotamerBits;
  static const int Divisor = _Divisor;
  static const int ScoreBits = sizeof(_Data) * 8 - _RotamerBits;
  static const Data one = 1;
  static const Data RotamerMask = ((one << RotamerBits) - one);
  static const bool UseSat = false;

  Data data_;
  RotamerScore() : data_(RotamerMask) {}
  RotamerScore(Data data) : data_(data) {}
  RotamerScore(Data rot, float score) {
    Data sdat = score * _Divisor;
    assert(sdat < (one << ScoreBits));
    assert(rot < (one << RotamerBits));
    data_ = rot | (sdat << RotamerBits);
  }
  float score() const { return data2float(get_score_data()); }
  Data rotamer() const { return data_ & RotamerMask; }
  void set_score(float score) { set_score_data(float2data(score)); }
  void set_rotamer(Data rot) {
    assert(rot < one << RotamerBits);
    data_ = (data_ & (~RotamerMask)) | rot;
  }

  void set_score_data(Data sd) {
    assert(sd < (one << ScoreBits));
    data_ = rotamer() | (sd << RotamerBits);
  }
  Data get_score_data() const { return data_ >> RotamerBits; }
  static float divisor() { return _Divisor; }

  static float data2float(Data data) { return float(data) / _Divisor; }
  static Data float2data(float f) { return Data(f * _Divisor); }

  bool operator<(THIS const& that) const {
    return data_ > that.data_;
  }  // reverse so low score is low
  bool operator==(THIS const& that) const { return data_ == that.data_; }
  bool operator!=(THIS const& that) const { return data_ != that.data_; }
  bool operator==(Data const& that) const { return data_ == that; }

  bool empty() const { return data_ == RotamerMask; }

  void set_or_merge(THIS const& that) {
    if (that < *this) {
      *this = that;
    }
  }

  static std::string name() {
    static std::string const name =
        "RotamerScore< " + std::string("data_size ") +
        boost::lexical_cast<std::string>(sizeof(Data)) + ", " +
        boost::lexical_cast<std::string>(THIS::RotamerBits) + ", " +
        boost::lexical_cast<std::string>(THIS::Divisor) + " >";
    return name;
  }

} __attribute__((packed));
template <class Data, int RBits, int Div>
int const RotamerScore<Data, RBits, Div>::RotamerBits;
template <class Data, int RBits, int Div>
int const RotamerScore<Data, RBits, Div>::Divisor;
template <class Data, int RBits, int Div>
int const RotamerScore<Data, RBits, Div>::ScoreBits;
template <class Data, int RBits, int Div>
bool const RotamerScore<Data, RBits, Div>::UseSat;
template <class Data, int RBits, int Div>
std::ostream& operator<<(std::ostream& out,
                         RotamerScore<Data, RBits, Div> const& val) {
  out << val.rotamer() << "<" << val.score() << ">";
  return out;
}

template <int RotBits = 0>
struct SatisfactionDatum {
  uint8_t data_;
  SatisfactionDatum() : data_(255) {}
  SatisfactionDatum(uint8_t d) : data_(d) {}
  bool empty() const { return data_ == 255; }
  bool not_empty() const { return data_ != 255; }
  int target_sat_num() const { return (int)data_; }
  int rotamer_sat_num() const { return 0; }
  bool operator==(SatisfactionDatum const& o) const { return data_ == o.data_; }
} __attribute__((packed));
template <int RotBits>
std::ostream& operator<<(std::ostream& out,
                         SatisfactionDatum<RotBits> const& val) {
  out << (int)val.data_;
  return out;
}

template <class _Data = uint16_t, int _RotamerBits = 9, int _Divisor = -13,
          class _SatDatum = SatisfactionDatum<>, int _NSat = 2>
struct RotamerScoreSat : public RotamerScore<_Data, _RotamerBits, _Divisor> {
  typedef RotamerScoreSat<_Data, _RotamerBits, _Divisor, _SatDatum, _NSat> THIS;
  typedef RotamerScore<_Data, _RotamerBits, _Divisor> BASE;
  typedef _Data Data;
  typedef _SatDatum SatDatum;
  static const bool UseSat = true;
  static const int RotamerBits = _RotamerBits;
  static const int Divisor = _Divisor;
  static const int NSat = _NSat;
  util::SimpleArray<NSat, SatDatum> sat_data_;
  RotamerScoreSat() : BASE() {}
  RotamerScoreSat(Data data) : BASE(data) {}
  RotamerScoreSat(Data rot, float score, int sat1 = -1, int sat2 = -1)
      : BASE(rot, score) {
    if (sat1 < 0 || NSat < 1) return;
    ALWAYS_ASSERT(sat1 < 256);
    sat_data_[0].data_ = (uint8_t)sat1;
    if (sat2 < 0 || NSat < 2) return;
    ALWAYS_ASSERT(sat2 < 256);
    sat_data_[1].data_ = (uint8_t)sat2;
    // std::cout << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << " " <<
    // sat1 << "/" << sat_data_[0] << " " << sat2 << "/" << sat_data_[1] <<
    // std::endl;
  }
  static std::string name() {
    static std::string const name =
        std::string("RotamerScoreSat< dsize=") +
        boost::lexical_cast<std::string>(sizeof(Data)) + ", " + "rbit=" +
        boost::lexical_cast<std::string>(BASE::RotamerBits) + ", " + "div=" +
        boost::lexical_cast<std::string>(BASE::Divisor) + ", " +
        std::string("nsat=") +
        boost::lexical_cast<std::string>(sizeof(SatDatum) * NSat) + " >";
    return name;
  }
  template <class Array>
  void get_sat_groups_raw(Array& sat_groups_out) const {
    for (int isat = 0; isat < NSat; ++isat) {
      sat_groups_out[isat] = sat_data_[isat].target_sat_num();
    }
  }
  void get_sat_groups(std::vector<int>& sat_groups_out) const {
    for (int isat = 0; isat < NSat; ++isat) {
      if (sat_data_[isat].not_empty()) {
        sat_groups_out.push_back(sat_data_[isat].target_sat_num());
      }
    }
  }
  void mark_sat_groups(std::vector<bool>& sat_groups_mask) const {
    for (int isat = 0; isat < NSat; ++isat) {
      if (sat_data_[isat].not_empty()) {
        sat_groups_mask[sat_data_[isat].target_sat_num()] = true;
      }
    }
  }
  bool is_new_sat(SatDatum const& sd) const {
    if (sd.empty()) return false;
    for (int i = 0; i < NSat; ++i) {
      if (sat_data_[i] == sd) return false;
    }
    return true;
  }
  void set_or_merge(THIS const& that) {
    if (this->empty()) {
      *this = that;
    } else {
      // merge sat data iff same rotamer, mostly for bounding grids
      if (that.rotamer() == this->rotamer()) {
        int osat = 0;
        for (int isat = 0; isat < NSat; ++isat) {
          if (sat_data_[isat].empty()) {
            while (osat < NSat && !is_new_sat(that.sat_data_[osat])) ++osat;
            if (osat >= NSat) break;
            sat_data_[isat] = that.sat_data_[osat];
            ++osat;
          }
        }
      }
      BASE::set_or_merge(that);
    }
  }

  bool operator==(THIS const& o) const {
    return this->data_ == o.data_ && this->sat_data_ == o.sat_data_;
  }
  bool operator!=(THIS const& o) const {
    return this->data_ != o.data_ || this->sat_data_ != o.sat_data_;
  }
  bool operator<(THIS const& that) const {
    int nsat = 0, onsat = 0;
    for (int i = 0; i < NSat; ++i) {
      nsat += this->sat_data_[i].not_empty();
      onsat += that.sat_data_[i].not_empty();
    }
    if (nsat == onsat) return this->data_ > that.data_;
    return nsat > onsat;
  }  // reverse so low score is low

} __attribute__((packed));
template <class Data, int RBits, int Div, class Sat, int N>
int const RotamerScoreSat<Data, RBits, Div, Sat, N>::RotamerBits;
template <class Data, int RBits, int Div, class Sat, int N>
int const RotamerScoreSat<Data, RBits, Div, Sat, N>::Divisor;
template <class Data, int RBits, int Div, class Sat, int N>
bool const RotamerScoreSat<Data, RBits, Div, Sat, N>::UseSat;
template <class Data, int RBits, int Div, class Sat, int N>
int const RotamerScoreSat<Data, RBits, Div, Sat, N>::NSat;
template <class Data, int RBits, int Div, class Sat, int N>
std::ostream& operator<<(std::ostream& out,
                         RotamerScoreSat<Data, RBits, Div, Sat, N> const& val) {
  out << val.rotamer() << "<" << val.score() << ">";
  std::vector<int> sat;
  val.get_sat_groups(sat);
  for (int j = 0; j < sat.size(); ++j) {
    out << "," << sat[j];
  }
  out << "  ";

  return out;
}

template <int _N, class _RotamerScore = RotamerScore<> >
struct RotamerScores {
  BOOST_STATIC_ASSERT((_N > 0));
  BOOST_STATIC_ASSERT((_N < 256));  // arbitrary

  typedef _RotamerScore RotScore;
  typedef typename RotScore::Data Data;
  typedef RotamerScores<_N, RotScore> THIS;

  static int const N = _N;
  util::SimpleArray<N, RotScore> rotscores_;

  RotamerScores() { rotscores_.fill(RotScore::RotamerMask); }

  // void add_rotamer( Data rot, float score ){
  // add_rotamer( RotScore(rot,score) );
  // }
  void add_rotamer(Data rot, float score, int sat1 = -1, int sat2 = -1) {
    add_rotamer_impl<RotScore::UseSat>(rot, score, sat1, sat2);
  }
  void rotamer_sat_groups(int irot, std::vector<int>& sat_groups_out) const {
    rotamer_sat_groups_impl<RotScore::UseSat>(irot, sat_groups_out);
  }
  void mark_sat_groups(int irot, std::vector<bool>& sat_groups_mask) const {
    mark_sat_groups_impl<RotScore::UseSat>(irot, sat_groups_mask);
  }
  template <class Array>
  void get_sat_groups_raw(int irot, Array& a) const {
    get_sat_groups_raw_impl<RotScore::UseSat, Array>(irot, a);
  }
  void add_rotamer(RotScore to_insert) {
    Data irot = to_insert.rotamer();
    int insert_pos = 0;
    RotScore worst(std::numeric_limits<Data>::max());
    for (int i = 0; i < N; ++i) {
      // std::cout << " iter " << i << " " << rotscores_[i].score() << " " <<
      // rotscores_[i].rotamer() << " "
      // << "cur " << rotscores_[i].data_ << " low " << worst.data_ <<
      // std::endl;
      // if rot already stored, this is the position we check
      if (rotscores_[i].rotamer() == irot) {
        insert_pos = i;
        // std::cout << "rotamer equal at " << i << std::endl;
        break;
      }
      // else we take the worst position
      if (worst < rotscores_[i]) {
        // std::cout << "worst is " << i << std::endl;
        worst = rotscores_[i];
        insert_pos = i;
      }
    }
    // std::cout << "insert_pos " << insert_pos << std::endl;
    // now insert if new val is better than worst stored val
    rotscores_[insert_pos].set_or_merge(to_insert);
  }
  template <int N2>
  void merge(RotamerScores<N2, RotScore> const& that) {
    for (int i = 0; i < N2; ++i) {
      if (that.empty(i)) break;
      add_rotamer(that.rotscores_[i]);
    }
  }
  float score_of_rotamer(int irot) const {
    for (int i = 0; i < N; ++i) {
      if (rotscores_[i].rotamer() == irot) {
        return rotscores_[i].score();
      }
    }
    return 0.0f;
  }

  float score(int i) const {
    assert(i < N);
    return rotscores_[i].score();
  }
  Data rotamer(int i) const {
    assert(i < N);
    return rotscores_[i].rotamer();
  }

  bool empty(int i) const { return rotscores_[i].empty(); }

  static int maxsize() { return _N; }

  int size() const {
    int i;
    for (i = 0; i < _N; ++i)
      if (rotscores_[i].empty()) break;
    return i;
  }

  void sort_rotamers() { std::sort(rotscores_.begin(), rotscores_.end()); }

  bool is_sorted() const {
    for (int i = 1; i < _N; ++i)
      if (rotscores_[i] < rotscores_[i - 1]) return false;
    return true;
  }

  static std::string name() {
    static std::string const name = std::string("RotamerScores< N=") +
                                    boost::lexical_cast<std::string>(_N) +
                                    ", " + RotScore::name() + " >";
    return name;
  }

  bool operator==(THIS const& o) const { return rotscores_ == o.rotscores_; }
  bool operator!=(THIS const& o) const { return rotscores_ != o.rotscores_; }

  template <bool UseSat>
  typename boost::enable_if_c<UseSat, void>::type add_rotamer_impl(Data rot,
                                                                   float score,
                                                                   int sat1,
                                                                   int sat2) {
    add_rotamer(RotScore(rot, score, sat1, sat2));
  }
  template <bool UseSat>
  typename boost::disable_if_c<UseSat, void>::type add_rotamer_impl(Data rot,
                                                                    float score,
                                                                    int sat1,
                                                                    int sat2) {
    add_rotamer(RotScore(rot, score));
  }
  template <bool UseSat>
  typename boost::enable_if_c<UseSat, void>::type rotamer_sat_groups_impl(
      int irot, std::vector<int>& sat_groups_out) const {
    rotscores_[irot].get_sat_groups(sat_groups_out);
  }
  template <bool UseSat>
  typename boost::disable_if_c<UseSat, void>::type rotamer_sat_groups_impl(
      int irot, std::vector<int>& sat_groups_out) const {
    return;
  }
  template <bool UseSat>
  typename boost::enable_if_c<UseSat, void>::type mark_sat_groups_impl(
      int irot, std::vector<bool>& sat_groups_mask) const {
    rotscores_[irot].mark_sat_groups(sat_groups_mask);
  }
  template <bool UseSat>
  typename boost::disable_if_c<UseSat, void>::type mark_sat_groups_impl(
      int irot, std::vector<bool>& sat_groups_mask) const {
    return;
  }
  template <bool UseSat, class Array>
  typename boost::enable_if_c<UseSat, void>::type get_sat_groups_raw_impl(
      int irot, Array& a) const {
    rotscores_[irot].get_sat_groups_raw(a);
  }
  template <bool UseSat, class Array>
  typename boost::disable_if_c<UseSat, void>::type get_sat_groups_raw_impl(
      int irot, Array&) const {
    return;
  }
};

template <int N, class R>
std::ostream& operator<<(std::ostream& out, RotamerScores<N, R> const& val) {
  out << val.name() << "( ";
  for (int i = 0; i < val.size(); ++i) {
    out << val.rotscores_[i] << " ";
  }
  out << ")";
  return out;
}
}
}
}

#endif
