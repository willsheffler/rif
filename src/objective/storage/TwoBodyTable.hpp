#ifndef INCLUDED_objective_storage_TwoBodyTable_HH
#define INCLUDED_objective_storage_TwoBodyTable_HH

#include "util/SimpleArray.hpp"
#include "util/assert.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>

#include <set>

namespace scheme {
namespace objective {
namespace storage {

// MUST do things in this order:
// fill in onebody
// call init_onebody_filter
// fill in twobody
// NOTE: global/local rotamer number mapping is done here
// global/local residue numbering MUST be handled in the client code
template <class _Data = float>
struct TwoBodyTable {
  typedef _Data Data;
  typedef boost::multi_array<Data, 2> Array2D;
  typedef boost::multi_array<Array2D, 2> TwoBody;

  size_t nres_, nrot_;
  Array2D onebody_;
  boost::multi_array<int, 2> all2sel_, sel2all_;
  std::vector<int> nsel_;
  TwoBody twobody_;

  TwoBodyTable() {}  // for use with load()

  TwoBodyTable(size_t nres, size_t nrots) { init(nres, nrots); }

  void init(size_t nres, size_t nrots) {
    nres_ = nres;
    nrot_ = nrots;
    onebody_.resize(boost::extents[nres][nrots]);
    all2sel_.resize(boost::extents[nres][nrots]);
    sel2all_.resize(boost::extents[nres][nrots]);
    nsel_.resize(nres, 0);
    twobody_.resize(boost::extents[nres][nres]);
  }

  Data const &onebody(int ires, int irot) const { return onebody_[ires][irot]; }

  void bounds_check_1b(int ires, int irot) const {
    ALWAYS_ASSERT_MSG(ires >= 0, "ires < 0!");
    ALWAYS_ASSERT_MSG(irot >= 0, "irot < 0!");
    ALWAYS_ASSERT_MSG(ires < nres_, "ires >= nres!");
    ALWAYS_ASSERT_MSG(irot < nrot_, "irot >= nrot!");
  }

  Data const &onebody_at(int ires, int irot) const {
    bounds_check_1b(ires, irot);
    return onebody_[ires][irot];
  }

  void set_onebody(int ires, int irot, Data const &val) {
    bounds_check_1b(ires, irot);
    onebody_[ires][irot] = val;
  }

  Data twobody(int ires, int jres, int irot, int jrot) const {
    int const ir = ires > jres ? ires : jres;
    int const jr = ires > jres ? jres : ires;
    if (twobody_[ir][jr].num_elements() > 0) {
      int const irotlocal = all2sel_[ires][irot];
      int const jrotlocal = all2sel_[jres][jrot];
      if (irotlocal < 0 || jrotlocal < 0) {
        return Data(9e9);  // at least one of the onebodies is bad
      }
      // swap if jres > ires
      int const irl = ires > jres ? irotlocal : jrotlocal;
      int const jrl = ires > jres ? jrotlocal : irotlocal;
      return twobody_[ir][jr][irl][jrl];
    } else {
      return Data(0.0);
    }
  }

  Data twobody_rotlocalnumbering(int ires, int jres, int irotlocal,
                                 int jrotlocal) const {
    int const ir = ires > jres ? ires : jres;
    int const jr = ires > jres ? jres : ires;
    int const irl = ires > jres ? irotlocal : jrotlocal;
    int const jrl = ires > jres ? jrotlocal : irotlocal;
    if (twobody_[ir][jr].num_elements() > 0) {
      return twobody_[ir][jr][irl][jrl];
    } else {
      return Data(0.0);
    }
  }

  // assumes onebody energies have been filled in at this point!
  void init_onebody_filter(float thresh) {
    all2sel_.resize(boost::extents[nres_][nrot_]);
    sel2all_.resize(boost::extents[nres_][nrot_]);
    nsel_.resize(nres_, 0);
    for (int ires = 0; ires < nres_; ++ires) {
      int i = 0;
      for (int irot = 0; irot < nrot_; ++irot) {
        if (onebody_[ires][irot] <= thresh) {
          all2sel_[ires][irot] = i;
          sel2all_[ires][i] = irot;
          ++i;
        } else {
          all2sel_[ires][irot] = -1;
        }
        nsel_[ires] = i;
        for (int j = i; j < nrot_; ++j) {
          sel2all_[ires][j] = -1;  // past end of selected
        }
      }
    }
  }
  void init_twobody(int ires, int jres) {
    twobody_[ires][jres].resize(boost::extents[nsel_[ires]][nsel_[jres]]);
  }
  void clear_twobody(int ires, int jres) {
    twobody_[ires][jres].resize(boost::extents[0][0]);
  }
  int twobody_mem_use() const {
    int memuse = 0;
    for (int i = 0; i < nres_; ++i) {
      for (int j = 0; j < nres_; ++j) {
        memuse += twobody_[i][j].num_elements() * sizeof(_Data);
        // std::cout << i << " " << j << " " << twobody_[i][j].num_elements() <<
        // std::endl;
      }
    }
    return memuse;
  }

  bool check_equal(TwoBodyTable<Data> const &that) const {
    bool iseq = true;
    iseq &= nres_ == that.nres_;
    iseq &= nrot_ == that.nrot_;
    if (!iseq) return false;
    for (int i = 0; i < nres_ * nrot_; ++i) {
      iseq &= sel2all_.data()[i] == that.sel2all_.data()[i];
      iseq &= all2sel_.data()[i] == that.all2sel_.data()[i];
      iseq &= onebody_.data()[i] == that.onebody_.data()[i];
    }
    for (int i = 0; i < nres_; ++i) {
      iseq &= nsel_[i] == that.nsel_[i];
    }
    for (int ir = 0; ir < nres_; ++ir) {
      for (int jr = 0; jr < nres_; ++jr) {
        iseq &= twobody_[ir][jr].num_elements() ==
                that.twobody_[ir][jr].num_elements();
        if (!iseq) return false;
        for (int k = 0; k < twobody_[ir][jr].num_elements(); ++k) {
          iseq &= twobody_[ir][jr].data()[k] == that.twobody_[ir][jr].data()[k];
        }
      }
    }
    return iseq;
  }

  void save(std::ostream &out, std::string const &description) const {
    ALWAYS_ASSERT(onebody_.num_elements() == nres_ * nrot_);
    ALWAYS_ASSERT(all2sel_.num_elements() == nres_ * nrot_);
    ALWAYS_ASSERT(sel2all_.num_elements() == nres_ * nrot_);
    ALWAYS_ASSERT(nsel_.size() == nres_);
    size_t const dsrcrize = description.size();
    out.write((char *)&dsrcrize, sizeof(size_t));
    out.write(description.c_str(), description.size() * sizeof(char));
    out.write((char *)&nres_, sizeof(size_t));
    out.write((char *)&nrot_, sizeof(size_t));
    for (int i = 0; i < nres_ * nrot_; ++i) {
      out.write((char *)&(onebody_.data()[i]), sizeof(Data));
      out.write((char *)&(all2sel_.data()[i]), sizeof(int));
      out.write((char *)&(sel2all_.data()[i]), sizeof(int));
    }
    for (int i = 0; i < nres_; ++i) {
      out.write((char *)&(nsel_[i]), sizeof(int));
    }
    for (int ir = 0; ir < nres_; ++ir) {
      for (int jr = 0; jr < nres_; ++jr) {
        size_t const N = twobody_[ir][jr].num_elements();
        if (N != 0 && N != nsel_[ir] * nsel_[jr]) {
          std::cout << "bad N: " << N << " should be 0 or "
                    << nsel_[ir] * nsel_[jr] << std::endl;
          ALWAYS_ASSERT(N == 0 || N == nsel_[ir] * nsel_[jr]);
        }
        out.write((char *)&N, sizeof(size_t));
        for (int k = 0; k < N; ++k) {
          out.write((char *)&(twobody_[ir][jr].data()[k]), sizeof(Data));
        }
      }
    }
  }
  void load(std::istream &in, std::string &description) {
    // std::cout << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << " " <<
    // "" << std::endl;
    size_t dsrcrize;
    in.read((char *)&dsrcrize, sizeof(size_t));
    // std::cout << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << " " <<
    // dsrcrize << std::endl;
    char *buf = new char[dsrcrize];
    in.read(buf, dsrcrize * sizeof(char));
    description.resize(dsrcrize);
    for (int i = 0; i < dsrcrize; ++i) description[i] = buf[i];
    delete[] buf;
    in.read((char *)&nres_, sizeof(size_t));
    in.read((char *)&nrot_, sizeof(size_t));
    onebody_.resize(boost::extents[nres_][nrot_]);
    all2sel_.resize(boost::extents[nres_][nrot_]);
    sel2all_.resize(boost::extents[nres_][nrot_]);
    for (int i = 0; i < nres_ * nrot_; ++i) {
      in.read((char *)&(onebody_.data()[i]), sizeof(Data));
      in.read((char *)&(all2sel_.data()[i]), sizeof(int));
      in.read((char *)&(sel2all_.data()[i]), sizeof(int));
    }
    nsel_.resize(nres_);
    for (int i = 0; i < nres_; ++i) {
      in.read((char *)&(nsel_[i]), sizeof(int));
    }
    twobody_.resize(boost::extents[nres_][nres_]);
    for (int ir = 0; ir < nres_; ++ir) {
      for (int jr = 0; jr < nres_; ++jr) {
        size_t const N = twobody_[ir][jr].num_elements();
        ALWAYS_ASSERT(N == 0 || N == nsel_[ir] * nsel_[jr]);
        in.read((char *)&N, sizeof(size_t));
        if (N == 0) {
          twobody_[ir][jr].resize(boost::extents[0][0]);
        } else {
          twobody_[ir][jr].resize(boost::extents[nsel_[ir]][nsel_[jr]]);
          for (int k = 0; k < N; ++k) {
            in.read((char *)&(twobody_[ir][jr].data()[k]), sizeof(Data));
            // static bool doit = true;
            // if( doit ){
            // 	twobody_[ir][jr].data()[k] += 1.0;
            // 	doit = false;
            // }
          }
        }
      }
    }
  }

  shared_ptr<TwoBodyTable<Data>> create_subtable(
      std::vector<bool> const &res_selection,
      std::vector<std::vector<float>> const
          &new1b,  // always in global numbering
      float filter1bthresh) const {
    shared_ptr<TwoBodyTable<Data>> newt_p = make_shared<TwoBodyTable<Data>>();
    TwoBodyTable &newt(*newt_p);
    // size_t nres_, nrot_;
    // Array2D onebody_;
    // boost::multi_array< int, 2 > all2sel_, sel2all_;
    // std::vector<int> nsel_;
    // TwoBody twobody_;
    // ALWAYS_ASSERT_MSG( &newt != this, "don't call create_subtable with
    // self!!!" );
    ALWAYS_ASSERT(res_selection.size() == nres_);
    newt.nrot_ = nrot_;  // assume same rot_index / global rotamer numbering
    newt.nres_ = 0;
    std::vector<int> res_g2l, res_l2g;
    for (int i = 0; i < nres_; ++i) {
      if (res_selection[i] == true) {
        newt.nres_++;
        ALWAYS_ASSERT(new1b[i].size() == nrot_);
        res_l2g.push_back(i);
        res_g2l.push_back(res_g2l.size());
      } else {
        res_g2l.push_back(-1);
      }
    }
    newt.onebody_.resize(boost::extents[newt.nres_][newt.nrot_]);
    for (int ilocal = 0; ilocal < newt.nres_; ++ilocal) {
      int iglobal = res_l2g[ilocal];
      for (int irot = 0; irot < newt.nrot_; ++irot) {
        newt.onebody_[ilocal][irot] = new1b[iglobal][irot];
      }
    }
    // std::cout << "create_subtable: new nres: " << newt.nres_ << std::endl;
    newt.init_onebody_filter(
        filter1bthresh);  // inits & fills all2sel_, sel2all_, and nsel_
    newt.twobody_.resize(boost::extents[newt.nres_][newt.nres_]);
    for (int ilocal = 0; ilocal < newt.nres_; ++ilocal) {
      for (int jlocal = 0; jlocal < newt.nres_; ++jlocal) {
        int iglobal = res_l2g[ilocal];
        int jglobal = res_l2g[jlocal];
        newt.init_twobody(ilocal, jlocal);
        Data minscore = 9e9, maxscore = -9e9;
        if (twobody_[iglobal][jglobal].num_elements() == 0) {
          // no table in old table, subtable will also have nothing
          newt.clear_twobody(ilocal, jlocal);
        } else {
          for (int ilocalrot = 0; ilocalrot < newt.nsel_[ilocal]; ++ilocalrot) {
            int iglobalrot = newt.sel2all_[ilocal][ilocalrot];
            int ioldrot = all2sel_[iglobal][iglobalrot];
            for (int jlocalrot = 0; jlocalrot < newt.nsel_[jlocal];
                 ++jlocalrot) {
              int jglobalrot = newt.sel2all_[jlocal][jlocalrot];
              int joldrot = all2sel_[jglobal][jglobalrot];
              Data score = 9e9;
              if (ioldrot >= 0 && joldrot >= 0) {
                score = this->twobody_[iglobal][jglobal][ioldrot][joldrot];
              }
              newt.twobody_[ilocal][jlocal][ilocalrot][jlocalrot] = score;
              minscore = std::min(minscore, score);
              maxscore = std::max(maxscore, score);
            }
          }
          if (minscore > -0.01 && maxscore < 0.01) {
            ALWAYS_ASSERT(0 <= ilocal && ilocal < newt.nres_);
            ALWAYS_ASSERT(0 <= jlocal && jlocal < newt.nres_);
            newt.clear_twobody(ilocal, jlocal);
          }
        }
      }
    }
    return newt_p;
  }

  // void deepcopy( TwoBodyTable<Data> const & that ) {
  // 	// Array2D onebody_;
  // 	// boost::multi_array< int, 2 > all2sel_, sel2all_;
  // 	// std::vector<int> nsel_;
  // 	// TwoBody twobody_;
  // 	nres_ = that.nres_;
  // 	nrot_ = that.nrot_;
  // 	onebody_.resize(
  // boost::extents[that.onebody_.shape()[0]][that.onebody_.shape()[1]] );//
  // this seems crappy...
  // 	all2sel_.resize(
  // boost::extents[that.all2sel_.shape()[0]][that.all2sel_.shape()[1]] );
  // 	sel2all_.resize(
  // boost::extents[that.sel2all_.shape()[0]][that.sel2all_.shape()[1]] );
  // 	ALWAYS_ASSERT( all2sel_.num_elements() == onebody_.num_elements() );
  // 	ALWAYS_ASSERT( sel2all_.num_elements() == onebody_.num_elements() );
  // 	for( int i = 0; i < that.onebody_.num_elements(); ++i ){
  // 		onebody_.data()[i] = that.onebody_.data()[i];
  // 		all2sel_.data()[i] = that.all2sel_.data()[i];
  // 		sel2all_.data()[i] = that.sel2all_.data()[i];
  // 	}
  // 	nsel_.resize( that.nsel_.size() );
  // 	for( int i = 0; i < nsel_.size(); ++i ) nsel_[i] = that.nsel_[i];
  // 	twobody_.resize(
  // boost::extents[that.twobody_.shape()[0]][that.twobody_.shape()[1]] );
  // 	for( int i = 0; i < nres_; ++i ){
  // 	for( int j = 0; j < nres_; ++j ){
  // 		twobody_[i][j].resize(
  // boost::extents[that.twobody_[i][j].shape()[0]][that.twobody_[i][j].shape()[1]]
  // );
  // 		ALWAYS_ASSERT( 0 == twobody_[i][j].num_elements() ||
  // nsel_[i]*nsel_[j] == twobody_[i][j].num_elements() );
  // 		for( int k = 0; k < twobody_[i][j].num_elements(); ++k ){
  // 			twobody_[i][j].data()[k] =
  // that.twobody_[i][j].data()[k];
  // 		}
  // 	}}

  // }
};
}
}
}

#endif
