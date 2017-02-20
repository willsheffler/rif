#ifndef INCLUDED_scheme_nest_maps_ScaleMap_HH
#define INCLUDED_scheme_nest_maps_ScaleMap_HH

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/make_signed.hpp>
#include <iostream>
#include <vector>
#include "util/SimpleArray.hpp"
#include "util/template_loop.hpp"

namespace scheme {
namespace nest {
namespace pmap {

///@brief Parameter Mapping Policy for cartesian grids
///@tparam DIM the dimension number of the input parameter space
///@tparam Value the output value type, default
///@tparam Index index type, default size_t
///@tparam Float float type, default double
///@note NEST num_cells MUST agree with cell_sizes_
///@note bounds and cell indices are represented as SimpleArrays (like params)
///NOT Value Types
template <int DIM, class Value = util::SimpleArray<DIM, double>,
          class Index = uint64_t, class Float = typename Value::Scalar>
struct ScaleMap {
  static int const DIMENSION = DIM;
  typedef ScaleMap<DIM, Value, Index, Float> ThisType;
  typedef Value ValueType;
  typedef Float FloatType;
  typedef Index IndexType;
  typedef typename boost::make_signed<Index>::type SignedIndex;
  typedef util::SimpleArray<DIM, Index> Indices;
  typedef util::SimpleArray<DIM, SignedIndex> SignedIndices;
  typedef util::SimpleArray<DIM, Float> Params;
  BOOST_STATIC_ASSERT_MSG(DIM > 0, "ScaleMap DIM must be > 0");

 private:
  ///@brief lower bound on value space
  Params lower_bound_;
  ///@brief upper bound on value space base size 1
  Params upper_bound_;
  ///@brief distributes cell_index accross dimensions
  Params cell_width_;
  Indices cell_sizes_;
  Indices cell_sizes_pref_sum_;
  Index num_cells_;

 public:
  Params const& lower_bound() const { return lower_bound_; }
  Params const& upper_bound() const { return upper_bound_; }
  Params const& cell_width() const { return cell_width_; }
  Indices const& cell_sizes() const { return cell_sizes_; }

  static std::string pmap_name() {
    return "ScaleMap<" + boost::lexical_cast<std::string>(DIM) + ">";
  }

  ///@brief construct with default lb, ub, bs
  ScaleMap() {
    cell_sizes_.fill(1);
    lower_bound_.fill(0);
    upper_bound_.fill(1);
    init();
  }
  ///@brief construct with default lb, ub
  template <class I>
  ScaleMap(I const& bs) : cell_sizes_(bs) {
    lower_bound_.fill(0);
    upper_bound_.fill(1);
    init();
  }
  ///@brief construct with default bs
  template <class P>
  ScaleMap(P const& lb, P const& ub) : lower_bound_(lb), upper_bound_(ub) {
    cell_sizes_.fill(1);
    init();
  }
  ///@brief construct with specified lb, ub and bs
  template <class P, class I>
  ScaleMap(P const& lb, P const& ub, I const& bs)
      : lower_bound_(lb), upper_bound_(ub), cell_sizes_(bs) {
    init();
  }

  template <class I>
  void init(I const& bs) {
    cell_sizes_ = bs;
    lower_bound_.fill(0);
    upper_bound_.fill(1);
    init();
  }

  template <class P>
  void init(P const& lb, P const& ub) {
    lower_bound_ = lb;
    upper_bound_ = ub;
    cell_sizes_.fill(1);
    init();
  }

  template <class P, class I>
  void init(P const& lb, P const& ub, I const& bs) {
    lower_bound_ = lb;
    upper_bound_ = ub;
    cell_sizes_ = bs;
    init();
  }

  ///@brief sets up cell_size_pref_sum
  void init() {
    num_cells_ = cell_sizes_.prod();
    for (size_t i = 0; i < DIM; ++i) {
      cell_sizes_pref_sum_[i] = cell_sizes_.prod(i);
      cell_width_[i] =
          (upper_bound_[i] - lower_bound_[i]) / (Float)cell_sizes_[i];
      assert(upper_bound_[i] > lower_bound_[i]);
    }
  }

  ///@brief sets value based on cell_index and parameters using geometric bounds
  ///@return false iff invalid parameters
  bool params_to_value(Params const& params, Index cell_index, Index resl,
                       Value& value) const {
    for (size_t i = 0; i < DIM; ++i) {
      assert(cell_sizes_[i] > 0);
      assert(cell_sizes_[i] < 100000);
      assert(lower_bound_[i] < upper_bound_[i]);
      Float bi = (cell_index / cell_sizes_pref_sum_[i]) % cell_sizes_[i];
      value[i] = lower_bound_[i] + cell_width_[i] * (bi + params[i]);
    }
    return true;
  }

  ///@brief sets params/cell_index from value
  bool value_to_params(Value const& value, Index resl, Params& params,
                       Index& cell_index) const {
    value_to_params_for_cell(value, resl, params, 0);
    cell_index = 0;
    for (size_t i = 0; i < DIM; ++i) {
      assert(cell_sizes_[i] > 0);
      assert(cell_sizes_[i] < 100000);
      assert(lower_bound_[i] < upper_bound_[i]);
      // Index cell_size_pref_sum = cell_sizes_.head(i).prod();
      Float ci = (Index)params[i];
      cell_index += cell_sizes_pref_sum_[i] * ci;
      params[i] -= (Float)ci;
      assert(0.0 <= params[i] && params[i] <= 1.0);
    }
    return true;
  }

  ///@brief sets params/cell_index from value
  void value_to_params_for_cell(Value const& value, Index resl, Params& params,
                                Index cell_index) const {
    Indices cell_indices = cellindex_to_indices(cell_index);
    // std::cout << "CIDX " << cell_indices.transpose() << std::endl;
    for (size_t i = 0; i < DIM; ++i) {
      assert(cell_sizes_[i] > 0);
      assert(cell_sizes_[i] < 100000);
      assert(lower_bound_[i] < upper_bound_[i]);
      // Index cell_size_pref_sum = cell_sizes_.head(i).prod();
      params[i] = (value[i] - lower_bound_[i]) / cell_width_[i];
      params[i] -= (Float)cell_indices[i];
    }
  }

  Index indices_to_cellindex(Indices const& indices) const {
    Index index = 0;
    for (size_t i = 0; i < DIM; ++i) {
      assert(indices[i] < cell_sizes_[i]);
      index += indices[i] * cell_sizes_pref_sum_[i];
    }
    return index;
  }

  Indices cellindex_to_indices(Index index) const {
    Indices indices;
    for (size_t i = 0; i < DIM; ++i) {
      indices[i] = (index / cell_sizes_pref_sum_[i]) % cell_sizes_[i];
      assert(indices[i] < cell_sizes_[i]);
    }
    return indices;
  }

  template <class OutIter>
  void push_cell_index(SignedIndices const& indices, OutIter out) const {
    *(out++) = indices_to_cellindex(indices.template cast<size_t>());
  }

  // template<class OutIter>
  // void get_neighbors(Indices const & indices, Index cell_index, Index resl,
  // OutIter out)  {
  // 	// std::cout << indices.transpose() << std::endl;
  // 	SignedIndices lb = ((indices.template cast<int>()-1).max(     0     ));
  // 	SignedIndices ub = ((indices.template cast<int>()+1).min((1<<resl)-1));
  // 	// std::cout << "IX " << indices.transpose() << " cell " << cell_index
  // << std::endl;
  // 	// std::cout << "LB " << lb.transpose() << std::endl;
  // 	// std::cout << "UB " << ub.transpose() << std::endl;
  // 	boost::function<void(SignedIndices)> functor;
  // 	functor = boost::bind( & ThisType::template push_index<OutIter>, this,
  // _1, cell_index, resl, out );
  // 	util::NESTED_FOR<DIM>(lb,ub,functor);
  // }

  ///@brief return the cell_index of neighboring cells within delta of value
  ///@note delta parameter is in "Parameter Space"
  template <class OutIter>
  void get_neighboring_cells(Value const& value, Index resl, Float param_delta,
                             OutIter out) const {
    // Float param_delta = 1.0 / (Float)(1<<resl);
    assert(param_delta > 0);
    // convert to value space, decided against this
    // Params delta_param;
    // for(size_t i = 0; i < DIM; ++i) delta_param[i] = delta / cell_width_[i];
    Params params;
    value_to_params_for_cell(value, resl, params, 0);
    SignedIndex const BIG = 12345678;
    SignedIndices lb = (params - param_delta + (Float)BIG)
                           .max((Float)BIG)
                           .template cast<SignedIndex>() -
                       BIG;
    SignedIndices ub = (params + param_delta)
                           .template cast<Index>()
                           .min(cell_sizes_ - (Index)1)
                           .template cast<SignedIndex>();
    // std::cout << "PM " << params.transpose() << std::endl;
    // std::cout << "DL " << delta_param.transpose() << std::endl;
    // std::cout << "LB " << lb.transpose() << std::endl;
    // std::cout << "UB " << ub.transpose() << std::endl;
    boost::function<void(SignedIndices)> functor;
    functor = boost::bind(&ThisType::template push_cell_index<OutIter>, this,
                          _1, out);
    util::NESTED_FOR<DIM>(lb, ub, functor);
  }

  ///@brief aka covering radius max distance from bin center to any value within
  ///bin
  Float bin_circumradius(Index resl) const {
    Params width =
        (upper_bound_ - lower_bound_) / cell_sizes_.template cast<Float>();
    return 0.5 / (Float)(1 << resl) *
           sqrt((width * width).sum());  // squaredNorm
  }

  ///@brief maximum distance from the bin center which must be within the bin
  Float bin_inradius(Index resl) const {
    Params width =
        (upper_bound_ - lower_bound_) / cell_sizes_.template cast<Float>();
    return 1.5 / (Float)(1 << resl) * width.minCoeff();  // norm
  }

  ///@brief cell size
  Index num_cells() const { return num_cells_; }
  virtual ~ScaleMap() {}
};

template <int DIM, class Value, class Index, class Float>
std::ostream& operator<<(std::ostream& out,
                         ScaleMap<DIM, Value, Index, Float> const& sm) {
  out << "ScaleMap cell_sizes = " << sm.cell_sizes()
      << " cell_widths = " << sm.cell_width();
  return out;
}
}
}
}

#endif
