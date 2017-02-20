#ifndef INCLUDED_scheme_nest_maps_DiscreteChoiceMap_HH
#define INCLUDED_scheme_nest_maps_DiscreteChoiceMap_HH

#include <boost/bind.hpp>
#include <boost/static_assert.hpp>
#include <iostream>
#include <vector>
#include "util/template_loop.hpp"

namespace scheme {
namespace nest {
namespace pmap {

/// @brief Parameter Mapping Policy class which represents a discrete set of
/// choices for 0 dimensional Nests
/// @note NEST num_cells MUST agree with choices.size()
template <int DIM, class Value, class Index, class Float>
struct DiscreteChoiceMap {
  BOOST_STATIC_ASSERT_MSG(DIM == 0, "DiscreteChoiceMap DIM must be == 0");
  static int const DIMENSION = DIM;
  typedef char Params;
  std::vector<Value> choices;
  DiscreteChoiceMap(std::vector<Value> const &_choices)
      : choices(_choices), num_cells_(choices.size()) {}
  ///@brief sets value based only on cell_index
  ///@note params has no meaning for zero-dimensional nests, only cell_index
  ///@return false iff invalid parameters
  bool params_to_value(Params const & /*params*/, Index cell_index, Index resl,
                       Value &value) const {
    if (cell_index >= choices.size()) return false;
    value = choices[cell_index];
    return true;
  }
  Index num_cells() const { return num_cells_; }
  virtual ~DiscreteChoiceMap() {}

 private:
  Index num_cells_;
};
}
}
}

#endif
