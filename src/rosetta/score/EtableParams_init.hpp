#ifndef INCLUDED_rosetta_objective_EtableParams_init_hh
#define INCLUDED_rosetta_objective_EtableParams_init_hh

#include <vector>
#include "rosetta/score/EtableParams.hpp"

namespace scheme {
namespace rosetta {
namespace score {

struct EtableParamsInit {
  ///@brief horrible function to fill horrible rosetta datastructure of LJ/LK
  ///params
  static void init_EtableParams(
      std::vector<EtableParamsOnePair<float> >& analytic_parameters);
};
}
}
}

#endif
