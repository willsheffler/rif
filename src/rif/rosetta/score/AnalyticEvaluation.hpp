// -*-
// mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t
// -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// under license.
// (c) The Rosetta software is developed by the contributing members of the
// Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about
// this can be
// (c) addressed to University of Washington UW TechTransfer, email:
// license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin AnalyticEvaluation.hh
///
/// @brief
/// Analytic lj and lk-solvation functions used for analytic evaluation of
/// fa_atr/fa_rep and fa_sol.
///
/// @detailsed
/// These functions are called by the Etable class when analytic evaluation is
/// specified. The intention of these functions is to make available a modular
/// implemenation of
/// the lj and lk-solvation potentials for use in alternate scoring systems.
///
///
/// @authors
/// Andrew Leaver-Fay (Analytic Functions)
/// Alex Ford (Refactoring)
///
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_rosetta_objective_AnalyticEvaluation_hh
#define INCLUDED_rosetta_objective_AnalyticEvaluation_hh

#include <cmath>

#ifndef assert
#include <cassert>
#endif

#include "rosetta/score/CubicPolynomial.hpp"

namespace rif {
namespace rosetta {
namespace score {

// /// LK-Solvation Implementation

// template <typename Real>
// struct LKParamsBase
// {
//  Real lj_radius_1;
//  Real lj_radius_2;

//  Real fasol_cubic_poly_close_flat;

//  CubicPolynomialParamsBase<Real> fasol_cubic_poly_close;
//  Real fasol_cubic_poly_close_end;
//  Real fasol_cubic_poly_close_start;

//  Real lk_coeff1;
//  Real lk_coeff2;
//  Real lk_inv_lambda2_1;
//  Real lk_inv_lambda2_2;

//  CubicPolynomialParamsBase<Real> fasol_cubic_poly_far;
//  Real fasol_cubic_poly_far_end;
//  Real fasol_cubic_poly_far_start;

//  Real fasol_final_weight;
// };

/// @brief Evaluate the mututal desolvation energy as atom 1 and atom 2
/// approach.
/// Combine the desolvation of atom 1 by atom 2 into the same value as the
/// desolvation
/// of atom 2 by atom 1.
template <typename Real, typename LKParams>
void lk_evaluation(LKParams const &p, Real const dis, Real const inv_dis2,
                   Real &fa_solE) {
  /// a) At distances below p.fasol_cubic_poly_close_start, the value of fasol
  /// is held constant.
  /// b) Then there's a cubic_poly to smooth between this constant region and
  /// the exponential region.
  /// c) Then the exponential region.
  /// d) Then the cubic_poly to smooth between the exponential region and where
  /// the energy goes to zero.
  if (dis < p.fasol_cubic_poly_close_start) {
    fa_solE = p.fasol_cubic_poly_close_flat * p.fasol_final_weight;
  } else if (dis < p.fasol_cubic_poly_close_end) {
    fa_solE = eval_cubic_polynomial(dis, p.fasol_cubic_poly_close);
    fa_solE *= p.fasol_final_weight;
  } else if (dis < p.fasol_cubic_poly_far_start) {
    /// exponential evaluation
    Real const dis_rad1 = dis - p.lj_radius_1;
    Real const x1 = (dis_rad1 * dis_rad1) * p.lk_inv_lambda2_1;
    Real const dis_rad2 = dis - p.lj_radius_2;
    Real const x2 = (dis_rad2 * dis_rad2) * p.lk_inv_lambda2_2;

    fa_solE =
        inv_dis2 * (std::exp(-x1) * p.lk_coeff1 + std::exp(-x2) * p.lk_coeff2);
    fa_solE *= p.fasol_final_weight;

  } else if (dis < p.fasol_cubic_poly_far_end) {
    fa_solE = eval_cubic_polynomial(dis, p.fasol_cubic_poly_far);
    fa_solE *= p.fasol_final_weight;
  } else {
    fa_solE = 0;
  }
}

/// @brief Evaluate the lk solvation energy derivative of given atom pair.
template <typename Real, typename LKParams>
void lk_derivative(LKParams const &p, Real const dis, Real const inv_dis,
                   Real const inv_dis2, Real &dfasolE_ddis) {
  /// a) At distances below p.fasol_cubic_poly_close_start, the value of fasol
  /// is held constant.
  /// b) Then there's a cubic_poly to smooth between this constant region and
  /// the exponential region.
  /// c) Then the exponential region.
  /// d) Then the cubic_poly to smooth between the exponential region and where
  /// the energy goes to zero.
  if (dis < p.fasol_cubic_poly_close_start) {
    dfasolE_ddis = 0;
  } else if (dis < p.fasol_cubic_poly_close_end) {
    dfasolE_ddis = cubic_polynomial_deriv(dis, p.fasol_cubic_poly_close);
    dfasolE_ddis *= p.fasol_final_weight;
  } else if (dis < p.fasol_cubic_poly_far_start) {
    /// exponential evaluation
    /// assert( atype1 <= atype2 ), which is accomplished at the top of this
    /// function.
    Real const dis_rad1 = dis - p.lj_radius_1;
    Real const x1 = (dis_rad1 * dis_rad1) * p.lk_inv_lambda2_1;
    Real const dis_rad2 = dis - p.lj_radius_2;
    Real const x2 = (dis_rad2 * dis_rad2) * p.lk_inv_lambda2_2;

    Real const solvE1 = std::exp(-x1) * p.lk_coeff1 * inv_dis2;
    Real const solvE2 = std::exp(-x2) * p.lk_coeff2 * inv_dis2;
    dfasolE_ddis = -2 * ((dis_rad1 * p.lk_inv_lambda2_1 + inv_dis) * solvE1 +
                         (dis_rad2 * p.lk_inv_lambda2_2 + inv_dis) * solvE2);
    dfasolE_ddis *= p.fasol_final_weight;

  } else if (dis < p.fasol_cubic_poly_far_end) {
    dfasolE_ddis = cubic_polynomial_deriv(dis, p.fasol_cubic_poly_far) *
                   p.fasol_final_weight;
  } else {
    dfasolE_ddis = 0;
  }
}

// /// LJ Implementation
// template <typename Real>
// struct LJQuadraticRepulsion
// {
//  Real xlo;
//  Real xhi;
//  Real slope;
//  Real extrapolated_slope;
//  Real ylo;
// };

// template <typename Real>
// struct LJParamsBase
// {
//  Real ljrep_linear_ramp_d2_cutoff;
//  Real lj_switch_intercept;
//  Real lj_switch_slope;

//  Real lj_r12_coeff;
//  Real lj_r6_coeff;

//  Real ljatr_cubic_poly_xhi;
//  Real ljatr_cubic_poly_xlo;
//  CubicPolynomialParamsBase<Real> ljatr_cubic_poly_parameters;

//  Real lj_minimum;
//  Real lj_val_at_minimum;

//  LJQuadraticRepulsion<Real> ljrep_extra_repulsion;

//  bool ljrep_from_negcrossing;
//  Real ljatr_final_weight;
// };

/// @details only to be called when the distance, dis, is less than the
/// switch-point for
/// the lj_switch_dis2sigma value.
// template < typename Real, typename LJParams >
// Real
// ljrep_linearized(
//  Real const dis,
//  LJParams const & p
// )
// {
//  assert( dis * dis < p.ljrep_linear_ramp_d2_cutoff );
//  return  dis*p.lj_switch_slope + p.lj_switch_intercept;
// }

/// @details: evaluate the standard Lennard-Jones 6-12 functional form.
/// Only call this function if the square distance is in the range
/// sqrt( p.ljrep_linear_ramp_d2_cutoff ) < dis < p.ljatr_cubic_poly_xhi
template <typename Real, typename LJParams>
Real lj_generic_form(Real const dis2, Real const inv_dis2, LJParams const &p) {
  assert(dis2 >= p.ljrep_linear_ramp_d2_cutoff);
  assert(dis2 <= p.ljatr_cubic_poly_xhi * p.ljatr_cubic_poly_xhi);
  Real const inv_dis6 = inv_dis2 * inv_dis2 * inv_dis2;
  // Real const inv_dis12 = inv_dis6 * inv_dis6;

  return (p.lj_r12_coeff * inv_dis6 + p.lj_r6_coeff) * inv_dis6;
}

/// @details: evaluate the attractive component of the LJ term as it
/// ramps to zero.
/// Only call this function if the square distance is in the range
/// p.ljatr_cubic_poly_xlo < dis < p.ljatr_cubic_poly_xhi
template <typename Real, typename LJParams>
Real ljatr_cubic_poly_ramp_to_zero(Real const dis, LJParams const &p) {
  assert(dis >= p.ljatr_cubic_poly_xlo);
  assert(dis <= p.ljatr_cubic_poly_xhi);
  return eval_cubic_polynomial(dis, p.ljatr_cubic_poly_parameters);
}

template <typename Real, typename LJParams>
void lj_evaluation(LJParams const &p, Real const dis, Real const dis2,
                   Real const inv_dis2, Real &lj_atrE, Real &lj_repE) {
  // locals
  Real ljE;

  Real atrE = 0.;
  Real repE = 0.;

  lj_atrE = 0.;
  lj_repE = 0.;

  if (dis2 < p.ljrep_linear_ramp_d2_cutoff) {  // dis * p.inv_lj_sigma <
                                               // lj_switch_dis2sigma ) {
    //  ctsa - use linear ramp instead of lj when the dis/sigma  ratio drops
    //  below theshold
    // ljE = ljrep_linearized( dis, p ); // sometimes got assertion errors due
    // to numerical issues
    ljE = dis * p.lj_switch_slope + p.lj_switch_intercept;
  } else if (dis < p.ljatr_cubic_poly_xlo) {
    //  ctsa - calc regular lennard-jones
    ljE = lj_generic_form(
        dis2, inv_dis2,
        p);  //  p.lj_r12_coeff * inv_dis12 + p.lj_r6_coeff * inv_dis6;
  } else if (dis < p.ljatr_cubic_poly_xhi) {
    ljE = ljatr_cubic_poly_ramp_to_zero(dis, p);
  } else {
    return;
  }

  /// Divvy up the lennard-jones energies into attractive and repulsive
  /// components;
  /// for most atom pairs, the attractive component goes smoothly to a constant,
  /// as the atoms approach, and then the repulsive component takes over from
  /// there.
  if (p.ljrep_from_negcrossing) {
    // only for the REPLS and HREPS atom types: start repelling when the
    // lennard-jones term
    // goes from being attractive to repulsive.
    if (ljE < 0)
      atrE = ljE;
    else
      repE = ljE;
  } else if (dis < p.lj_minimum) {
    atrE = p.lj_val_at_minimum;
    repE = ljE - p.lj_val_at_minimum;
  } else {
    atrE = ljE;
  }

  /// Some atom pairs include extra repulsion that's modeled as a quadratic
  /// function
  /// in some region and then linearized outside of that region.  Specifically,
  /// this code
  /// exists for the OCbb / Ocbb interactions.
  if (dis < p.ljrep_extra_repulsion.xhi) {
    if (dis < p.ljrep_extra_repulsion.xlo) {
      repE += (dis - p.ljrep_extra_repulsion.xlo) *
                  p.ljrep_extra_repulsion.extrapolated_slope +
              p.ljrep_extra_repulsion.ylo;
    } else {
      repE += (p.ljrep_extra_repulsion.xhi - dis) *
              (p.ljrep_extra_repulsion.xhi - dis) *
              p.ljrep_extra_repulsion.slope;
    }
  }

  /// Now zero out the attractive energy if necessry; this is done for hydrogen
  /// interactions
  lj_atrE = atrE * p.ljatr_final_weight;
  lj_repE = repE;
}

/// @details: evaluate the derivative for the attractive component
/// of the LJ term as it ramps to zero.
/// Only call this function if the square distance is in the range
/// p.ljatr_cubic_poly_xlo < dis < p.ljatr_cubic_poly_xhi
template <typename Real, typename LJParams>
Real ljatr_cubic_poly_ramp_to_zero_deriv(Real const dis, LJParams const &p) {
  assert(dis >= p.ljatr_cubic_poly_xlo);
  assert(dis <= p.ljatr_cubic_poly_xhi);
  return cubic_polynomial_deriv(dis, p.ljatr_cubic_poly_parameters);
}

template <typename Real, typename LJParams>
void lj_derivatives(LJParams const &p, Real const dis, Real const inv_dis,
                    Real const dis2, Real const inv_dis2, Real &dljatrE_ddis,
                    Real &dljrepE_ddis) {
  // locals
  Real dljE(1), inv_dis6(1);

  if (dis2 < p.ljrep_linear_ramp_d2_cutoff) {  // dis * p.inv_lj_sigma <
                                               // lj_switch_dis2sigma ) {
    //  ctsa - use linear ramp instead of lj when the dis/sigma  ratio drops
    //  below theshold
    dljE = p.lj_switch_slope;
  } else if (dis < p.ljatr_cubic_poly_xlo) {
    //  ctsa - calc regular lennard-jones
    inv_dis6 = inv_dis2 * inv_dis2 * inv_dis2;
    Real const inv_dis7 = inv_dis * inv_dis6;

    dljE = inv_dis7 * (-12. * p.lj_r12_coeff * inv_dis6 - 6. * p.lj_r6_coeff);
  } else if (dis < p.ljatr_cubic_poly_xhi) {
    dljE = ljatr_cubic_poly_ramp_to_zero_deriv(dis, p);
  } else {
    /// assuming ljatr_cubic_poly_xhi == LK distance cutoff, since this will
    /// skip lk evaluation
    return;
  }

  /// Divvy up the lennard-jones energies into attractive and repulsive
  /// components;
  /// for most atom pairs, the attractive component goes smoothly to a constant,
  /// as the atoms approach, and then the repulsive component takes over from
  /// there.
  if (p.ljrep_from_negcrossing) {
    // only for the REPLS and HREPS atom types: start repelling when the
    // lennard-jones term
    // goes from being attractive to repulsive.
    // so, calculate the energy to decide whether it's attractive or repulsive
    if (dis2 < p.ljrep_linear_ramp_d2_cutoff ||
        inv_dis6 * (p.lj_r12_coeff * inv_dis6 + p.lj_r6_coeff) > 0.0) {
      dljrepE_ddis = dljE;
    }
  } else if (dis < p.lj_minimum) {
    dljatrE_ddis = 0;
    dljrepE_ddis = dljE;
  } else {
    dljatrE_ddis = dljE;
  }

  /// Some atom pairs include extra repulsion that's modeled as a quadratic
  /// function
  /// in some region and then linearized outside of that region.  Specifically,
  /// this code
  /// exists for the OCbb / Ocbb interactions.
  if (dis < p.ljrep_extra_repulsion.xhi) {
    if (dis < p.ljrep_extra_repulsion.xlo) {
      dljrepE_ddis += p.ljrep_extra_repulsion.extrapolated_slope;
    } else {
      dljrepE_ddis += -2 * (p.ljrep_extra_repulsion.xhi - dis) *
                      p.ljrep_extra_repulsion.slope;
    }
  }

  /// Finally, zero out the attractive energy derivatives if necessry; this is
  /// done for hydrogen interactions
  dljatrE_ddis *= p.ljatr_final_weight;
}
}
}
}

#endif
