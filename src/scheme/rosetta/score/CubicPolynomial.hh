// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin CubicPolynomial.hh
///
/// @brief
/// Cubic polynomial functions for use in analytic energy functions.
///
/// @authors
/// Andrew Leaver-Fay (Analytic Functions)
/// Alex Ford (Refactoring)
///
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_scheme_rosetta_objective_CubicPolynomial_hh
#define INCLUDED_scheme_rosetta_objective_CubicPolynomial_hh

namespace scheme { namespace rosetta { namespace score {

/// @brief Simple struct holding coefficients of a cubic polynomial of the form
//  'c3 * x^3 + c2 * x^2 + c1 * x + c0'.
template< typename Real >
struct CubicPolynomialParamsBase {
	Real c0, c1, c2, c3;
	CubicPolynomialParamsBase() : c0(0), c1(0), c2(0), c3(0) {}
	CubicPolynomialParamsBase(Real a, Real b, Real c, Real d) : c0(a),c1(b),c2(c),c3(d) {}
};

template< typename Real >
Real
eval_cubic_polynomial(
	Real const x,
	Real const c0,
	Real const c1,
	Real const c2,
	Real const c3
)
{
	return ((c3*x+c2)*x+c1)*x + c0;
}

template< typename Real, typename CubicPolynomialParams >
Real
eval_cubic_polynomial(
	Real const x,
	CubicPolynomialParams const & cp
)
{
	return ((cp.c3*x+cp.c2)*x+cp.c1)*x + cp.c0;
}

template< typename Real, typename CubicPolynomialParams >
Real
cubic_polynomial_deriv(
	Real const x,
	CubicPolynomialParams const & cp
)
{
	return (3*cp.c3*x + 2*cp.c2)*x + cp.c1;
}

template< typename Real >
Real
cubic_polynomial_deriv(
	Real const x,
	Real const c0,
	Real const c1,
	Real const c2,
	Real const c3
)
{
	return (3*c3*x + 2*c2)*x + c1;
}

}}}



#endif
