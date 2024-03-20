// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/VorticalFlow.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the single-material compressible flow
    equations
  \details   This file defines a Problem policy class for the single-material
    compressible flow equations, defined under PDE/CompFlow/. See
    PDE/CompFlow/Problem.h for general requirements on Problem policy classes
    for CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowProblemVorticalFlow_h
#define CompFlowProblemVorticalFlow_h

#include <string>
#include <unordered_set>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Problem.hpp"
#include "Inciter/InputDeck/New2InputDeck.hpp"
#include "EoS/GetMatProp.hpp"

namespace inciter {

extern ctr::New2InputDeck g_inputdeck;

//! CompFlow system of PDEs problem: vortical flow
//! \see Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
//!   equations with relevance to Inertial Confinement Fusion", Journal of
//!   Computational Physics 267 (2014) 196-209.
class CompFlowProblemVorticalFlow {

  private:
    using ncomp_t = tk::ncomp_t;
    using eq = newtag::compflow;

  public:
    //! Initialize numerical solution
    static tk::InitializeFn::result_type
    initialize( ncomp_t, const std::vector< EOS >&,
                tk::real x, tk::real y, tk::real z, tk::real );

    //! Evaluate analytical solution at (x,y,z) for all components
    static tk::InitializeFn::result_type
    analyticSolution( ncomp_t, const std::vector< EOS >&,
                      tk::real x, tk::real y, tk::real z, tk::real );

    //! Compute and return source term for vortical flow manufactured solution
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \param[in,out] sv Source term vector
    //! \note The function signature must follow tk::SrcFn
    static tk::SrcFn::result_type
    src( ncomp_t, const std::vector< EOS >& mat_blk,
         tk::real x, tk::real y, tk::real z, tk::real,
         std::vector< tk::real >& sv )
    {
      Assert(sv.size() == 5, "Incorrect source vector size");

      // manufactured solution parameters
      auto a = g_inputdeck.get< eq, newtag::alpha >();
      auto b = g_inputdeck.get< eq, newtag::beta >();
      // ratio of specific heats
      tk::real g = getmatprop< newtag::gamma >();
      // evaluate solution at x,y,z
      auto s = initialize( 5, mat_blk, x, y, z, 0.0 );

      // density source
      sv[0] = 0.0;
      // momentum source
      sv[1] = a*s[1]/s[0] - b*s[2]/s[0];
      sv[2] = b*s[1]/s[0] + a*s[2]/s[0];
      sv[3] = 0.0;
      // energy source
      sv[4] = (sv[1]*s[1] + sv[2]*s[2])/s[0] + 8.0*a*a*a*z*z/(g-1.0);
    }

    //! Return analytic field names to be output to file
    std::vector< std::string > analyticFieldNames( ncomp_t ) const;

    //! Return names of integral variables to be output to diagnostics file
    std::vector< std::string > names( ncomp_t ) const;

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::VORTICAL_FLOW; }
};

} // inciter::

#endif // CompFlowProblemVorticalFlow_h
