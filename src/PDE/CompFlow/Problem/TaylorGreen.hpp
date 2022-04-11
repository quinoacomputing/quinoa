// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/TaylorGreen.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problem.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowProblemTaylorGreen_h
#define CompFlowProblemTaylorGreen_h

#include <string>
#include <unordered_set>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"
#include "EoS/EoS_Base.hpp"

namespace inciter {

//! CompFlow system of PDEs problem: Taylor-Green
//! \see G.I. Taylor, A.E. Green, "Mechanism of the Production of Small Eddies
//!   from Large Ones", Proc. R. Soc. Lond. A 1937 158 499-521; DOI:
//!   10.1098/rspa.1937.0036. Published 3 February 1937
//! \see Waltz, et. al, "Verification of a three-dimensional unstructured finite
//!   element method using analytic and manufactured solutions", Computers and
//!   Fluids, 2013, Vol.81, pp.57-67.
class CompFlowProblemTaylorGreen {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::compflow;
    static constexpr ncomp_t m_ncomp = 5;    //!< Number of scalar components

  public:
    //! Initialize numerical solution
    static tk::InitializeFn::result_type
    initialize( ncomp_t system, ncomp_t,
                tk::real x, tk::real y, tk::real, tk::real );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static tk::InitializeFn::result_type
    analyticSolution( ncomp_t system, ncomp_t,
                      std::vector< EoS_Base* >, tk::real x, tk::real y,
                      tk::real, tk::real );

    //! Compute and return source term for Rayleigh-Taylor manufactured solution
    //! \param[in] x X coordinate where to evaluate the source
    //! \param[in] y Y coordinate where to evaluate the source
    //! \param[in,out] sv Source term vector
    //! \note The function signature must follow tk::SrcFn
    static tk::SrcFn::result_type
    src( ncomp_t, ncomp_t, tk::real x, tk::real y, tk::real, tk::real,
         std::vector< tk::real >& sv )
    {
      Assert(sv.size() == 5, "Incorrect source vector size");
      sv[0] = sv[1] = sv[2] = sv[3] = 0.0;
      sv[4] = 3.0*M_PI/8.0*( cos(3.0*M_PI*x)*cos(M_PI*y)
                        - cos(3.0*M_PI*y)*cos(M_PI*x) );
    }

    //! Return field names to be output to file
    std::vector< std::string > analyticFieldNames( ncomp_t ) const;

    //! Return names of integral variables to be output to diagnostics file
    std::vector< std::string > names( ncomp_t ) const;

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::TAYLOR_GREEN; }
};

} // inciter::

#endif // CompFlowProblemTaylorGreen_h
