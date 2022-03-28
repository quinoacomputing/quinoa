// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/RayleighTaylor.hpp
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
#ifndef CompFlowProblemRayleighTaylor_h
#define CompFlowProblemRayleighTaylor_h

#include <string>
#include <unordered_set>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "EoS/EoS.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! CompFlow system of PDEs problem: Rayleigh-Taylor
//! \see Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
//!   equations with relevance to Inertial Confinement Fusion", Journal of
//!   Computational Physics 267 (2014) 196-209.
class CompFlowProblemRayleighTaylor {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::compflow;

  public:
    //! Initialize numerical solution
    static tk::InitializeFn::result_type
    initialize( ncomp_t system, ncomp_t, tk::real x, tk::real y,
                tk::real z, tk::real t );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static tk::InitializeFn::result_type
    analyticSolution( ncomp_t system, ncomp_t, tk::real x, tk::real y,
                      tk::real z, tk::real t );

    //! Compute and return source term for Rayleigh-Taylor manufactured solution
    //! \param[in] system Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \param[in] t Physical time at which to evaluate the source
    //! \param[in,out] sv Source term vector
    //! \note The function signature must follow tk::SrcFn
    static tk::SrcFn::result_type
    src( ncomp_t system, ncomp_t, tk::real x, tk::real y, tk::real z, tk::real t,
         std::vector< tk::real >& sv )
    {
      Assert(sv.size() == 5, "Incorrect source vector size");
      using tag::param; using std::sin; using std::cos;

      // manufactured solution parameters
      auto a = g_inputdeck.get< param, eq, tag::alpha >()[system];
      auto bx = g_inputdeck.get< param, eq, tag::betax >()[system];
      auto by = g_inputdeck.get< param, eq, tag::betay >()[system];
      auto bz = g_inputdeck.get< param, eq, tag::betaz >()[system];
      auto k = g_inputdeck.get< param, eq, tag::kappa >()[system];
      auto p0 = g_inputdeck.get< param, eq, tag::p0 >()[system];
      // ratio of specific heats
      auto g = gamma< tag::compflow >(system);

      // evaluate solution at x,y,z,t
      auto s = initialize( system, 5, x, y, z, t );

      // density, velocity, energy, pressure
      auto rho = s[0];
      auto u = s[1]/s[0];
      auto v = s[2]/s[0];
      auto w = s[3]/s[0];
      auto E = s[4]/s[0];
      auto p = p0 + a*(bx*x*x + by*y*y + bz*z*z);

      // spatial gradients
      std::array< tk::real, 3 > drdx{{ -2.0*bx*x, -2.0*by*y, -2.0*bz*z }};
      std::array< tk::real, 3 > dpdx{{ 2.0*a*bx*x, 2.0*a*by*y, 2.0*a*bz*z }};
      tk::real ft = cos(k*M_PI*t);
      std::array< tk::real, 3 > dudx{{ ft*M_PI*z*cos(M_PI*x),
                                       0.0,
                                       ft*sin(M_PI*x) }};
      std::array< tk::real, 3 > dvdx{{ 0.0,
                                       -ft*M_PI*z*sin(M_PI*y),
                                       ft*cos(M_PI*y) }};
      std::array< tk::real, 3 > dwdx{{ ft*M_PI*0.5*M_PI*z*z*sin(M_PI*x),
                                       ft*M_PI*0.5*M_PI*z*z*cos(M_PI*y),
                                      -ft*M_PI*z*(cos(M_PI*x) - sin(M_PI*y)) }};
      std::array< tk::real, 3 > dedx{{
        dpdx[0]/rho/(g-1.0) - p/(g-1.0)/rho/rho*drdx[0]
        + u*dudx[0] + v*dvdx[0] + w*dwdx[0],
        dpdx[1]/rho/(g-1.0) - p/(g-1.0)/rho/rho*drdx[1]
        + u*dudx[1] + v*dvdx[1] + w*dwdx[1],
        dpdx[2]/rho/(g-1.0) - p/(g-1.0)/rho/rho*drdx[2]
        + u*dudx[2] + v*dvdx[2] + w*dwdx[2] }};

      // time derivatives
      auto dudt = -k*M_PI*sin(k*M_PI*t)*z*sin(M_PI*x);
      auto dvdt = -k*M_PI*sin(k*M_PI*t)*z*cos(M_PI*y);
      auto dwdt =  k*M_PI*sin(k*M_PI*t)/2*M_PI*z*z*(cos(M_PI*x) - sin(M_PI*y));
      auto dedt = u*dudt + v*dvdt + w*dwdt;

      // density source
      sv[0] = u*drdx[0] + v*drdx[1] + w*drdx[2];
      // momentum source
      sv[1] = rho*dudt+u*sv[0]+dpdx[0] + s[1]*dudx[0]+s[2]*dudx[1]+s[3]*dudx[2];
      sv[2] = rho*dvdt+v*sv[0]+dpdx[1] + s[1]*dvdx[0]+s[2]*dvdx[1]+s[3]*dvdx[2];
      sv[3] = rho*dwdt+w*sv[0]+dpdx[2] + s[1]*dwdx[0]+s[2]*dwdx[1]+s[3]*dwdx[2];
      // energy source
      sv[4] = rho*dedt + E*sv[0] + s[1]*dedx[0]+s[2]*dedx[1]+s[3]*dedx[2]
           + u*dpdx[0]+v*dpdx[1]+w*dpdx[2];
    }

    //! Return field names to be output to file
    std::vector< std::string > analyticFieldNames( ncomp_t ) const;

    //! Return names of integral variables to be output to diagnostics file
    std::vector< std::string > names( ncomp_t /*ncomp*/ ) const;

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::RAYLEIGH_TAYLOR; }
};

} // inciter::

#endif // CompFlowProblemRayleighTaylor_h
