// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/NLEnergyGrowth.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file declares a problem policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problem.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowProblemNLEnergyGrowth_h
#define CompFlowProblemNLEnergyGrowth_h

#include <string>
#include <unordered_set>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "EoS/GetMatProp.hpp"
#include "EoS/EOS.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! CompFlow system of PDEs problem: nonlinear energy growth (NLEG)
//! \see Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
//!   equations with relevance to Inertial Confinement Fusion", Journal of
//!   Computational Physics 267 (2014) 196-209.
class CompFlowProblemNLEnergyGrowth {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::compflow;

    //! Compute internal energy parameter
    static tk::real hx( tk::real bx, tk::real by, tk::real bz,
                        tk::real x, tk::real y, tk::real z );

    //! Compute a power of the internal energy
    static tk::real ec( tk::real ce, tk::real kappa, tk::real t, tk::real h,
                        tk::real p );

  public:
    //! Initialize numerical solution
    static tk::InitializeFn::result_type
    initialize( ncomp_t, const std::vector< EOS >&,
                tk::real x, tk::real y, tk::real z, tk::real t );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static tk::InitializeFn::result_type
    analyticSolution( ncomp_t, const std::vector< EOS >&, tk::real x, tk::real y,
                      tk::real z, tk::real t );

    //! Compute and return source term for NLEG manufactured solution
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \param[in] t Physical time at which to evaluate the source
    //! \param[in,out] sv Source term vector
    //! \note The function signature must follow tk::SrcFn
    static tk::SrcFn::result_type
    src( ncomp_t, const std::vector< EOS >&, tk::real x,
         tk::real y, tk::real z, tk::real t, std::vector< tk::real >& sv )
    {
      Assert(sv.size() == 5, "Incorrect source vector size");
      using tag::param; using std::sin; using std::cos;
      // manufactured solution parameters
      const auto a = g_inputdeck.get< param, eq, tag::alpha >()[0];
      const auto bx = g_inputdeck.get< param, eq, tag::betax >()[0];
      const auto by = g_inputdeck.get< param, eq, tag::betay >()[0];
      const auto bz = g_inputdeck.get< param, eq, tag::betaz >()[0];
      const auto ce = g_inputdeck.get< param, eq, tag::ce >()[0];
      const auto kappa = g_inputdeck.get< param, eq, tag::kappa >()[0];
      const auto r0 = g_inputdeck.get< param, eq, tag::r0 >()[0];
      // ratio of specific heats
      const auto g = gamma< tag::compflow >(0);
      // spatial component of density field
      const auto gx = 1.0 - x*x - y*y - z*z;
      // derivative of spatial component of density field
      const std::array< tk::real, 3 > dg{{ -2.0*x, -2.0*y, -2.0*z }};
      // spatial component of energy field
      const auto h = hx( bx, by, bz, x, y, z );
      // derivative of spatial component of energy field
      std::array< tk::real, 3 >
        dh{{ -bx*M_PI*sin(bx*M_PI*x)*cos(by*M_PI*y)*cos(bz*M_PI*z),
             -by*M_PI*cos(bx*M_PI*x)*sin(by*M_PI*y)*cos(bz*M_PI*z),
             -bz*M_PI*cos(bx*M_PI*x)*cos(by*M_PI*y)*sin(bz*M_PI*z) }};
      // temporal function f and its derivative
      const auto ft = std::exp(-a*t);
      const auto dfdt = -a*ft;
      // density and its derivatives
      const auto rho = r0 + ft*gx;
      const std::array< tk::real, 3 > drdx{{ ft*dg[0], ft*dg[1], ft*dg[2] }};
      const auto drdt = gx*dfdt;
      // internal energy and its derivatives
      const auto ie = ec( ce, kappa, t, h, -1.0/3.0 );
      const std::array< tk::real, 3 > dedx{{
        2.0*std::pow(ie,4.0)*kappa*h*dh[0]*t,
        2.0*std::pow(ie,4.0)*kappa*h*dh[1]*t,
        2.0*std::pow(ie,4.0)*kappa*h*dh[2]*t }};
      const auto dedt = kappa*h*h*std::pow(ie,4.0);
      // density source
      sv[0] = drdt;
      // momentum source
      sv[1] = (g-1.0)*(rho*dedx[0] + ie*drdx[0]);
      sv[2] = (g-1.0)*(rho*dedx[1] + ie*drdx[1]);
      sv[3] = (g-1.0)*(rho*dedx[2] + ie*drdx[2]);
      // energy source
      sv[4] = rho*dedt + ie*drdt;
    }

    //! Return analytic field names to be output to file
    std::vector< std::string > analyticFieldNames( ncomp_t ) const;

    //! Return names of integral variables to be output to diagnostics file
    std::vector< std::string > names( ncomp_t ) const;

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::NL_ENERGY_GROWTH; }
};

} // inciter::

#endif // CompFlowProblemNLEnergyGrowth_h
