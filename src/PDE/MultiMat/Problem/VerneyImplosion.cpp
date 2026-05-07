// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/VerneyImplosion.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the multi-material flow equations
  \details   This file defines a Problem policy class for the multi-material
    compressible flow equations, defined in PDE/MultiMat/DGMultiMat.hpp. See
    PDE/MultiMat/Problem.hpp for general requirements on Problem policy classes
    for MultiMat.
*/
// *****************************************************************************

#include "VerneyImplosion.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiMatProblemVerneyImplosion;

tk::InitializeFn::result_type
MultiMatProblemVerneyImplosion::initialize( ncomp_t ncomp,
                                            const std::vector< EOS >& mat_blk,
                                            tk::real x,
                                            tk::real y,
                                            tk::real z,
                                            tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
//! \details This function only initializes the Verney shell implosion
//!          problem
// *****************************************************************************
{
  auto nmat = g_inputdeck.get< eq, tag::nmat >();
  auto alphamin = g_inputdeck.get< eq, tag::min_volumefrac >();
  const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

  tk::real p;
  std::vector< tk::real > s(ncomp, 0.0);
  std::vector< tk::real > alpha(nmat, 0.0), rho(nmat, 0.0), vel(nmat, 0.0);

  // Compute radial position
  tk::real r = std::sqrt(x*x+y*y+z*z);
  // Inner region (air)
  if (r < 0.08) {
    alpha = {{ 1.0-2*alphamin, alphamin, alphamin }};
    rho = {{ 0.001, 7900.0, 1.0 }};
    vel = {{ 0.0, 0.0, 0.0 }};
    p = 101325.0;
  }
  else if (0.08 <= r && r <= 0.085) {
    alpha = {{ alphamin, 1.0-2*alphamin, alphamin }};
    rho = {{ 0.001, 7900.0, 1.0 }};
    tk::real r0 = 0.08;
    tk::real u0 = 1400.0;
    tk::real vel_mag = u0 * std::pow((r0*r0)/(r*r),2.0);
    vel = {{ -vel_mag*x/r, -vel_mag*y/r, -vel_mag*z/r }};
    p = 101325.0;
  }
  else {
    alpha = {{ alphamin, alphamin, 1.0-2*alphamin }};
    rho = {{ 0.001, 7900.0, 1.0 }};
    vel = {{ 0.0, 0.0, 0.0 }};
    p = 101325.0;
  }
  tk::real rhob = 0.0;
  for (std::size_t k=0; k<nmat; ++k) rhob += rho[k];
  s[momentumIdx(nmat, 0)] = rhob * vel[0];
  s[momentumIdx(nmat, 1)] = rhob * vel[1];
  s[momentumIdx(nmat, 2)] = rhob * vel[2];
  for (std::size_t k=0; k<nmat; ++k) {
    s[volfracIdx(nmat, k)] = alpha[k];
    s[densityIdx(nmat, k)] = s[volfracIdx(nmat, k)] * rho[k];
    std::array< std::array< tk::real, 3 >, 3 > g;
    if (solidx[k] > 0) {
      for (std::size_t i=0; i<3; ++i) {
        for (std::size_t j=0; j<3; ++j) {
          if (i==j) g[i][j] = 1.0;
          else g[i][j] = 0.0;
          s[deformIdx(nmat,solidx[k],i,j)] = g[i][j];
        }
      }
    }
    else {
      g = {{}};
    }
    // printf("k = %lu\n", k);
    // printf("alpha, rho = %e, %e\n", alpha[k], rho[k]);
    // printf("vel = %e, %e, %e\n", vel[0], vel[1], vel[2]);
    // printf("p = %e\n", p);
    s[energyIdx(nmat, k)] =
      mat_blk[k].compute< EOS::totalenergy >(
        alpha[k]*rho[k], vel[0], vel[1], vel[2],
        alpha[k]*p, alpha[k], g);
  }

  return s;
}
