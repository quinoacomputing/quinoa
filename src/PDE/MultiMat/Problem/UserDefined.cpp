// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/UserDefined.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a Problem policy class for the compressible flow
    equations, defined in PDE/MultiMat/MultiMat.h. See PDE/MultiMat/Problem.h
    for general requirements on Problem policy classes for MultiMat.
*/
// *****************************************************************************

#include <limits>

#include "UserDefined.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiMatProblemUserDefined;

tk::InitializeFn::result_type
MultiMatProblemUserDefined::initialize( ncomp_t ncomp,
                                        const std::vector< EOS >& mat_blk,
                                        tk::real,
                                        tk::real,
                                        tk::real,
                                        tk::real )
// *****************************************************************************
//! Set initial conditions
//! \param[in] ncomp Number of scalar components in this PDE system
//! \return Values of all components
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  tk::InitializeFn::result_type s( ncomp, 0.0 );

  auto nmat = g_inputdeck.get< tag::param, eq, tag::nmat >();
  const auto& solidx = g_inputdeck.get< tag::param, tag::multimat,
    tag::matidxmap >().template get< tag::solidx >();

  // Set background ICs
  const auto& ic = g_inputdeck.get< tag::param, eq, tag::ic >();
  const auto& bgmatid = ic.get< tag::materialid >();
  const auto& bgvelic = ic.get< tag::velocity >();
  const auto& bgpreic = ic.get< tag::pressure >();
  const auto& bgtempic = ic.get< tag::temperature >();

  Assert( bgtempic.size() > 0, "No background temperature IC" );
  Assert( bgpreic.size() > 0, "No background pressure IC" );

  auto alphamin = 1.0e-12;

  // initialize background material states
  for (std::size_t k=0; k<nmat; ++k) {
    if (k == bgmatid.at(0).at(0)-1) {
      s[volfracIdx(nmat,k)] = 1.0 - (static_cast< tk::real >(nmat-1))*alphamin;
    }
    else {
      s[volfracIdx(nmat,k)] = alphamin;
    }
  }

  tk::real u = bgvelic[0][0];
  tk::real v = bgvelic[0][1];
  tk::real w = bgvelic[0][2];

  auto rb = 0.0;
  for (std::size_t k=0; k<nmat; ++k) {
    // density
    auto rhok = mat_blk[k].compute< EOS::density >(bgpreic[0][0],
      bgtempic[0][0]);
    // partial density
    s[densityIdx(nmat,k)] = s[volfracIdx(nmat,k)] * rhok;
    // deformation gradients
    std::array< std::array< tk::real, 3 >, 3 > g;
    if (solidx[k] > 0) {
      for (std::size_t i=0; i<3; ++i) {
        for (std::size_t j=0; j<3; ++j) {
          if (i==j) g[i][j] = 1.0;
          else g[i][j] = 0.0;
          s[deformIdx(nmat,solidx[k],i,j)] = s[volfracIdx(nmat,k)]*g[i][j];
        }
      }
    }
    else {
      g = {{}};
    }
    // total specific energy
    s[energyIdx(nmat,k)] = s[volfracIdx(nmat,k)] *
      mat_blk[k].compute< EOS::totalenergy >(rhok, u, v, w, bgpreic[0][0],
      g);
    // bulk density
    rb += s[densityIdx(nmat,k)];
  }

  // bulk momentum
  s[momentumIdx(nmat,0)] = rb * u;
  s[momentumIdx(nmat,1)] = rb * v;
  s[momentumIdx(nmat,2)] = rb * w;

  if (bgpreic[0].empty() || bgtempic[0].empty())
    Throw("User must specify background pressure and temperature in IC.");

  return s;
}
