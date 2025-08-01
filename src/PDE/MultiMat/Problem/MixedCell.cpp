// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/MixedCell.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for a mixed cell with different pressures
  \details   This file defines a Problem policy class for the compressible flow
    equations, defined in PDE/MultiMat/MultiMat.h. See PDE/MultiMat/Problem.h
    for general requirements on Problem policy classes for MultiMat.
*/
// *****************************************************************************

#include <limits>

#include "MixedCell.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiMatProblemMixedCell;

tk::InitializeFn::result_type
MultiMatProblemMixedCell::initialize( ncomp_t ncomp,
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

  auto nmat = g_inputdeck.get< eq, tag::nmat >();
  auto alphamin = g_inputdeck.get< eq, tag::min_volumefrac >();
  const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

  // Set background ICs
  const auto& ic = g_inputdeck.get< tag::ic >();
  const auto& bgmatid = ic.get< tag::materialid >();
  const auto& bgvelic = ic.get< tag::velocity >();
  const auto& bgpreic = ic.get< tag::pressure >();
  const auto& bgtempic = ic.get< tag::temperature >();

  std::vector< tk::real > pressure(nmat, 0.0);

  if (bgtempic < 1e-12) Throw( "No background temperature IC" );
  if (bgpreic < 1e-12) Throw( "No background pressure IC" );

  // initialize background material states
  for (std::size_t k=0; k<nmat; ++k) {
    s[volfracIdx(nmat,k)] = 1.0/nmat;
    pressure[k] = bgpreic*std::pow(10.0,k);
  }

  tk::real u = bgvelic[0];
  tk::real v = bgvelic[1];
  tk::real w = bgvelic[2];

  auto rb = 0.0;
  for (std::size_t k=0; k<nmat; ++k) {
    // density
    auto rhok = mat_blk[k].compute< EOS::density >(bgpreic, bgtempic);
    // partial density
    s[densityIdx(nmat,k)] = s[volfracIdx(nmat,k)] * rhok;
    // deformation gradients
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
    // total specific energy
    s[energyIdx(nmat,k)] =
      mat_blk[k].compute< EOS::totalenergy >(s[volfracIdx(nmat,k)]*rhok, u, v,
      w, s[volfracIdx(nmat,k)]*pressure[k], s[volfracIdx(nmat,k)], g);
    // bulk density
    rb += s[densityIdx(nmat,k)];
  }

  // bulk momentum
  s[momentumIdx(nmat,0)] = rb * u;
  s[momentumIdx(nmat,1)] = rb * v;
  s[momentumIdx(nmat,2)] = rb * w;

  if (bgpreic< 1e-12 || bgtempic< 1e-12)
    Throw("User must specify background pressure and temperature in IC.");

  return s;
}
