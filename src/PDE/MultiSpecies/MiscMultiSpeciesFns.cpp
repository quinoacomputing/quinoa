// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/MiscMultiSpeciesFns.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Misc multi-species system functions
  \details   This file defines functions that required for multi-species
    compressible fluid dynamics.
*/
// *****************************************************************************

#include <iostream>

#include "MiscMultiSpeciesFns.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Integrate/Basis.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"
#include "EoS/GetMatProp.hpp"
#include "MultiSpecies/Mixture/Mixture.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

void initializeSpeciesEoS( std::vector< EOS >& mat_blk )
// *****************************************************************************
//  Initialize the species block with configured EOS
//! \param[in,out] mat_blk Material block that gets initialized
// *****************************************************************************
{
  // EoS initialization
  // ---------------------------------------------------------------------------
  auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();
  const auto& matprop = g_inputdeck.get< tag::material >();
  const auto& matidxmap = g_inputdeck.get< tag::matidxmap >();
  // assume only one type of species
  auto mateos = matprop[matidxmap.get< tag::eosidx >()[0]].get<tag::eos>();
  for (std::size_t k=0; k<nspec; ++k) {
    mat_blk.emplace_back(mateos, EqType::multispecies, k);
  }
}

tk::real
timeStepSizeMultiSpecies(
  const std::vector< EOS >& mat_blk,
  const std::vector< int >& esuf,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const std::size_t nelem,
  std::size_t nspec,
  const tk::Fields& U,
  const tk::Fields& P )
// *****************************************************************************
//  Time step restriction for multi species cell-centered schemes
//! \param[in] mat_blk EOS species block
//! \param[in] esuf Elements surrounding elements array
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] nelem Number of elements
//! \param[in] nspec Number of speciess in this PDE system
//! \param[in] U High-order solution vector
//! \param[in] P High-order vector of primitives
//! \return Maximum allowable time step based on cfl criterion
// *****************************************************************************
{
  const auto ndof = g_inputdeck.get< tag::ndof >();
  const auto rdof = g_inputdeck.get< tag::rdof >();
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

  tk::real u, v, w, a, vn, dSV_l, dSV_r;
  std::vector< tk::real > delt(U.nunk(), 0.0);
  std::vector< tk::real > ugp(ncomp, 0.0), pgp(nprim, 0.0);

  // compute maximum characteristic speed at all internal element faces
  for (std::size_t f=0; f<esuf.size()/2; ++f)
  {
    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    auto er = esuf[2*f+1];

    std::array< tk::real, 3 > fn;
    fn[0] = geoFace(f,1);
    fn[1] = geoFace(f,2);
    fn[2] = geoFace(f,3);

    // left element

    // Compute the basis function for the left element
    std::vector< tk::real > B_l(rdof, 0.0);
    B_l[0] = 1.0;

    // get conserved quantities
    ugp = eval_state(ncomp, rdof, ndof, el, U, B_l);
    // get primitive quantities
    pgp = eval_state(nprim, rdof, ndof, el, P, B_l);

    // initialize mixture
    Mixture mix(nspec, ugp, mat_blk);

    // mixture density
    auto rhob = mix.get_mix_density();

    // advection velocity
    u = ugp[multispecies::momentumIdx(nspec, 0)]/rhob;
    v = ugp[multispecies::momentumIdx(nspec, 1)]/rhob;
    w = ugp[multispecies::momentumIdx(nspec, 2)]/rhob;

    vn = u*geoFace(f,1) + v*geoFace(f,2) + w*geoFace(f,3);

    // acoustic speed
    a = mix.frozen_soundspeed( rhob, pgp[multispecies::temperatureIdx(nspec,0)],
      mat_blk );

    dSV_l = geoFace(f,0) * (std::fabs(vn) + a);

    // right element
    if (er > -1) {
      std::size_t eR = static_cast< std::size_t >( er );

      // Compute the basis function for the right element
      std::vector< tk::real > B_r(rdof, 0.0);
      B_r[0] = 1.0;

      // get conserved quantities
      ugp = eval_state(ncomp, rdof, ndof, eR, U, B_r);
      pgp = eval_state(nprim, rdof, ndof, eR, P, B_r);

      // initialize mixture
      Mixture mixr(nspec, ugp, mat_blk);

      // mixture density
      rhob = mixr.get_mix_density();

      // advection velocity
      u = ugp[multispecies::momentumIdx(nspec, 0)]/rhob;
      v = ugp[multispecies::momentumIdx(nspec, 1)]/rhob;
      w = ugp[multispecies::momentumIdx(nspec, 2)]/rhob;

      vn = u*geoFace(f,1) + v*geoFace(f,2) + w*geoFace(f,3);

      // acoustic speed
      a = mix.frozen_soundspeed( rhob,
        pgp[multispecies::temperatureIdx(nspec,0)], mat_blk );

      dSV_r = geoFace(f,0) * (std::fabs(vn) + a);

      delt[eR] += std::max( dSV_l, dSV_r );
    } else {
      dSV_r = dSV_l;
    }

    delt[el] += std::max( dSV_l, dSV_r );
  }

  tk::real mindt = std::numeric_limits< tk::real >::max();

  // compute allowable dt
  for (std::size_t e=0; e<nelem; ++e)
  {
    mindt = std::min( mindt, geoElem(e,0)/delt[e] );
  }

  return mindt;
}

} //inciter::
