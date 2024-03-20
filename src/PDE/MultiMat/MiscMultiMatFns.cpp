// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/MiscMultiMatFns.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Misc multi-material system functions
  \details   This file defines functions that required for multi-material
    compressible hydrodynamics.
*/
// *****************************************************************************

#include <iostream>

#include "MiscMultiMatFns.hpp"
#include "Inciter/InputDeck/New2InputDeck.hpp"
#include "Integrate/Basis.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::New2InputDeck g_newinputdeck;

void initializeMaterialEoS( std::vector< EOS >& mat_blk )
// *****************************************************************************
//  Initialize the material block with configured EOS
//! \param[in,out] mat_blk Material block that gets initialized
// *****************************************************************************
{
  // EoS initialization
  auto nmat = g_newinputdeck.get< newtag::multimat, newtag::nmat >();
  const auto& matprop = g_newinputdeck.get< newtag::material >();
  const auto& matidxmap = g_newinputdeck.get< newtag::matidxmap >();
  for (std::size_t k=0; k<nmat; ++k) {
    auto mateos = matprop[matidxmap.get< newtag::eosidx >()[k]].get<newtag::eos>();
    mat_blk.emplace_back(mateos, EqType::multimat, k);
  }
}

bool
cleanTraceMultiMat(
  tk::real t,
  std::size_t nelem,
  const std::vector< EOS >& mat_blk,
  const tk::Fields& geoElem,
  std::size_t nmat,
  tk::Fields& U,
  tk::Fields& P )
// *****************************************************************************
//  Clean up the state of trace materials for multi-material PDE system
//! \param[in] t Physical time
//! \param[in] nelem Number of elements
//! \param[in] mat_blk EOS material block
//! \param[in] geoElem Element geometry array
//! \param[in] nmat Number of materials in this PDE system
//! \param[in/out] U High-order solution vector which gets modified
//! \param[in/out] P High-order vector of primitives which gets modified
//! \return Boolean indicating if an unphysical material state was found
// *****************************************************************************
{
  const auto ndof = g_newinputdeck.get< newtag::ndof >();
  const auto rdof = g_newinputdeck.get< newtag::rdof >();
  std::size_t ncomp = U.nprop()/rdof;
  auto al_eps = 1.0e-02;
  auto neg_density = false;

  std::vector< tk::real > ugp(ncomp, 0.0);

  for (std::size_t e=0; e<nelem; ++e)
  {
    // find material in largest quantity, and determine if pressure
    // relaxation is needed. If it is, determine materials that need
    // relaxation, and the total volume of these materials.
    std::vector< int > relaxInd(nmat, 0);
    auto almax(0.0), relaxVol(0.0);
    std::size_t kmax = 0;
    for (std::size_t k=0; k<nmat; ++k)
    {
      auto al = U(e, volfracDofIdx(nmat, k, rdof, 0));
      if (al > almax)
      {
        almax = al;
        kmax = k;
      }
      else if (al < al_eps)
      {
        relaxInd[k] = 1;
        relaxVol += al;
      }
    }
    relaxInd[kmax] = 1;
    relaxVol += almax;

    // get conserved quantities
    std::vector< tk::real > B(rdof, 0.0);
    B[0] = 1.0;
    ugp = eval_state(ncomp, rdof, ndof, e, U, B);

    auto u = P(e, velocityDofIdx(nmat, 0, rdof, 0));
    auto v = P(e, velocityDofIdx(nmat, 1, rdof, 0));
    auto w = P(e, velocityDofIdx(nmat, 2, rdof, 0));
    auto pmax = P(e, pressureDofIdx(nmat, kmax, rdof, 0))/almax;
    auto agmax = getDeformGrad(nmat, kmax, ugp);
    auto tmax = mat_blk[kmax].compute< EOS::temperature >(
      U(e, densityDofIdx(nmat, kmax, rdof, 0)), u, v, w,
      U(e, energyDofIdx(nmat, kmax, rdof, 0)), almax, agmax );

    tk::real p_target(0.0), d_al(0.0), d_arE(0.0);
    //// get equilibrium pressure
    //std::vector< tk::real > kmat(nmat, 0.0);
    //tk::real ratio(0.0);
    //for (std::size_t k=0; k<nmat; ++k)
    //{
    //  auto arhok = U(e, densityDofIdx(nmat,k,rdof,0));
    //  auto alk = U(e, volfracDofIdx(nmat,k,rdof,0));
    //  auto apk = P(e, pressureDofIdx(nmat,k,rdof,0)3);
    //  auto ak = mat_blk[k].compute< EOS::soundspeed >(arhok, apk, alk, k);
    //  kmat[k] = arhok * ak * ak;

    //  p_target += alk * apk / kmat[k];
    //  ratio += alk * alk / kmat[k];
    //}
    //p_target /= ratio;
    //p_target = std::max(p_target, 1e-14);
    p_target = std::max(pmax, 1e-14);

    // 1. Correct minority materials and store volume/energy changes
    for (std::size_t k=0; k<nmat; ++k)
    {
      auto alk = U(e, volfracDofIdx(nmat, k, rdof, 0));
      auto pk = P(e, pressureDofIdx(nmat, k, rdof, 0)) / alk;
      // for positive volume fractions
      if (matExists(alk))
      {
        // check if volume fraction is lesser than threshold (al_eps) and
        // if the material (effective) pressure is negative. If either of
        // these conditions is true, perform pressure relaxation.
        if ((alk < al_eps) ||
          (pk < mat_blk[k].compute< EOS::min_eff_pressure >(1e-12,
          U(e, densityDofIdx(nmat, k, rdof, 0)), alk))
          /*&& (std::fabs((pk-pmax)/pmax) > 1e-08)*/)
        {
          //auto gk = gamma< newtag::multimat >(0, k);

          tk::real alk_new(0.0);
          //// volume change based on polytropic expansion/isentropic compression
          //if (pk > p_target)
          //{
          //  alk_new = std::pow((pk/p_target), (1.0/gk)) * alk;
          //}
          //else
          //{
          //  auto arhok = U(e, densityDofIdx(nmat, k, rdof, 0));
          //  auto ck = eos_soundspeed< newtag::multimat >(0, arhok, alk*pk,
          //    alk, k);
          //  auto kk = arhok * ck * ck;
          //  alk_new = alk - (alk*alk/kk) * (p_target-pk);
          //}
          alk_new = alk;

          // determine target relaxation pressure
          auto prelax = mat_blk[k].compute< EOS::min_eff_pressure >(1e-10,
            U(e, densityDofIdx(nmat, k, rdof, 0)), alk_new);
          prelax = std::max(prelax, p_target);

          // energy change
          auto rhomat = U(e, densityDofIdx(nmat, k, rdof, 0))
            / alk_new;
          auto gmat = getDeformGrad(nmat, k, ugp);
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              gmat[i][j] /= alk_new;
          auto rhoEmat = mat_blk[k].compute< EOS::totalenergy >(rhomat, u, v, w,
            prelax, gmat);

          // volume-fraction and total energy flux into majority material
          d_al += (alk - alk_new);
          d_arE += (U(e, energyDofIdx(nmat, k, rdof, 0))
            - alk_new * rhoEmat);

          // update state of trace material
          U(e, volfracDofIdx(nmat, k, rdof, 0)) = alk_new;
          U(e, energyDofIdx(nmat, k, rdof, 0)) = alk_new*rhoEmat;
          P(e, pressureDofIdx(nmat, k, rdof, 0)) = alk_new*prelax;
        }
      }
      // check for unbounded volume fractions
      else if (alk < 0.0 || !std::isfinite(alk))
      {
        auto rhok = mat_blk[k].compute< EOS::density >(p_target,
          std::max(1e-8,tmax));
        if (std::isfinite(alk)) d_al += (alk - 1e-14);
        // update state of trace material
        U(e, volfracDofIdx(nmat, k, rdof, 0)) = 1e-14;
        U(e, densityDofIdx(nmat, k, rdof, 0)) = 1e-14 * rhok;
        auto gmax = agmax;
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            gmax[i][j] /= almax;
        U(e, energyDofIdx(nmat, k, rdof, 0)) = 1e-14
          * mat_blk[k].compute< EOS::totalenergy >(rhok, u, v, w, p_target,
          gmax);
        P(e, pressureDofIdx(nmat, k, rdof, 0)) = 1e-14 *
          p_target;
        for (std::size_t i=1; i<rdof; ++i) {
          U(e, volfracDofIdx(nmat, k, rdof, i)) = 0.0;
          U(e, densityDofIdx(nmat, k, rdof, i)) = 0.0;
          U(e, energyDofIdx(nmat, k, rdof, i)) = 0.0;
          P(e, pressureDofIdx(nmat, k, rdof, i)) = 0.0;
        }
      }
      else {
        // determine target relaxation pressure
        auto prelax = mat_blk[k].compute< EOS::min_eff_pressure >(1e-10,
          U(e, densityDofIdx(nmat, k, rdof, 0)), alk);
        prelax = std::max(prelax, p_target);
        auto rhok = U(e, densityDofIdx(nmat, k, rdof, 0)) / alk;
        auto gk = getDeformGrad(nmat, k, ugp);
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            gk[i][j] /= alk;
        // update state of trace material
        U(e, energyDofIdx(nmat, k, rdof, 0)) = alk
          * mat_blk[k].compute< EOS::totalenergy >( rhok, u, v, w, prelax, gk );
        P(e, pressureDofIdx(nmat, k, rdof, 0)) = alk *
          prelax;
        for (std::size_t i=1; i<rdof; ++i) {
          U(e, energyDofIdx(nmat, k, rdof, i)) = 0.0;
          P(e, pressureDofIdx(nmat, k, rdof, i)) = 0.0;
        }
      }
    }

    // 2. Based on volume change in majority material, compute energy change
    //auto gmax = gamma< newtag::multimat >(0, kmax);
    //auto pmax_new = pmax * std::pow(almax/(almax+d_al), gmax);
    //auto rhomax_new = U(e, densityDofIdx(nmat, kmax, rdof, 0))
    //  / (almax+d_al);
    //auto rhoEmax_new = eos_totalenergy< newtag::multimat >(0, rhomax_new, u,
    //  v, w, pmax_new, kmax);
    //auto d_arEmax_new = (almax+d_al) * rhoEmax_new
    //  - U(e, energyDofIdx(nmat, kmax, rdof, 0));

    U(e, volfracDofIdx(nmat, kmax, rdof, 0)) += d_al;
    //U(e, energyDofIdx(nmat, kmax, rdof, 0)) += d_arEmax_new;

    // 2. Flux energy change into majority material
    U(e, energyDofIdx(nmat, kmax, rdof, 0)) += d_arE;
    P(e, pressureDofIdx(nmat, kmax, rdof, 0)) =
      mat_blk[kmax].compute< EOS::pressure >(
      U(e, densityDofIdx(nmat, kmax, rdof, 0)), u, v, w,
      U(e, energyDofIdx(nmat, kmax, rdof, 0)),
      U(e, volfracDofIdx(nmat, kmax, rdof, 0)), kmax, agmax );

    // enforce unit sum of volume fractions
    auto alsum = 0.0;
    for (std::size_t k=0; k<nmat; ++k)
      alsum += U(e, volfracDofIdx(nmat, k, rdof, 0));

    for (std::size_t k=0; k<nmat; ++k) {
      U(e, volfracDofIdx(nmat, k, rdof, 0)) /= alsum;
      U(e, densityDofIdx(nmat, k, rdof, 0)) /= alsum;
      U(e, energyDofIdx(nmat, k, rdof, 0)) /= alsum;
      P(e, pressureDofIdx(nmat, k, rdof, 0)) /= alsum;
    }

    //// bulk quantities
    //auto rhoEb(0.0), rhob(0.0), volb(0.0);
    //for (std::size_t k=0; k<nmat; ++k)
    //{
    //  if (relaxInd[k] > 0.0)
    //  {
    //    rhoEb += U(e, energyDofIdx(nmat,k,rdof,0));
    //    volb += U(e, volfracDofIdx(nmat,k,rdof,0));
    //    rhob += U(e, densityDofIdx(nmat,k,rdof,0));
    //  }
    //}

    //// 2. find mixture-pressure
    //tk::real pmix(0.0), den(0.0);
    //pmix = rhoEb - 0.5*rhob*(u*u+v*v+w*w);
    //for (std::size_t k=0; k<nmat; ++k)
    //{
    //  auto gk = gamma< newtag::multimat >(0, k);
    //  auto Pck = pstiff< newtag::multimat >(0, k);

    //  pmix -= U(e, volfracDofIdx(nmat,k,rdof,0)) * gk * Pck *
    //    relaxInd[k] / (gk-1.0);
    //  den += U(e, volfracDofIdx(nmat,k,rdof,0)) * relaxInd[k]
    //    / (gk-1.0);
    //}
    //pmix /= den;

    //// 3. correct energies
    //for (std::size_t k=0; k<nmat; ++k)
    //{
    //  if (relaxInd[k] > 0.0)
    //  {
    //    auto alk_new = U(e, volfracDofIdx(nmat,k,rdof,0));
    //    U(e, energyDofIdx(nmat,k,rdof,0)) = alk_new *
    //      eos_totalenergy< newtag::multimat >(0, rhomat[k], u, v, w, pmix,
    //      k);
    //    P(e, pressureDofIdx(nmat, k, rdof, 0)) = alk_new * pmix;
    //  }
    //}

    pmax = P(e, pressureDofIdx(nmat, kmax, rdof, 0)) /
      U(e, volfracDofIdx(nmat, kmax, rdof, 0));

    // check for unphysical state
    for (std::size_t k=0; k<nmat; ++k)
    {
      auto alpha = U(e, volfracDofIdx(nmat, k, rdof, 0));
      auto arho = U(e, densityDofIdx(nmat, k, rdof, 0));
      auto apr = P(e, pressureDofIdx(nmat, k, rdof, 0));

      // lambda for screen outputs
      auto screenout = [&]()
      {
        std::cout << "Physical time:     " << t << std::endl;
        std::cout << "Element centroid:  " << geoElem(e,1) << ", "
          << geoElem(e,2) << ", " << geoElem(e,3) << std::endl;
        std::cout << "Material-id:       " << k << std::endl;
        std::cout << "Volume-fraction:   " << alpha << std::endl;
        std::cout << "Partial density:   " << arho << std::endl;
        std::cout << "Partial pressure:  " << apr << std::endl;
        std::cout << "Major pressure:    " << pmax << " (mat " << kmax << ")"
          << std::endl;
        std::cout << "Major temperature: " << tmax << " (mat " << kmax << ")"
          << std::endl;
        std::cout << "Velocity:          " << u << ", " << v << ", " << w
          << std::endl;
      };

      if (arho < 0.0)
      {
        neg_density = true;
        screenout();
      }
    }
  }
  return neg_density;
}

tk::real
timeStepSizeMultiMat(
  const std::vector< EOS >& mat_blk,
  const std::vector< int >& esuf,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const std::size_t nelem,
  std::size_t nmat,
  const tk::Fields& U,
  const tk::Fields& P )
// *****************************************************************************
//  Time step restriction for multi material cell-centered schemes
//! \param[in] mat_blk EOS material block
//! \param[in] esuf Elements surrounding elements array
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] nelem Number of elements
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] U High-order solution vector
//! \param[in] P High-order vector of primitives
//! \return Maximum allowable time step based on cfl criterion
// *****************************************************************************
{
  const auto ndof = g_newinputdeck.get< newtag::ndof >();
  const auto rdof = g_newinputdeck.get< newtag::rdof >();
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

    // advection velocity
    u = pgp[velocityIdx(nmat, 0)];
    v = pgp[velocityIdx(nmat, 1)];
    w = pgp[velocityIdx(nmat, 2)];

    vn = u*geoFace(f,1) + v*geoFace(f,2) + w*geoFace(f,3);

    // acoustic speed
    a = 0.0;
    for (std::size_t k=0; k<nmat; ++k)
    {
      if (ugp[volfracIdx(nmat, k)] > 1.0e-04) {
        auto gk = getDeformGrad(nmat, k, ugp);
        gk = tk::rotateTensor(gk, fn);
        a = std::max( a, mat_blk[k].compute< EOS::soundspeed >(
          ugp[densityIdx(nmat, k)],
          pgp[pressureIdx(nmat, k)], ugp[volfracIdx(nmat, k)], k, gk ) );
      }
    }

    dSV_l = geoFace(f,0) * (std::fabs(vn) + a);

    // right element
    if (er > -1) {
      std::size_t eR = static_cast< std::size_t >( er );

      // Compute the basis function for the right element
      std::vector< tk::real > B_r(rdof, 0.0);
      B_r[0] = 1.0;

      // get conserved quantities
      ugp = eval_state( ncomp, rdof, ndof, eR, U, B_r);
      // get primitive quantities
      pgp = eval_state( nprim, rdof, ndof, eR, P, B_r);

      // advection velocity
      u = pgp[velocityIdx(nmat, 0)];
      v = pgp[velocityIdx(nmat, 1)];
      w = pgp[velocityIdx(nmat, 2)];

      vn = u*geoFace(f,1) + v*geoFace(f,2) + w*geoFace(f,3);

      // acoustic speed
      a = 0.0;
      for (std::size_t k=0; k<nmat; ++k)
      {
        if (ugp[volfracIdx(nmat, k)] > 1.0e-04) {
          auto gk = getDeformGrad(nmat, k, ugp);
          gk = tk::rotateTensor(gk, fn);
          a = std::max( a, mat_blk[k].compute< EOS::soundspeed >(
            ugp[densityIdx(nmat, k)],
            pgp[pressureIdx(nmat, k)], ugp[volfracIdx(nmat, k)], k, gk ) );
        }
      }

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

tk::real
timeStepSizeMultiMatFV(
  const std::vector< EOS >& mat_blk,
  const tk::Fields& geoElem,
  std::size_t nelem,
  std::size_t nmat,
  const tk::Fields& U,
  const tk::Fields& P,
  std::vector< tk::real >& local_dte )
// *****************************************************************************
//  Time step restriction for multi material cell-centered FV scheme
//! \param[in] mat_blk Material EOS block
//! \param[in] geoElem Element geometry array
//! \param[in] nelem Number of elements
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] U High-order solution vector
//! \param[in] P High-order vector of primitives
//! \param[in,out] local_dte Time step size for each element (for local
//!   time stepping)
//! \return Maximum allowable time step based on cfl criterion
// *****************************************************************************
{
  const auto ndof = g_newinputdeck.get< newtag::ndof >();
  const auto rdof = g_newinputdeck.get< newtag::rdof >();
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

  std::vector< tk::real > ugp(ncomp, 0.0), pgp(nprim, 0.0);
  tk::real mindt = std::numeric_limits< tk::real >::max();

  // loop over all elements
  for (std::size_t e=0; e<nelem; ++e)
  {
    // basis function at centroid
    std::vector< tk::real > B(rdof, 0.0);
    B[0] = 1.0;

    // get conserved quantities
    ugp = eval_state(ncomp, rdof, ndof, e, U, B);
    // get primitive quantities
    pgp = eval_state(nprim, rdof, ndof, e, P, B);

    // magnitude of advection velocity
    auto u = pgp[velocityIdx(nmat, 0)];
    auto v = pgp[velocityIdx(nmat, 1)];
    auto w = pgp[velocityIdx(nmat, 2)];
    auto vmag = std::sqrt(tk::dot({u, v, w}, {u, v, w}));

    // acoustic speed
    tk::real a = 0.0;
    for (std::size_t k=0; k<nmat; ++k)
    {
      if (ugp[volfracIdx(nmat, k)] > 1.0e-04) {
        a = std::max( a, mat_blk[k].compute< EOS::soundspeed >(
          ugp[densityIdx(nmat, k)], pgp[pressureIdx(nmat, k)],
          ugp[volfracIdx(nmat, k)], k ) );
      }
    }

    // characteristic wave speed
    auto v_char = vmag + a;

    // characteristic length (radius of insphere)
    auto dx = std::min(std::cbrt(geoElem(e,0)), geoElem(e,4))
      /std::sqrt(24.0);

    // element dt
    local_dte[e] = dx/v_char;
    mindt = std::min(mindt, local_dte[e]);
  }

  return mindt;
}

std::array< std::array< tk::real, 3 >, 3 >
getDeformGrad(
  std::size_t nmat,
  std::size_t k,
  const std::vector< tk::real >& state )
// *****************************************************************************
//  Get the inverse deformation gradient tensor for a material at given location
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] k Material id whose deformation gradient is required
//! \param[in] state State vector at location
//! \return Inverse deformation gradient tensor (alpha_k * g_k) of material k
// *****************************************************************************
{
  const auto& solidx = g_newinputdeck.get< newtag::matidxmap, newtag::solidx >();
  std::array< std::array< tk::real, 3 >, 3 > agk;

  if (solidx[k] > 0) {
    // deformation gradient for solids
    for (std::size_t i=0; i<3; ++i) {
      for (std::size_t j=0; j<3; ++j)
        agk[i][j] = state[deformIdx(nmat,solidx[k],i,j)];
    }
  }
  else {
    // empty vector for fluids
    agk = {{}};
  }

  return agk;
}

std::array< std::array< tk::real, 3 >, 3 >
getCauchyStress(
  std::size_t nmat,
  std::size_t k,
  std::size_t ncomp,
  const std::vector< tk::real >& state )
// *****************************************************************************
//  Get the elastic Cauchy stress tensor for a material at given location
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] k Material id whose deformation gradient is required
//! \param[in] ncomp Number of components in the PDE system
//! \param[in] state State vector at location
//! \return Elastic Cauchy stress tensor (alpha * \sigma_ij) of material k
// *****************************************************************************
{
  const auto& solidx = g_newinputdeck.get< newtag::matidxmap, newtag::solidx >();
  std::array< std::array< tk::real, 3 >, 3 >
    asigk{{ {{0,0,0}}, {{0,0,0}}, {{0,0,0}} }};

  // elastic Cauchy stress for solids
  if (solidx[k] > 0) {
    for (std::size_t i=0; i<3; ++i) {
      for (std::size_t j=0; j<3; ++j)
        asigk[i][j] = state[ncomp +
          stressIdx(nmat, solidx[k], stressCmp[i][j])];
    }
  }

  return asigk;
}

bool
haveSolid(
  std::size_t nmat,
  const std::vector< std::size_t >& solidx )
// *****************************************************************************
//  Check whether we have solid materials in our problem
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] solidx Material index indicator
//! \return true if we have at least one solid, false otherwise.
// *****************************************************************************
{
  bool haveSolid = false;
  for (std::size_t k=0; k<nmat; ++k)
    if (solidx[k] > 0) haveSolid = true;

  return haveSolid;
}

std::size_t numSolids(
  std::size_t nmat,
  const std::vector< std::size_t >& solidx )
// *****************************************************************************
//  Count total number of solid materials in the problem
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] solidx Material index indicator
//! \return Total number of solid materials in the problem
// *****************************************************************************
{
  // count number of solid materials
  std::size_t nsld(0);
  for (std::size_t k=0; k<nmat; ++k)
    if (solidx[k] > 0) ++nsld;

  return nsld;
}

} //inciter::
