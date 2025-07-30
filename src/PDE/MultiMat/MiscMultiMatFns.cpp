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
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Integrate/Basis.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "EoS/GetMatProp.hpp"
#include "Kokkos_Core.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

void initializeMaterialEoS( std::vector< EOS >& mat_blk )
// *****************************************************************************
//  Initialize the material block with configured EOS
//! \param[in,out] mat_blk Material block that gets initialized
// *****************************************************************************
{
  // EoS initialization
  // ---------------------------------------------------------------------------
  auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
  const auto& matprop = g_inputdeck.get< tag::material >();
  const auto& matidxmap = g_inputdeck.get< tag::matidxmap >();
  for (std::size_t k=0; k<nmat; ++k) {
    auto mateos = matprop[matidxmap.get< tag::eosidx >()[k]].get<tag::eos>();
    mat_blk.emplace_back(mateos, EqType::multimat, k);
  }

  // set rho0 for all materials
  // ---------------------------------------------------------------------------
  std::vector< tk::real > rho0mat( nmat, 0.0 );
  const auto& ic = g_inputdeck.get< tag::ic >();
  const auto& icbox = ic.get< tag::box >();
  const auto& icmbk = ic.get< tag::meshblock >();
  // Get background properties
  std::size_t k = ic.get< tag::materialid >() - 1;
  tk::real pr = ic.get< tag::pressure >();
  tk::real tmp = ic.get< tag::temperature >();
  rho0mat[k] = mat_blk[k].compute< EOS::density >(pr, tmp);

  // Check inside used defined box
  if (!icbox.empty())
  {
    for (const auto& b : icbox) {   // for all boxes
      k = b.get< tag::materialid >() - 1;
      auto boxmas = b.get< tag::mass >();
      // if mass and volume are given, compute density as mass/volume
      if (boxmas > 0.0) {
        std::vector< tk::real > box
          { b.get< tag::xmin >(), b.get< tag::xmax >(),
            b.get< tag::ymin >(), b.get< tag::ymax >(),
            b.get< tag::zmin >(), b.get< tag::zmax >() };
        auto V_ex = (box[1]-box[0]) * (box[3]-box[2]) * (box[5]-box[4]);
        rho0mat[k] = boxmas / V_ex;
      }
      // else compute density from eos
      else {
        pr = b.get< tag::pressure >();
        tmp = b.get< tag::temperature >();
        rho0mat[k] = mat_blk[k].compute< EOS::density >(pr, tmp);
      }
    }
  }

  // Check inside user-specified mesh blocks
  if (!icmbk.empty())
  {
    for (const auto& b : icmbk) { // for all blocks
      k = b.get< tag::materialid >() - 1;
      auto boxmas = b.get< tag::mass >();
      // if mass and volume are given, compute density as mass/volume
      if (boxmas > 0.0) {
        auto V_ex = b.get< tag::volume >();
        rho0mat[k] = boxmas / V_ex;
      }
      // else compute density from eos
      else {
        pr = b.get< tag::pressure >();
        tmp = b.get< tag::temperature >();
        rho0mat[k] = mat_blk[k].compute< EOS::density >(pr, tmp);
      }
    }
  }

  // Store computed rho0's in the EOS-block
  for (std::size_t i=0; i<nmat; ++i) {
    mat_blk[i].set< EOS::setRho0 >(rho0mat[i]);
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
  const auto ndof = g_inputdeck.get< tag::ndof >();
  const auto rdof = g_inputdeck.get< tag::rdof >();
  const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();
  std::size_t ncomp = U.nprop()/rdof;
  auto neg_density = false;

  std::vector< tk::real > ugp(ncomp, 0.0);

  for (std::size_t e=0; e<nelem; ++e)
  {
    // find material in largest quantity.
    auto almax(0.0);
    std::size_t kmax = 0;
    for (std::size_t k=0; k<nmat; ++k)
    {
      auto al = U(e, volfracDofIdx(nmat, k, rdof, 0));
      if (al > almax)
      {
        almax = al;
        kmax = k;
      }
    }

    // get conserved quantities
    std::vector< tk::real > B(rdof, 0.0);
    B[0] = 1.0;
    ugp = eval_state(ncomp, rdof, ndof, e, U, B);

    auto u = P(e, velocityDofIdx(nmat, 0, rdof, 0));
    auto v = P(e, velocityDofIdx(nmat, 1, rdof, 0));
    auto w = P(e, velocityDofIdx(nmat, 2, rdof, 0));
    auto pmax = P(e, pressureDofIdx(nmat, kmax, rdof, 0))/almax;
    auto gmax = getDeformGrad(nmat, kmax, ugp);
    auto tmax = mat_blk[kmax].compute< EOS::temperature >(
      U(e, densityDofIdx(nmat, kmax, rdof, 0)), u, v, w,
      U(e, energyDofIdx(nmat, kmax, rdof, 0)), almax, gmax );

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
      if (solidx[k] == 0 && solidx[kmax] == 0 && matExists(alk))
      {
        // check if volume fraction is lesser than threshold (volfracPRelaxLim)
        // and if the material (effective) pressure is negative. If either of
        // these conditions is true, perform pressure relaxation.
        if ((alk < volfracPRelaxLim()) ||
          (pk < mat_blk[k].compute< EOS::min_eff_pressure >(1e-12,
          U(e, densityDofIdx(nmat, k, rdof, 0)), alk))
          /*&& (std::fabs((pk-pmax)/pmax) > 1e-08)*/)
        {
          // determine target relaxation pressure
          auto prelax = mat_blk[k].compute< EOS::min_eff_pressure >(1e-10,
            U(e, densityDofIdx(nmat, k, rdof, 0)), alk);
          prelax = std::max(prelax, p_target);

          // energy change
          auto arhomat = U(e, densityDofIdx(nmat, k, rdof, 0));
          auto gmat = getDeformGrad(nmat, k, ugp);
          auto arhoEmat = mat_blk[k].compute< EOS::totalenergy >(arhomat, u, v, w,
            alk*prelax, alk, gmat);

          // total energy flux into majority material
          d_arE += (U(e, energyDofIdx(nmat, k, rdof, 0))
            - arhoEmat);

          // update state of trace material
          U(e, volfracDofIdx(nmat, k, rdof, 0)) = alk;
          U(e, energyDofIdx(nmat, k, rdof, 0)) = arhoEmat;
          P(e, pressureDofIdx(nmat, k, rdof, 0)) = alk*prelax;
        }
      }
      else if (!matExists(alk)) {  // condition so that else-branch not exec'ed for solids
        // determine target relaxation pressure
        auto prelax = mat_blk[k].compute< EOS::min_eff_pressure >(1e-10,
          U(e, densityDofIdx(nmat, k, rdof, 0)), alk);
        prelax = std::max(prelax, p_target);
        auto arhok = U(e, densityDofIdx(nmat, k, rdof, 0));
        auto gk = std::array< std::array< tk::real, 3 >, 3 >
          {{ {{1, 0, 0}},
             {{0, 1, 0}},
             {{0, 0, 1}} }};
        // update state of trace material
        U(e, energyDofIdx(nmat, k, rdof, 0)) =
          mat_blk[k].compute< EOS::totalenergy >( arhok, u, v, w, alk*prelax,
          alk, gk );
        P(e, pressureDofIdx(nmat, k, rdof, 0)) = alk *
          prelax;
        resetSolidTensors(nmat, k, e, U, P);
        for (std::size_t i=1; i<rdof; ++i) {
          U(e, energyDofIdx(nmat, k, rdof, i)) = 0.0;
          P(e, pressureDofIdx(nmat, k, rdof, i)) = 0.0;
        }
      }
    }

    U(e, volfracDofIdx(nmat, kmax, rdof, 0)) += d_al;

    // 2. Flux energy change into majority material
    U(e, energyDofIdx(nmat, kmax, rdof, 0)) += d_arE;
    P(e, pressureDofIdx(nmat, kmax, rdof, 0)) =
      mat_blk[kmax].compute< EOS::pressure >(
      U(e, densityDofIdx(nmat, kmax, rdof, 0)), u, v, w,
      U(e, energyDofIdx(nmat, kmax, rdof, 0)),
      U(e, volfracDofIdx(nmat, kmax, rdof, 0)), kmax, gmax );

    // 3. enforce unit sum of volume fractions
    auto alsum = 0.0;
    for (std::size_t k=0; k<nmat; ++k)
      alsum += U(e, volfracDofIdx(nmat, k, rdof, 0));

    for (std::size_t k=0; k<nmat; ++k) {
      U(e, volfracDofIdx(nmat, k, rdof, 0)) /= alsum;
      U(e, densityDofIdx(nmat, k, rdof, 0)) /= alsum;
      U(e, energyDofIdx(nmat, k, rdof, 0)) /= alsum;
      P(e, pressureDofIdx(nmat, k, rdof, 0)) /= alsum;
      if (solidx[k] > 0) {
        for (std::size_t i=0; i<6; ++i)
          P(e, stressDofIdx(nmat, solidx[k], i, rdof, 0)) /= alsum;
      }
    }

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
  const auto ndof = g_inputdeck.get< tag::ndof >();
  const auto rdof = g_inputdeck.get< tag::rdof >();
  const auto use_mass_avg =
    g_inputdeck.get< tag::multimat, tag::dt_sos_massavg >();
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

    auto vmag = std::sqrt(u*u + v*v + w*w); 

    // acoustic speed
    tk::real a = 0.0;
    tk::real mixtureDensity = 0.0;
    for (std::size_t k=0; k<nmat; ++k)
    {
      if (use_mass_avg > 0)
      {
        // mass averaging SoS
        a += ugp[densityIdx(nmat,k)]*mat_blk[k].compute< EOS::soundspeed >(
          ugp[densityIdx(nmat, k)], pgp[pressureIdx(nmat, k)],
          ugp[volfracIdx(nmat, k)], k );

        mixtureDensity += ugp[densityIdx(nmat,k)];
      }
      else
      {
        if (ugp[volfracIdx(nmat, k)] > 1.0e-04)
        {
          a = std::max( a, mat_blk[k].compute< EOS::soundspeed >(
            ugp[densityIdx(nmat, k)], pgp[pressureIdx(nmat, k)],
            ugp[volfracIdx(nmat, k)], k ) );
        }
      }
    }
    if (use_mass_avg > 0) a /= mixtureDensity;

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

tk::real
timeStepSizeViscousFV(
  const tk::Fields& geoElem,
  std::size_t nelem,
  std::size_t nmat,
  const tk::Fields& U )
// *****************************************************************************
//  Compute the time step size restriction based on this physics
//! \param[in] geoElem Element geometry array
//! \param[in] nelem Number of elements
//! \param[in] nmat Number of materials
//! \param[in] U Solution vector
//! \return Maximum allowable time step based on viscosity
// *****************************************************************************
{
  const auto& ndof = g_inputdeck.get< tag::ndof >();
  const auto& rdof = g_inputdeck.get< tag::rdof >();
  std::size_t ncomp = U.nprop()/rdof;

  auto mindt = std::numeric_limits< tk::real >::max();

  for (std::size_t e=0; e<nelem; ++e)
  {
    // get conserved quantities at centroid
    std::vector< tk::real > B(rdof, 0.0);
    B[0] = 1.0;
    auto ugp = eval_state(ncomp, rdof, ndof, e, U, B);

    // Kinematic viscosity
    tk::real nu(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      auto alk = ugp[volfracIdx(nmat, k)];
      auto mu = getmatprop< tag::mu >(k);
      if (alk > 1.0e-04) nu = std::max(nu, alk*mu/ugp[densityIdx(nmat,k)]);
    }

    // characteristic length (radius of insphere)
    auto dx = std::min(std::cbrt(geoElem(e,0)), geoElem(e,4))
      /std::sqrt(24.0);

    // element dt
    mindt = std::min(mindt, dx*dx/nu);
  }

  return mindt;
}

void
resetSolidTensors(
  std::size_t nmat,
  std::size_t k,
  std::size_t e,
  tk::Fields& U,
  tk::Fields& P )
// *****************************************************************************
//  Reset the solid tensors
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] k Material id whose deformation gradient is required
//! \param[in] e Id of element whose solution is to be limited
//! \param[in/out] U High-order solution vector which gets modified
//! \param[in/out] P High-order vector of primitives which gets modified
// *****************************************************************************
{
  const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();
  const auto rdof = g_inputdeck.get< tag::rdof >();

  if (solidx[k] > 0) {
    for (std::size_t i=0; i<3; ++i) {
      for (std::size_t j=0; j<3; ++j) {
        // deformation gradient reset
        if (i==j) U(e, deformDofIdx(nmat, solidx[k], i, j, rdof, 0)) = 1.0;
        else U(e, deformDofIdx(nmat, solidx[k], i, j, rdof, 0)) = 0.0;

        // elastic Cauchy-stress reset
        P(e, stressDofIdx(nmat, solidx[k], stressCmp[i][j], rdof, 0)) = 0.0;

        // high-order reset
        for (std::size_t l=1; l<rdof; ++l) {
          U(e, deformDofIdx(nmat, solidx[k], i, j, rdof, l)) = 0.0;
          P(e, stressDofIdx(nmat, solidx[k], stressCmp[i][j], rdof, l)) = 0.0;
        }
      }
    }
  }
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
//! \return Inverse deformation gradient tensor (g_k) of material k
// *****************************************************************************
{
  const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();
  std::array< std::array< tk::real, 3 >, 3 > gk;

  if (solidx[k] > 0) {
    // deformation gradient for solids
    for (std::size_t i=0; i<3; ++i) {
      for (std::size_t j=0; j<3; ++j)
        gk[i][j] = state[deformIdx(nmat,solidx[k],i,j)];
    }
  }
  else {
    // empty vector for fluids
    gk = {{}};
  }

  return gk;
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
  const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();
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

std::size_t numSolids(
  std::size_t nmat,
  Kokkos::View<const size_t*, memory_space> solidx )
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
    if (solidx(k) > 0) ++nsld;

  return nsld;
}

} //inciter::
