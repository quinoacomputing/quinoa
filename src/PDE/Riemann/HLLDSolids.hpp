// *****************************************************************************
/*!
  \file      src/PDE/Riemann/HLLDSolids.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2023 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     HLLD Riemann flux function for solids
  \details   This file implements the HLLD Riemann solver for solids.
             Ref. Barton, P. T. (2019). An interface-capturing Godunov method
             for the simulation of compressible solid-fluid problems. Journal of
             Computational Physics, 390, 25-50.
*/
// *****************************************************************************
#ifndef HLLDSolids_h
#define HLLDSolids_h

#include <vector>

#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EOS.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"

namespace inciter {

//! HLLD approximate Riemann solver for solids
struct HLLDSolids {

  //! HLLD approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann flux solution according to HLLD, appended by
  //!   Riemann velocities and volume-fractions.
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::vector< EOS >& mat_blk,
        const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& )
  {
    auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
    const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

    auto nsld = numSolids(nmat, solidx);
    auto ncomp = u[0].size()-(3+nmat+nsld*6);
    std::vector< tk::real > flx(ncomp,0), fl(ncomp,0), fr(ncomp, 0);
    // Primitive quantities
    tk::real arhol(0.0), arhor(0.0);
    tk::real amatl(0.0), amatr(0.0), asmatl(0.0), asmatr(0.0),
      ac_l(0.0), ac_r(0.0), asc_l(0.0), asc_r(0.0);
    std::vector< tk::real > al_l(nmat, 0.0), al_r(nmat, 0.0),
                            pml(nmat, 0.0), pmr(nmat, 0.0);
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > g_l, g_r,
      asig_l, asig_r;
    std::vector< std::array< tk::real, 3 > > asign_l, asign_r;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gn_l, gn_r,
      asignn_l, asignn_r;
    std::array< tk::real, 3 > sign_l {{0, 0, 0}}, sign_r {{0, 0, 0}};

    // Independently limited velocities for advection
    auto ul = u[0][ncomp+velocityIdx(nmat, 0)];
    auto vl = u[0][ncomp+velocityIdx(nmat, 1)];
    auto wl = u[0][ncomp+velocityIdx(nmat, 2)];
    auto ur = u[1][ncomp+velocityIdx(nmat, 0)];
    auto vr = u[1][ncomp+velocityIdx(nmat, 1)];
    auto wr = u[1][ncomp+velocityIdx(nmat, 2)];

    // Rotated velocities from advective velocities
    auto vnl = tk::rotateVector({ul, vl, wl}, fn);
    auto vnr = tk::rotateVector({ur, vr, wr}, fn);

    for (std::size_t k=0; k<nmat; ++k)
    {
      // Left outer state
      // -----------------------------------------------------------------------
      al_l[k] = u[0][volfracIdx(nmat, k)];
      pml[k] = u[0][ncomp+pressureIdx(nmat, k)];
      arhol = u[0][densityIdx(nmat, k)];
      brhol += arhol;

      // inv deformation gradient and Cauchy stress tensors
      g_l.push_back(getDeformGrad(nmat, k, u[0]));
      asig_l.push_back(getCauchyStress(nmat, k, ncomp, u[0]));
      for (std::size_t i=0; i<3; ++i) asig_l[k][i][i] -= pml[k];

      // normal stress (traction) vector
      asign_l.push_back(tk::matvec(asig_l[k], fn));
      for (std::size_t i=0; i<3; ++i)
        sign_l[i] += asign_l[k][i];

      // rotate stress vector
      asignn_r.push_back(tk::rotateTensor(asig_r[k], fn));

      // rotate deformation gradient tensor
      gn_l.push_back(tk::rotateTensor(g_l[k], fn));
      amatl = mat_blk[k].compute< EOS::soundspeed >(
        u[0][densityIdx(nmat, k)], pml[k], al_l[k], k, gn_l[k] );

      // Right outer state
      // -----------------------------------------------------------------------
      al_r[k] = u[1][volfracIdx(nmat, k)];
      pmr[k] = u[1][ncomp+pressureIdx(nmat, k)];
      arhor = u[1][densityIdx(nmat, k)];
      brhor += arhor;

      // inv deformation gradient and Cauchy stress tensors
      g_r.push_back(getDeformGrad(nmat, k, u[1]));
      asig_r.push_back(getCauchyStress(nmat, k, ncomp, u[1]));
      for (std::size_t i=0; i<3; ++i) asig_r[k][i][i] -= pmr[k];

      // normal stress (traction) vector
      asign_r.push_back(tk::matvec(asig_r[k], fn));
      for (std::size_t i=0; i<3; ++i)
        sign_r[i] += asign_r[k][i];

      // rotate stress vector
      asignn_r.push_back(tk::rotateTensor(asig_r[k], fn));
      
      // rotate deformation gradient tensor
      gn_r.push_back(tk::rotateTensor(g_r[k], fn));
      amatr = mat_blk[k].compute< EOS::soundspeed >(
        u[1][densityIdx(nmat, k)], pmr[k], al_r[k], k, gn_r[k] );

      // Signal velocities
      auto Sl_minus = std::min((vnl[0]-amatl), (vnr-amatr));
      auto Sl_plus = std::max((vnl[0]+amatl), (vnr+amatr));
      auto Si = (rhol*vnl[0]*(Sl_minus-vnl[0]) - rhor*vnr[0]*(Sl_plus-vnr[0])
        + signn_l[0][0]-signn_r[0][0])
        / (rhol*(Sl_minus-vnl[0]) - rhor*(Sl_plus-vnr[0]));
      
      // Left tilde state
      // -----------------------------------------------------------------------
      al_t_l[k] = al_l[k]*(Sl_minus-vnl[0])/(Sl_minus-Si);
      pm_t_l[k] = pml[k]*(Sl_minus-vnl[0])/(Sl_minus-Si);
      arho_t_l = arhol*(Sl_minus-vnl[0])/(Sl_minus-Si);
      brho_r_l += arho_t_l;

      // inv deformation gradient and Cauchy stress tensors
      gn_t_l.push_back(getDeformGrad(gn_l));
      asig_l.push_back(getCauchyStress(nmat, k, ncomp, u[0]));
      for (std::size_t i=0; i<3; ++i) asig_l[k][i][i] -= pml[k];

      // normal stress (traction) vector
      asign_l.push_back(tk::matvec(asig_l[k], fn));
      for (std::size_t i=0; i<3; ++i)
        sign_l[i] += asign_l[k][i];

      // rotate stress vector
      asignn_r.push_back(tk::rotateTensor(asig_r[k], fn));

      // Compute shear speed of sounds with tilde variables
      amatl = mat_blk[k].compute< EOS::soundspeed >(
        arho_t_l, pm_t_l[k], al_t_l[k], k, gn_t_l[k] );

      // Right tilde state
      // -----------------------------------------------------------------------
      al_r[k] = u[1][volfracIdx(nmat, k)];
      pmr[k] = u[1][ncomp+pressureIdx(nmat, k)];
      arhor = u[1][densityIdx(nmat, k)];
      brhor += arhor;

      // inv deformation gradient and Cauchy stress tensors
      g_r.push_back(getDeformGrad(nmat, k, u[1]));
      asig_r.push_back(getCauchyStress(nmat, k, ncomp, u[1]));
      for (std::size_t i=0; i<3; ++i) asig_r[k][i][i] -= pmr[k];

      // normal stress (traction) vector
      asign_r.push_back(tk::matvec(asig_r[k], fn));
      for (std::size_t i=0; i<3; ++i)
        sign_r[i] += asign_r[k][i];

      // rotate stress vector
      asignn_r.push_back(tk::rotateTensor(asig_r[k], fn));
      
      // rotate deformation gradient tensor
      gn_r.push_back(tk::rotateTensor(g_r[k], fn));
      amatr = mat_blk[k].compute< EOS::soundspeed >(
        u[1][densityIdx(nmat, k)], pmr[k], al_r[k], k, gn_r[k] );
      
    }


    
    Assert( flx.size() == (ncomp+nmat+1+3*nsld), "Size of "
            "multi-material flux vector incorrect" );

    return flx;
  }

  ////! Flux type accessor
  ////! \return Flux type
  //static ctr::FluxType type() noexcept {
  //  return ctr::FluxType::HLLDSolids; }
};

} // inciter::

#endif // HLLDSolids_h
