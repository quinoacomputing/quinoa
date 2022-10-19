// *****************************************************************************
/*!
  \file      src/PDE/Riemann/Rusanov.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Rusanov Riemann flux function
  \details   This file implements the Rusanov Riemann solver, specific to ALECG.
*/
// *****************************************************************************
#ifndef Rusanov_h
#define Rusanov_h

#include <vector>

#include "Types.hpp"
#include "Fields.hpp"
#include "Tags.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EoS.hpp"

namespace inciter {

//! Rusanov approximate Riemann solver
struct Rusanov {

  using real = tk::real;

  //! Rusanov approximate Riemann solver flux function
  //! \param[in] mat_blk EOS material block
  //! \param[in] nx X component of the surface normal
  //! \param[in] ny Y component of the surface normal
  //! \param[in] nz Z component of the surface normal
  //! \param[in] mx X component of the weighted surface normal on chare
  //!   boundary, weighted by the number of contributions to the edge
  //! \param[in] my Y component of the weighted surface normal on chare
  //!   boundary, weighted by the number of contributions to the edge
  //! \param[in] mz Z component of the weighted surface normal on chare
  //!   boundary, weighted by the number of contributions to the edge
  //! \param[in] rL Left density
  //! \param[in] ruL Left X momentum
  //! \param[in] rvL Left Y momentum
  //! \param[in] rwL Left Z momentum
  //! \param[in] reL Left total specific energy
  //! \param[in] rR Right density
  //! \param[in] ruR Right X momentum
  //! \param[in] rvR Right Y momentum
  //! \param[in] rwR Right Z momentum
  //! \param[in] reR Right total specific energy
  //! \param[in] w1L Left X mesh velocity
  //! \param[in] w2L Left Y mesh velocity
  //! \param[in] w3L Left Z mesh velocity
  //! \param[in] w1R Right X mesh velocity
  //! \param[in] w2R Right Y mesh velocity
  //! \param[in] w3R Right Z mesh velocity
  //! \param[in] pL Left pressure
  //! \param[in] pR Right pressure
  //! \param[in,out] fr Riemann solution for density according to Rusanov
  //! \param[in,out] fru Riemann solution for X momenutm according to Rusanov
  //! \param[in,out] frv Riemann solution for Y momenutm according to Rusanov
  //! \param[in,out] frw Riemann solution for Z momenutm according to Rusanov
  //! \param[in,out] fre Riemann solution for specific total energy according
  //!   to Rusanov
  #pragma omp declare simd
  static void
  flux( const std::vector< EOS >& mat_blk,
        real nx, real ny, real nz,
        real mx, real my, real mz,
        real rL, real ruL, real rvL, real rwL, real reL,
        real rR, real ruR, real rvR, real rwR, real reR,
        real w1L, real w2L, real w3L, real w1R, real w2R, real w3R,
        real pL, real pR,
        real& fr, real& fru, real& frv, real& frw, real& fre )
  {
    auto ul = ruL/rL - w1L;
    auto vl = rvL/rL - w2L;
    auto wl = rwL/rL - w3L;

    auto ur = ruR/rR - w1R;
    auto vr = rvR/rR - w2R;
    auto wr = rwR/rR - w3R;

    auto al = mat_blk[0].eosCall< EOS::eos_soundspeed >( rL, pL );
    auto ar = mat_blk[0].eosCall< EOS::eos_soundspeed >( rR, pR );

    // dissipation
    real len = tk::length( {mx,my,mz} );
    real vml = ul*mx + vl*my + wl*mz;
    real vmr = ur*mx + vr*my + wr*mz;
    auto sl = std::abs(vml) + al*len;
    auto sr = std::abs(vmr) + ar*len;
    auto smax = std::max( sl, sr );

    // face-normal velocities
    real vnl = ul*nx + vl*ny + wl*nz;
    real vnr = ur*nx + vr*ny + wr*nz;

    // numerical fluxes
    fr  = 0.5*(rL*vnl + rR*vnr - smax*(rR - rL));
    fru = 0.5*(ruL*vnl + pL*nx + ruR*vnr + pR*nx - smax*(ruR - ruL));
    frv = 0.5*(rvL*vnl + pL*ny + rvR*vnr + pR*ny - smax*(rvR - rvL));
    frw = 0.5*(rwL*vnl + pL*nz + rwR*vnr + pR*nz - smax*(rwR - rwL));
    fre = 0.5*(reL*vnl + reR*vnr
               + pL*(ruL*nx + rvL*ny + rwL*nz)/rL
               + pR*(ruR*nx + rvR*ny + rwR*nz)/rR
               - smax*(reR - reL));
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::Rusanov; }
};

} // inciter::

#endif // Rusanov_h
