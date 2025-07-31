// *****************************************************************************
/*!
  \file      src/PDE/Riemann/AUSMMultiSpecies.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Advection Upstream Splitting Method (AUSM+) Riemann flux function
             for multi-species fluid dynamics.
  \details   This file implements the Advection Upstream Splitting Method (AUSM)
             Riemann solver, with the all-speed corrections.
             Ref. Liou, M. S. (2006). A sequel to AUSM, Part II: AUSM+-up for
             all speeds. Journal of computational physics, 214(1), 137-170.
*/
// *****************************************************************************
#ifndef AUSMMultiSpecies_h
#define AUSMMultiSpecies_h

#include <vector>

#include "Fields.hpp"
#include "SplitMachFns.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EOS.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"
#include "MultiSpecies/Mixture/Mixture.hpp"

namespace inciter {

//! AUSMMultiSpecies+up approximate Riemann solver
struct AUSMMultiSpecies {

  //! AUSM+up approximate Riemann solver flux function for multi-species flow
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann flux solution according to AUSM+up, appended by Riemann
  //!   velocities and volume-fractions.
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::vector< EOS >& mat_blk,
        const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& = {} )
  {
    auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

    // All-speed parameters
    // These parameters control the amount of all-speed diffusion necessary for
    // low-Mach flows. Setting k_u and k_p to zero does not add any all-speed
    // diffusion, whereas setting k_u and k_p to 1 adds maximum recommended
    // all-speed diffusion. See "Liou, M. S. (2006). A sequel to AUSM, Part II:
    // AUSM+-up for all speeds. Journal of computational physics, 214(1),
    // 137-170" for more mathematical explanation. k_u is the velocity diffusion
    // term and k_p is the pressure diffusion term. These two terms reduce
    // pressure-velocity decoupling (chequerboarding/odd-even oscillations).
    auto k_u = g_inputdeck.get< tag::lowspeed_ku >();
    auto k_p = g_inputdeck.get< tag::lowspeed_kp >();

    auto ncomp = u[0].size()-1;
    std::vector< tk::real > flx( ncomp, 0 );

    tk::real rhol(0.0), rhor(0.0), pl(0.0), pr(0.0), hl(0.0), hr(0.0),
      Tl(0.0), Tr(0.0), al(0.0), ar(0.0), a12(0.0), rho12(0.0);

    // initialize mixtures
    Mixture mixl(nspec, u[0], mat_blk);
    Mixture mixr(nspec, u[1], mat_blk);

    // Mixture densities
    rhol = mixl.get_mix_density();
    rhor = mixr.get_mix_density();
    Tl = u[0][ncomp+multispecies::temperatureIdx(nspec, 0)];
    Tr = u[1][ncomp+multispecies::temperatureIdx(nspec, 0)];

    // Velocities
    auto ul = u[0][multispecies::momentumIdx(nspec, 0)]/rhol;
    auto vl = u[0][multispecies::momentumIdx(nspec, 1)]/rhol;
    auto wl = u[0][multispecies::momentumIdx(nspec, 2)]/rhol;
    auto ur = u[1][multispecies::momentumIdx(nspec, 0)]/rhor;
    auto vr = u[1][multispecies::momentumIdx(nspec, 1)]/rhor;
    auto wr = u[1][multispecies::momentumIdx(nspec, 2)]/rhor;

    pl = mixl.pressure( rhol, Tl );
    hl = u[0][multispecies::energyIdx(nspec, 0)] + pl;
    al = mixl.frozen_soundspeed( rhol, Tl, mat_blk );

    pr = mixr.pressure( rhor, Tr );
    hr = u[1][multispecies::energyIdx(nspec, 0)] + pr;
    ar = mixr.frozen_soundspeed( rhor, Tr, mat_blk );

    // Average states for mixture speed of sound
    a12 = 0.5*(al+ar);
    rho12 = 0.5*(rhol+rhor);

    // Face-normal velocities from advective velocities
    auto vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];
    auto vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];

    // Mach numbers
    auto ml = vnl/a12;
    auto mr = vnr/a12;

    tk::real f_a(1.0);

    // Split Mach polynomials
    auto msl = splitmach_ausm( ml, f_a );
    auto msr = splitmach_ausm( mr, f_a );

    // Riemann Mach number
    auto m0 = 1.0 - (0.5*(vnl*vnl + vnr*vnr)/(a12*a12));
    auto mp = -k_p* std::max(m0, 0.0) * (pr-pl) / (f_a*rho12*a12*a12);
    auto m12 = msl[0] + msr[1] + mp;
    auto vriem = a12 * m12;

    // Riemann pressure
    auto pu = -k_u* msl[2] * msr[3] * f_a * rho12 * a12 * (vnr-vnl);
    auto p12 = msl[2]*pl + msr[3]*pr + pu;

    // Flux vector splitting
    auto l_plus = 0.5 * (vriem + std::fabs(vriem));
    auto l_minus = 0.5 * (vriem - std::fabs(vriem));

    // Conservative fluxes
    for (std::size_t k=0; k<nspec; ++k)
    {
      flx[multispecies::densityIdx(nspec, k)] =
        l_plus *u[0][multispecies::densityIdx(nspec, k)] +
        l_minus*u[1][multispecies::densityIdx(nspec, k)];
    }

    flx[multispecies::energyIdx(nspec, 0)] = l_plus*hl + l_minus*hr;

    for (std::size_t idir=0; idir<3; ++idir)
    {
      flx[multispecies::momentumIdx(nspec, idir)] =
        l_plus *u[0][multispecies::momentumIdx(nspec, idir)] +
        l_minus*u[1][multispecies::momentumIdx(nspec, idir)] + p12*fn[idir];
    }

    return flx;
  }

  //! AUSM+up approximate Riemann solver flux Jacbobian function for
  //! multi-species flow
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann flux Jacobian according to AUSM+up.
  //! \note The function signature must follow tk::RiemannFluxJacFn
  static tk::RiemannFluxJacFn::result_type
  fluxJac( const std::vector< EOS >& mat_blk,
           const std::array< tk::real, 3 >& fn,
           const std::array< std::vector< tk::real >, 2 >& u,
           const std::vector< std::array< tk::real, 3 > >& = {} )
  {
    auto k_u = g_inputdeck.get< tag::lowspeed_ku >();
    auto k_p = g_inputdeck.get< tag::lowspeed_kp >();
    auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();
    auto ncomp = u[0].size()-1;
    tk::real rhol(0.0), rhor(0.0), pl(0.0), pr(0.0), Tl(0.0), Tr(0.0), al(0.0),
             ar(0.0), a12(0.0), rho12(0.0), f_a(1.0);
    std::size_t uid, vid, wid, Tid;

    // Variables with partials of the left and right states that can be
    // overwritten each loop iteration
    std::array< tk::real, 2 > dmldP, dmrdP, dmbardP, dmpdP, da12dP,
                              dm12dP, dpudP, dvriemdP;

    // Factory for variables that have partials with respect to both the left
    // and right states and need to be saved
    auto drl_fact = [ncomp] {
      return std::array{ std::vector< tk::real >(ncomp),
                         std::vector< tk::real >(ncomp) };
    };
    auto dl_plusdP = drl_fact(), dl_minusdP = drl_fact(), dp12dP = drl_fact(),
         drho12dP = drl_fact();

    // Flux derivatives with respect to primitive and conserved variables
    std::array dFdP{ std::vector(ncomp, std::vector< tk::real >(ncomp)),
                     std::vector(ncomp, std::vector< tk::real >(ncomp)) },
               dFdU{ std::vector(ncomp, std::vector< tk::real >(ncomp)),
                     std::vector(ncomp, std::vector< tk::real >(ncomp)) };

    // initialize mixtures
    Mixture mixl(nspec, u[0], mat_blk);
    Mixture mixr(nspec, u[1], mat_blk);

    // Mixture densities
    rhol = mixl.get_mix_density();
    rhor = mixr.get_mix_density();

    // Create vectors of conserved and primitive variables for ease of access
    auto U = drl_fact(), P = drl_fact();
    uid = nspec; vid = nspec + 1; wid = nspec + 2; Tid = ncomp - 1;
    for (std::size_t k=0; k<nspec; ++k)
    {
      U[0][k] = u[0][multispecies::densityIdx(nspec, k)];
      U[1][k] = u[1][multispecies::densityIdx(nspec, k)];
      P[0][k] = u[0][multispecies::densityIdx(nspec, k)];
      P[1][k] = u[1][multispecies::densityIdx(nspec, k)];
    }

    for (std::size_t idir=0; idir<3; ++idir)
    {
      U[0][idir + uid] = u[0][multispecies::momentumIdx(nspec, idir)];
      U[1][idir + uid] = u[1][multispecies::momentumIdx(nspec, idir)];
      P[0][idir + uid] = u[0][multispecies::momentumIdx(nspec, idir)] / rhol;
      P[1][idir + uid] = u[1][multispecies::momentumIdx(nspec, idir)] / rhor;
    }

    U[0][Tid] = u[0][multispecies::energyIdx(nspec, 0)];
    U[1][Tid] = u[1][multispecies::energyIdx(nspec, 0)];
    P[0][Tid] = u[0][ncomp+multispecies::temperatureIdx(nspec, 0)];
    P[1][Tid] = u[1][ncomp+multispecies::temperatureIdx(nspec, 0)];

    // Quantities to be calculated before looping over primitives
    for (std::size_t k=0; k<nspec; ++k)
    {
      drho12dP[0][k] = 0.5;
      drho12dP[1][k] = 0.5;
    }

    pl = mixl.pressure( rhol, Tl );
    al = mixl.frozen_soundspeed( rhol, Tl, mat_blk );

    pr = mixr.pressure( rhor, Tr );
    ar = mixr.frozen_soundspeed( rhor, Tr, mat_blk );

    // Average states for mixture speed of sound
    a12 = 0.5*(al+ar);
    rho12 = 0.5*(rhol+rhor);

    auto dUdP = conservedPrimitiveJac(mat_blk, fn, u);

    // Face-normal velocities from advective velocities
    auto vnl = P[0][uid] * fn[0]
             + P[0][vid] * fn[1]
             + P[0][wid] * fn[2];
    auto vnr = P[1][uid] * fn[0]
             + P[1][vid] * fn[1]
             + P[1][wid] * fn[2];

    // Mach numbers
    auto ml = vnl/a12;
    auto mr = vnr/a12;

    // Flux split vectors and their derivatives w.r.t. mach
    auto msl = splitmach_ausm( ml, f_a );
    auto msr = splitmach_ausm( mr, f_a );
    auto dmsldm = splitmach_derivs(ml, f_a);
    auto dmsrdm = splitmach_derivs(mr, f_a);

    // Pressure, sound speed, and normal velocity derivatives
    auto daldP = mixl.soundspeed_prim_partials(rhol, P[0][Tid], mat_blk);
    auto dardP = mixr.soundspeed_prim_partials(rhor, P[1][Tid], mat_blk);
    auto dpldP = mixl.pressure_prim_partials(rhol, P[0][Tid], mat_blk);
    auto dprdP = mixr.pressure_prim_partials(rhor, P[1][Tid], mat_blk);
    std::vector< tk::real > dvnldP(ncomp, 0.0), dvnrdP(ncomp, 0.0);
    std::transform(fn.begin(), fn.end(), dvnldP.begin() + uid,
                   [](tk::real ni){ return ni; });
    std::transform(fn.begin(), fn.end(), dvnrdP.begin() + uid,
                   [](tk::real ni){ return ni; });

    // M_1/2 derivative pre-calculations
    auto mbar = std::sqrt(0.5*(vnl*vnl + vnr*vnr)/(a12*a12));
    auto idenom = 1. / (f_a * rho12 * a12 * a12);
    auto max_mbar = std::max(1. - mbar*mbar, 0.0);
    auto mp = -k_p* max_mbar * (pr - pl) * idenom;
    auto C = mbar*mbar > 1 ? 0.0 : 2. * k_p * mbar * (pr - pl) * idenom;
    auto m12 = msl[0] + msr[1] + mp;
    auto vriem = a12 * m12;
    auto l_plus = 0.5 * (vriem + std::fabs(vriem));
    auto l_minus = 0.5 * (vriem - std::fabs(vriem));

    // P_1/2 derivative pre-caluclations
    auto pu = -k_u* msl[2] * msr[3] * f_a * rho12 * a12 * (vnr-vnl);

    // Loop over primitives that fluxes are taken derivatives with
    for (std::size_t k=0; k<ncomp; ++k)
    {
      // Mach number derivatives and face SoS and density derivatives
      dmldP[0] = -ml / (2. * a12) * daldP[k];
      dmldP[1] = -ml / (2. * a12) * dardP[k];
      dmrdP[0] = -mr / (2. * a12) * daldP[k];
      dmrdP[1] = -mr / (2. * a12) * dardP[k];
      da12dP[0] = 0.5 * daldP[k];
      da12dP[1] = 0.5 * dardP[k];

      // Mach derivatives pick up an extra term when taken with velocity on the
      // same side
      if (k >= uid && k <= wid)
      {
        dmldP[0] += fn[k - uid]/a12;
        dmrdP[1] += fn[k - uid]/a12;
      }

      // M^bar has a singularity when equal to zero, handle explicitly
      if (std::fabs(mbar) < 1e-12){
        dmbardP[0] = 0.5 * (dmldP[0] + dmrdP[0]);
        dmbardP[1] = 0.5 * (dmldP[1] + dmrdP[1]);
      } else {
        dmbardP[0] = 0.5 /  mbar * (ml * dmldP[0] + mr * dmrdP[0]);
        dmbardP[1] = 0.5 /  mbar * (ml * dmldP[1] + mr * dmrdP[1]);
      }

      dmpdP[0] = C * dmbardP[0]
               + k_p * max_mbar * idenom * dpldP[k]
               - mp / rho12 *drho12dP[0][k]
               - 2. * mp / a12 * da12dP[0];
      dmpdP[1] = C * dmbardP[1]
               - k_p * max_mbar * idenom * dprdP[k]
               - mp / rho12 *drho12dP[1][k]
               - 2. * mp / a12 * da12dP[1];

      dm12dP[0] = dmsldm[0] * dmldP[0]
                + dmsrdm[1] * dmrdP[0]
                + dmpdP[0];
      dm12dP[1] = dmsldm[0] * dmldP[1]
                + dmsrdm[1] * dmrdP[1]
                + dmpdP[1];

      // P_1/2 derivatives
      dpudP[0] = -k_u * msr[3] * f_a * rho12 * a12 * (vnr-vnl) * dmsldm[2]
                 * dmldP[0]
               + -k_u * msl[2] * f_a * rho12 * a12 * (vnr-vnl) * dmsrdm[3]
                 * dmrdP[0]
               + pu / rho12 * drho12dP[0][k]
               + pu / a12 * da12dP[0]
               + k_u* msl[2] * msr[3] * f_a * rho12 * a12 * dvnldP[k];
      dpudP[1] = -k_u * msr[3] * f_a * rho12 * a12 * (vnr-vnl) * dmsldm[2]
                 * dmldP[1]
               + -k_u * msl[2] * f_a * rho12 * a12 * (vnr-vnl) * dmsrdm[3]
                 * dmrdP[1]
               + pu / rho12 * drho12dP[1][k]
               + pu / a12 * da12dP[1]
               - k_u* msl[2] * msr[3] * f_a * rho12 * a12 * dvnrdP[k];

      dp12dP[0][k] = dmsldm[2] * dmldP[0] * pl
                   + msl[2] * dpldP[k]
                   + dmsrdm[3] * dmrdP[0] * pr
                   + dpudP[0];
      dp12dP[1][k] = dmsldm[2] * dmldP[1] * pl
                   + dmsrdm[3] * dmrdP[1] * pr
                   + msr[3] * dprdP[k]
                   + dpudP[1];

      dvriemdP[0] = da12dP[0] * m12 + a12 * dm12dP[0];
      dvriemdP[1] = da12dP[1] * m12 + a12 * dm12dP[1];

      if (vriem > 0) {
        dl_plusdP[0][k]  = dvriemdP[0];
        dl_plusdP[1][k]  = dvriemdP[1];
        dl_minusdP[0][k] = 0.;
        dl_minusdP[1][k] = 0.;
      } else {
        dl_plusdP[0][k]  = 0.;
        dl_plusdP[1][k]  = 0.;
        dl_minusdP[0][k] = dvriemdP[0];
        dl_minusdP[1][k] = dvriemdP[1];
      }
    }

    // Loop over flux derivatives are taken on
    for (std::size_t k=0; k<ncomp; ++k)
    {
      // Loop over primitives derivatives are taken with respect to
      for (std::size_t l=0; l<ncomp; ++l)
      {
        dFdP[0][k][l] = U[0][k] * dl_plusdP[0][l]
                      + U[1][k] * dl_minusdP[0][l]
                      + dUdP[0][k][l] * l_plus;

        dFdP[1][k][l] = U[0][k] * dl_plusdP[1][l]
                      + U[1][k] * dl_minusdP[1][l]
                      + dUdP[1][k][l] * l_minus;

        // Pressure split terms: only added to momentum fluxes
        if ( k >= uid && k <= wid )
        {
          dFdP[0][k][l] += dp12dP[0][l] * fn[k - uid];
          dFdP[1][k][l] += dp12dP[1][l] * fn[k - uid];
        }

        // Add pressure to energy to recover enthalpy
        if ( k == Tid )
        {
          dFdP[0][k][l] += pl * dl_plusdP[0][l]
                         + dpldP[l] * l_plus
                         + pr * dl_minusdP[0][l];

          dFdP[1][k][l] += pl * dl_plusdP[1][l]
                         + dprdP[l] * l_minus
                         + pr * dl_minusdP[1][l];
        }
      }
    }

    // Compute dUdP^-1 to get dPdU
    double dPldU[ncomp*ncomp], dPrdU[ncomp*ncomp];
    for (std::size_t i=0; i<ncomp; ++i)
      for (std::size_t j=0; j<ncomp; ++j)
      {
        dPldU[ncomp*i+j] = dUdP[0][i][j];
        dPrdU[ncomp*i+j] = dUdP[1][i][j];
      }

    lapack_int ipiv[ncomp];
    #ifndef NDEBUG
    lapack_int ierr =
    #endif
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, ncomp, ncomp, dPldU, ncomp, ipiv);
    Assert(ierr==0, "Lapack error in LU factorization of dUdPl");
    #ifndef NDEBUG
    lapack_int jerr =
    #endif
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, ncomp, dPldU, ncomp, ipiv);
    Assert(jerr==0, "Lapack error in inverting dUdPl");
    #ifndef NDEBUG
    lapack_int ierr =
    #endif
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, ncomp, ncomp, dPrdU, ncomp, ipiv);
    Assert(ierr==0, "Lapack error in LU factorization of dUdPr");
    #ifndef NDEBUG
    lapack_int jerr =
    #endif
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, ncomp, dPrdU, ncomp, ipiv);
    Assert(jerr==0, "Lapack error in inverting dUdPr");

    // Final matrix multiplication
    // We have dF_i/dP_j, (flux in first index, primitive derivative in the
    // second) and dP_i/dU_j (primitive in first index, conserved derivative in
    // the second). We want dF_i/dP_k dP_k/dU_j = dF_i/dU_j
    for (std::size_t i=0; i<ncomp; ++i)
      for (std::size_t j=0; j<ncomp; ++j)
        for (std::size_t k=0; k<ncomp; ++k)
        {
          dFdU[0][i][j] += dFdP[0][i][k] * dPldU[ncomp*j+k];
          dFdU[1][i][j] += dFdP[1][i][k] * dPrdU[ncomp*j+k];
        }

    return dFdU;
  }

  //! NOTE: this is a fairly general function, and likely shouldn't be here,
  //! but I'm not sure the best place for it.
  //! Calculates the Jacobian of the converved variables with respect to the
  //! primitive variables
  //! \param[in] u Left and right unknown/state vector
  //! \return Derivatives of conserved variables with respect to primitive
  //! variables for the left and right states
  //! \note The function signature must follow tk::RiemannFluxJacFn
  static tk::RiemannFluxJacFn::result_type
  conservedPrimitiveJac(
    const std::vector< EOS >& mat_blk,
    const std::array< tk::real, 3 >& ,
    const std::array< std::vector< tk::real >, 2 >& u,
    const std::vector< std::array< tk::real, 3 > >& = {} )
  {
    auto ncomp = u[0].size()-1;
    auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();
    std::array dUdP{ std::vector(ncomp, std::vector< tk::real >(ncomp)),
                     std::vector(ncomp, std::vector< tk::real >(ncomp)) };
    tk::real rhol(0.0), rhor(0.0), Tl(0.0), Tr(0.0);

    // Initialize mixtures
    Mixture mixl(nspec, u[0], mat_blk);
    Mixture mixr(nspec, u[1], mat_blk);

    Tl = u[0][ncomp+multispecies::temperatureIdx(nspec, 0)];
    Tr = u[1][ncomp+multispecies::temperatureIdx(nspec, 0)];
    rhol = mixl.get_mix_density();
    rhor = mixr.get_mix_density();

    // Velocities
    auto ul = u[0][multispecies::momentumIdx(nspec, 0)]/rhol;
    auto vl = u[0][multispecies::momentumIdx(nspec, 1)]/rhol;
    auto wl = u[0][multispecies::momentumIdx(nspec, 2)]/rhol;
    auto ur = u[1][multispecies::momentumIdx(nspec, 0)]/rhor;
    auto vr = u[1][multispecies::momentumIdx(nspec, 1)]/rhor;
    auto wr = u[1][multispecies::momentumIdx(nspec, 2)]/rhor;

    // Partials of species density w.r.t. species density
    for (std::size_t k=0; k<nspec; ++k)
    {
      dUdP[0][multispecies::densityIdx(nspec, k)]
          [multispecies::densityIdx(nspec, k)] = 1.0;
      dUdP[1][multispecies::densityIdx(nspec, k)]
          [multispecies::densityIdx(nspec, k)] = 1.0;
    }

    // Partials of momentum...
    // ... w.r.t. species density
    for (std::size_t k=0; k<nspec; ++k)
    {
      dUdP[0][multispecies::momentumIdx(nspec, 0)]
          [multispecies::densityIdx(nspec, k)] = ul;
      dUdP[0][multispecies::momentumIdx(nspec, 1)]
          [multispecies::densityIdx(nspec, k)] = vl;
      dUdP[0][multispecies::momentumIdx(nspec, 2)]
          [multispecies::densityIdx(nspec, k)] = wl;
      dUdP[1][multispecies::momentumIdx(nspec, 0)]
          [multispecies::densityIdx(nspec, k)] = ur;
      dUdP[1][multispecies::momentumIdx(nspec, 1)]
          [multispecies::densityIdx(nspec, k)] = vr;
      dUdP[1][multispecies::momentumIdx(nspec, 2)]
          [multispecies::densityIdx(nspec, k)] = wr;
    }

    // ... w.r.t. velocity
    for (std::size_t idir=0; idir<3; ++idir)
    {
      dUdP[0][multispecies::momentumIdx(nspec, idir)]
          [multispecies::momentumIdx(nspec, idir)] = rhol;
      dUdP[1][multispecies::momentumIdx(nspec, idir)]
          [multispecies::momentumIdx(nspec, idir)] = rhor;
    }

    // Partials of energy...
    // ... w.r.t. species density
    for (std::size_t k=0; k<nspec; ++k)
    {
      // This assumes de_s/drho_s = 0
      dUdP[0][multispecies::energyIdx(nspec, 0)]
          [multispecies::densityIdx(nspec, k)] =
          mat_blk[k].compute< EOS::internalenergy >(Tl) +
          0.5 * (ul*ul + vl*vl + wl*wl);
      dUdP[1][multispecies::energyIdx(nspec, 0)]
          [multispecies::densityIdx(nspec, k)] =
          mat_blk[k].compute< EOS::internalenergy >(Tr) +
          0.5 * (ur*ur + vr*vr + wr*wr);
    }
    // ... w.r.t. velocity
    dUdP[0][multispecies::energyIdx(nspec, 0)]
        [multispecies::momentumIdx(nspec, 0)] = rhol * ul;
    dUdP[0][multispecies::energyIdx(nspec, 0)]
        [multispecies::momentumIdx(nspec, 1)] = rhol * vl;
    dUdP[0][multispecies::energyIdx(nspec, 0)]
        [multispecies::momentumIdx(nspec, 2)] = rhol * wl;
    dUdP[1][multispecies::energyIdx(nspec, 0)]
        [multispecies::momentumIdx(nspec, 0)] = rhor * ur;
    dUdP[1][multispecies::energyIdx(nspec, 0)]
        [multispecies::momentumIdx(nspec, 1)] = rhor * vr;
    dUdP[1][multispecies::energyIdx(nspec, 0)]
        [multispecies::momentumIdx(nspec, 2)] = rhor * wr;

    // ... w.r.t. temperature
    dUdP[0][multispecies::energyIdx(nspec, 0)]
        [multispecies::energyIdx(nspec, 0)] = rhol * mixl.mix_Cv(Tl, mat_blk);
    dUdP[1][multispecies::energyIdx(nspec, 0)]
        [multispecies::energyIdx(nspec, 0)] = rhor * mixl.mix_Cv(Tr, mat_blk);

    return dUdP;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::AUSM; }
};

} // inciter::

#endif // AUSMMultiSpecies_h
