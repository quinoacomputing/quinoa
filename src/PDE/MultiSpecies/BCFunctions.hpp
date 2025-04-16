// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/BCFunctions.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions specifying boundary conditions.
  \details   Functions that return boundary state when the interior state at
             at the boundary location is provided.
*/
// *****************************************************************************
#ifndef BCFunctions_h
#define BCFunctions_h

#include "FunctionPrototypes.hpp"
#include "MiscMultiSpeciesFns.hpp"

namespace inciter {

  //! \brief Boundary state function providing the left and right state of a
  //!   face at symmetry boundaries
  //! \param[in] ncomp Number of scalar components in this PDE system
  //! \param[in] ul Left (domain-internal) state
  //! \param[in] fn Unit face normal
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \note The function signature must follow tk::StateFn.
  static tk::StateFn::result_type
  symmetry( [[maybe_unused]] ncomp_t ncomp,
            const std::vector< EOS >&,
            const std::vector< tk::real >& ul,
            tk::real, tk::real, tk::real, tk::real,
            const std::array< tk::real, 3 >& fn )
  {
    auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

    Assert( ul.size() == ncomp, "Incorrect size for appended "
            "internal state vector" );

    tk::real rho(0.0);
    for (std::size_t k=0; k<nspec; ++k)
      rho += ul[multispecies::densityIdx(nspec, k)];

    auto ur = ul;

    // Internal cell velocity components
    auto v1l = ul[multispecies::momentumIdx(nspec, 0)]/rho;
    auto v2l = ul[multispecies::momentumIdx(nspec, 1)]/rho;
    auto v3l = ul[multispecies::momentumIdx(nspec, 2)]/rho;
    // Normal component of velocity
    auto vnl = v1l*fn[0] + v2l*fn[1] + v3l*fn[2];
    // Ghost state velocity components
    auto v1r = v1l - 2.0*vnl*fn[0];
    auto v2r = v2l - 2.0*vnl*fn[1];
    auto v3r = v3l - 2.0*vnl*fn[2];
    // Boundary condition
    ur[multispecies::momentumIdx(nspec, 0)] = rho * v1r;
    ur[multispecies::momentumIdx(nspec, 1)] = rho * v2r;
    ur[multispecies::momentumIdx(nspec, 2)] = rho * v3r;

    Assert( ur.size() == ncomp, "Incorrect size for appended "
            "boundary state vector" );

    return {{ std::move(ul), std::move(ur) }};
  }

  //! \brief Boundary state function providing the left and right state of a
  //!   face at farfield boundaries
  //! \param[in] ncomp Number of scalar components in this PDE system
  //! \param[in] ul Left (domain-internal) state
  //! \param[in] fn Unit face normal
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \details The farfield boudary calculation, implemented here, is
  //!   based on the characteristic theory of hyperbolic systems.
  //! \note The function signature must follow tk::StateFn
  static tk::StateFn::result_type
  farfield( [[maybe_unused]] ncomp_t ncomp,
            const std::vector< EOS >& mat_blk,
            const std::vector< tk::real >& ul,
            tk::real, tk::real, tk::real, tk::real,
            const std::array< tk::real, 3 >& fn )
  {
    auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

    // Farfield primitive quantities
    auto fp =
      g_inputdeck.get< tag::bc >()[0].get< tag::pressure >();
    auto ft =
      g_inputdeck.get< tag::bc >()[0].get< tag::temperature >();
    auto fu =
      g_inputdeck.get< tag::bc >()[0].get< tag::velocity >();
    auto fspec =
      g_inputdeck.get< tag::bc >()[0].get< tag::mass_fractions >();

    Assert( ul.size() == ncomp, "Incorrect size for appended "
            "internal state vector" );

    auto ur = ul;

    Mixture mixl(nspec); // Initialize mixture class
    mixl.set_state(ul, mat_blk);
    tk::real rhol = mixl.get_mix_density();

    // Internal cell velocity components
    auto v1l = ul[multispecies::momentumIdx(nspec, 0)]/rhol;
    auto v2l = ul[multispecies::momentumIdx(nspec, 1)]/rhol;
    auto v3l = ul[multispecies::momentumIdx(nspec, 2)]/rhol;

    // Normal velocity
    auto vn = v1l*fn[0] + v2l*fn[1] + v3l*fn[2];

    // Acoustic speed
    auto pl = mixl.pressure(rhol, v1l, v2l, v3l,
      ul[multispecies::energyIdx(nspec,0)], mat_blk);
    auto a = mixl.frozen_soundspeed(rhol, pl, mat_blk);

    // Mach number
    auto Ma = vn / a;

    if (Ma <= -1) {  // Supersonic inflow
      // For supersonic inflow, all the characteristics are from outside.
      // Therefore, we calculate the ghost cell state using the primitive
      // variables from outside.
      Mixture mixr(nspec);
      mixr.set_massfrac(fspec, fp, ft, mat_blk);

      tk::real rhor = mixr.get_mix_density();
      for (std::size_t k=0; k<nspec; ++k) {
        ur[multispecies::densityIdx(nspec,k)] = fspec[k] * rhor;
      }
      ur[multispecies::energyIdx(nspec,0)] = mixr.totalenergy(rhor, fu[0],
        fu[1], fu[2], fp, mat_blk);
      for (std::size_t i=0; i<3; ++i) {
        ur[multispecies::momentumIdx(nspec,i)] = rhor * fu[i];
      }

    } else if (Ma > -1 && Ma < 0) {  // Subsonic inflow
      // For subsonic inflow, there is 1 outgoing characteristic and 4
      // incoming characteristics. Therefore, we calculate the ghost cell state
      // by taking pressure from the internal cell and other quantities from
      // the outside.
      Mixture mixr(nspec);
      mixr.set_massfrac(fspec, pl, ft, mat_blk);

      tk::real rhor = mixr.get_mix_density();
      for (std::size_t k=0; k<nspec; ++k) {
        ur[multispecies::densityIdx(nspec,k)] = fspec[k] * rhor;
      }
      ur[multispecies::energyIdx(nspec,0)] = mixr.totalenergy(rhor, fu[0],
        fu[1], fu[2], pl, mat_blk);
      for (std::size_t i=0; i<3; ++i) {
        ur[multispecies::momentumIdx(nspec,i)] = rhor * fu[i];
      }

    } else if (Ma >= 0 && Ma < 1) {  // Subsonic outflow
      // For subsonic outflow, there is 1 incoming characteristic and 4
      // outgoing characteristics. Therefore, we calculate the ghost cell state
      // by taking pressure from the outside and other quantities from the
      // internal cell.
      std::vector< tk::real > massfrac_l;
      for (std::size_t k = 0; k < nspec; k++)
        massfrac_l[k] = ul[multispecies::densityIdx(nspec, k)] / rhol;
      tk::real tl = mixl.temperature(rhol, v1l, v2l, v3l,
        ul[multispecies::energyIdx(nspec,0)], mat_blk);
      Mixture mixr(nspec);
      mixr.set_massfrac(massfrac_l, fp, tl, mat_blk);

      ur[multispecies::energyIdx(nspec,0)] = mixr.totalenergy(rhol, v1l,
        v2l, v3l, fp, mat_blk);
    }
    // Otherwise, for supersonic outflow, all the characteristics are from
    // internal cell. Therefore, we calculate the ghost cell state using the
    // conservative variables from internal cell (which is what ur is
    // initialized to).

    Assert( ur.size() == ncomp, "Incorrect size for appended "
            "boundary state vector" );

    return {{ std::move(ul), std::move(ur) }};
  }

  //! \brief Boundary state function providing the left and right state of a
  //!   face at extrapolation boundaries
  //! \param[in] ul Left (domain-internal) state
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \note The function signature must follow tk::StateFn.
  static tk::StateFn::result_type
  extrapolate( ncomp_t,
               const std::vector< EOS >&,
               const std::vector< tk::real >& ul,
               tk::real, tk::real, tk::real, tk::real,
               const std::array< tk::real, 3 >& )
  {
    return {{ ul, ul }};
  }

  //! \brief Boundary state function providing the left and right state of a
  //!   face at no-slip wall boundaries
  //! \param[in] ncomp Number of scalar components in this PDE system
  //! \param[in] ul Left (domain-internal) state
//  //! \param[in] fn Unit face normal
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \note The function signature must follow tk::StateFn.
  static tk::StateFn::result_type
  noslipwall( [[maybe_unused]] ncomp_t ncomp,
              const std::vector< EOS >&,
              const std::vector< tk::real >& ul,
              tk::real, tk::real, tk::real, tk::real,
              const std::array< tk::real, 3 >& /*fn*/ )
  {
    auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

    Assert( ul.size() == ncomp, "Incorrect size for appended "
            "internal state vector" );

    tk::real rho(0.0);
    for (std::size_t k=0; k<nspec; ++k)
      rho += ul[multispecies::densityIdx(nspec, k)];

    auto ur = ul;

    // Internal cell velocity components
    auto v1l = ul[multispecies::momentumIdx(nspec, 0)]/rho;
    auto v2l = ul[multispecies::momentumIdx(nspec, 1)]/rho;
    auto v3l = ul[multispecies::momentumIdx(nspec, 2)]/rho;
    // Ghost state velocity components
    auto v1r = -v1l;
    auto v2r = -v2l;
    auto v3r = -v3l;
    // Boundary condition
    ur[multispecies::momentumIdx(nspec, 0)] = rho * v1r;
    ur[multispecies::momentumIdx(nspec, 1)] = rho * v2r;
    ur[multispecies::momentumIdx(nspec, 2)] = rho * v3r;

    Assert( ur.size() == ncomp, "Incorrect size for appended "
            "boundary state vector" );

    return {{ std::move(ul), std::move(ur) }};
  }

  //----------------------------------------------------------------------------
  // Boundary Gradient functions
  //----------------------------------------------------------------------------

  //! \brief Boundary gradient function copying the left gradient to the right
  //!   gradient at a face
  //! \param[in] dul Left (domain-internal) state
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \note The function signature must follow tk::StateFn.
  static tk::StateFn::result_type
  noOpGrad( ncomp_t,
            const std::vector< EOS >&,
            const std::vector< tk::real >& dul,
            tk::real, tk::real, tk::real, tk::real,
            const std::array< tk::real, 3 >& )
  {
    return {{ dul, dul }};
  }

  //! \brief Boundary gradient function for the symmetry boundary condition
  //! \param[in] ncomp Number of variables whos gradients are needed
  //! \param[in] dul Left (domain-internal) gradients
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \note The function signature must follow tk::StateFn.
  static tk::StateFn::result_type
  symmetryGrad( ncomp_t ncomp,
                const std::vector< EOS >&,
                const std::vector< tk::real >& dul,
                tk::real, tk::real, tk::real, tk::real,
                const std::array< tk::real, 3 >& )
  {
    Assert(dul.size() == 3*ncomp, "Incorrect size of boundary gradient vector");

    auto dur = dul;

    for (std::size_t i=0; i<3*ncomp; ++i)
      dur[i] = -dul[i];

    return {{ std::move(dul), std::move(dur) }};
  }

} // inciter::

#endif // BCFunctions_h
