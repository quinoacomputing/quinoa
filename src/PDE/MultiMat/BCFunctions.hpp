// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/BCFunctions.hpp
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

namespace inciter {

  //! \brief Boundary state function providing the left and right state of a
  //!   face at symmetry boundaries
  //! \param[in] system Equation system index
  //! \param[in] ncomp Number of scalar components in this PDE system
  //! \param[in] ul Left (domain-internal) state
  //! \param[in] fn Unit face normal
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \note The function signature must follow tk::StateFn. For multimat, the
  //!   left or right state is the vector of conserved quantities, followed by
  //!   the vector of primitive quantities appended to it.
  static tk::StateFn::result_type
  symmetry( ncomp_t system, ncomp_t ncomp,
            const std::vector< EOS >&,
            const std::vector< tk::real >& ul,
            tk::real, tk::real, tk::real, tk::real,
            const std::array< tk::real, 3 >& fn )
  {
    const auto nmat =
      g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[system];

    Assert( ul.size() == ncomp+nmat+3, "Incorrect size for appended internal "
            "state vector" );

    tk::real rho(0.0);
    for (std::size_t k=0; k<nmat; ++k)
      rho += ul[densityIdx(nmat, k)];

    auto ur = ul;

    // Internal cell velocity components
    auto v1l = ul[ncomp+velocityIdx(nmat, 0)];
    auto v2l = ul[ncomp+velocityIdx(nmat, 1)];
    auto v3l = ul[ncomp+velocityIdx(nmat, 2)];
    // Normal component of velocity
    auto vnl = v1l*fn[0] + v2l*fn[1] + v3l*fn[2];
    // Ghost state velocity components
    auto v1r = v1l - 2.0*vnl*fn[0];
    auto v2r = v2l - 2.0*vnl*fn[1];
    auto v3r = v3l - 2.0*vnl*fn[2];
    // Boundary condition
    for (std::size_t k=0; k<nmat; ++k)
    {
      ur[volfracIdx(nmat, k)] = ul[volfracIdx(nmat, k)];
      ur[densityIdx(nmat, k)] = ul[densityIdx(nmat, k)];
      ur[energyIdx(nmat, k)] = ul[energyIdx(nmat, k)];
    }
    ur[momentumIdx(nmat, 0)] = rho * v1r;
    ur[momentumIdx(nmat, 1)] = rho * v2r;
    ur[momentumIdx(nmat, 2)] = rho * v3r;

    // Internal cell primitive quantities using the separately reconstructed
    // primitive quantities. This is used to get ghost state for primitive
    // quantities

    // velocity
    ur[ncomp+velocityIdx(nmat, 0)] = v1r;
    ur[ncomp+velocityIdx(nmat, 1)] = v2r;
    ur[ncomp+velocityIdx(nmat, 2)] = v3r;
    // material pressures
    for (std::size_t k=0; k<nmat; ++k)
      ur[ncomp+pressureIdx(nmat, k)] = ul[ncomp+pressureIdx(nmat, k)];

    Assert( ur.size() == ncomp+nmat+3, "Incorrect size for appended boundary "
            "state vector" );

    return {{ std::move(ul), std::move(ur) }};
  }

  //! \brief Boundary state function providing the left and right state of a
  //!   face at farfield outlet boundaries
  //! \param[in] system Equation system index
  //! \param[in] ncomp Number of scalar components in this PDE system
  //! \param[in] ul Left (domain-internal) state
  //! \param[in] fn Unit face normal
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \details The farfield outlet boudary calculation, implemented here, is
  //!   based on the characteristic theory of hyperbolic systems. For subsonic
  //!   outlet flow, there is 1 incoming characteristic per material.
  //!   Therefore, we calculate the ghost cell state by taking material
  //!   pressure from the outside and other quantities from the internal cell.
  //!   For supersonic outlet flow, all the characteristics are from internal
  //!   cell and we obtain the ghost cell state from the internal cell.
  //! \note The function signature must follow tk::StateFn
  static tk::StateFn::result_type
  farfieldOutlet( ncomp_t system,
                  ncomp_t ncomp,
                  const std::vector< EOS >& mat_blk,
                  const std::vector< tk::real >& ul,
                  tk::real, tk::real, tk::real, tk::real,
                  const std::array< tk::real, 3 >& fn )
  {
    const auto nmat =
      g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[system];

    auto fp =
      g_inputdeck.get< tag::param, tag::multimat, tag::farfield_pressure >()[ system ];

    Assert( ul.size() == ncomp+nmat+3, "Incorrect size for appended internal "
            "state vector" );

    auto ur = ul;

    // Internal cell velocity components
    auto v1l = ul[ncomp+velocityIdx(nmat, 0)];
    auto v2l = ul[ncomp+velocityIdx(nmat, 1)];
    auto v3l = ul[ncomp+velocityIdx(nmat, 2)];

    // Normal velocity
    auto vn = v1l*fn[0] + v2l*fn[1] + v3l*fn[2];

    // Acoustic speed
    tk::real a(0.0);
    for (std::size_t k=0; k<nmat; ++k)
      if (ul[volfracIdx(nmat, k)] > 1.0e-04)
        a = std::max( a, mat_blk[k].eosCall< EOS::soundspeed >(
          ul[densityIdx(nmat, k)], ul[ncomp+pressureIdx(nmat, k)],
          ul[volfracIdx(nmat, k)] ) );

    // Mach number
    auto Ma = vn / a;

    if(Ma >= 0 && Ma < 1) {         // Subsonic outflow
      for (std::size_t k=0; k<nmat; ++k)
        ur[energyIdx(nmat, k)] = ul[volfracIdx(nmat, k)] *
        mat_blk[k].eosCall< EOS::totalenergy >(
          ur[densityIdx(nmat, k)]/ul[volfracIdx(nmat, k)], v1l, v2l, v3l, fp );

      // Internal cell primitive quantities using the separately reconstructed
      // primitive quantities. This is used to get ghost state for primitive
      // quantities

      // material pressures
      for (std::size_t k=0; k<nmat; ++k)
        ur[ncomp+pressureIdx(nmat, k)] = ul[volfracIdx(nmat, k)] * fp;
    }

    Assert( ur.size() == ncomp+nmat+3, "Incorrect size for appended boundary "
            "state vector" );

    return {{ std::move(ul), std::move(ur) }};
  }

  //! \brief Boundary state function providing the left and right state of a
  //!   face at extrapolation boundaries
  //! \param[in] ul Left (domain-internal) state
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \note The function signature must follow tk::StateFn. For multimat, the
  //!   left or right state is the vector of conserved quantities, followed by
  //!   the vector of primitive quantities appended to it.
  static tk::StateFn::result_type
  extrapolate( ncomp_t, ncomp_t,
               const std::vector< EOS >&,
               const std::vector< tk::real >& ul,
               tk::real, tk::real, tk::real, tk::real,
               const std::array< tk::real, 3 >& )
  {
    return {{ ul, ul }};
  }

} // inciter::

#endif // BCFunctions_h
