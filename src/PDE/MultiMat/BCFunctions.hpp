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
#include "MiscMultiMatFns.hpp"

namespace inciter {

  //! \brief Boundary state function providing the left and right state of a
  //!   face at symmetry boundaries
  //! \param[in] ncomp Number of scalar components in this PDE system
  //! \param[in] ul Left (domain-internal) state
  //! \param[in] fn Unit face normal
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \note The function signature must follow tk::StateFn. For multimat, the
  //!   left or right state is the vector of conserved quantities, followed by
  //!   the vector of primitive quantities appended to it.
  static tk::StateFn::result_type
  symmetry( ncomp_t ncomp,
            const std::vector< EOS >&,
            const std::vector< tk::real >& ul,
            tk::real, tk::real, tk::real, tk::real,
            const std::array< tk::real, 3 >& fn )
  {
    auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
    const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

    [[maybe_unused]] auto nsld = numSolids(nmat, solidx);

    Assert( ul.size() == ncomp+nmat+3+nsld*6, "Incorrect size for appended "
            "internal state vector" );

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
      if (solidx[k] > 0) {
        // Internal inverse deformation tensor
        std::array< std::array< tk::real, 3 >, 3 > g;
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            g[i][j] = ul[deformIdx(nmat,solidx[k],i,j)];
        // Internal Cauchy stress tensor
        std::array< std::array< tk::real, 3 >, 3 > s;
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            s[i][j] = ul[ncomp+stressIdx(nmat,solidx[k],stressCmp[i][j])];
        // Make reflection matrix
        std::array< std::array< tk::real, 3 >, 3 >
        reflectionMat{{{1,0,0}, {0,1,0}, {0,0,1}}};
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            reflectionMat[i][j] -= 2*fn[i]*fn[j];
        // Reflect g
        g = tk::reflectTensor(g, reflectionMat);
        // Reflect s
        s = tk::reflectTensor(s, reflectionMat);
        // Copy g and s into ur
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j) {
            ur[deformIdx(nmat,solidx[k],i,j)] = g[i][j];
            ur[ncomp+stressIdx(nmat,solidx[k],stressCmp[i][j])] = s[i][j];
          }
      }
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

    Assert( ur.size() == ncomp+nmat+3+nsld*6, "Incorrect size for appended "
            "boundary state vector" );

    return {{ std::move(ul), std::move(ur) }};
  }

  //! \brief Boundary state function providing the left and right state of a
  //!   face at inlet boundaries
  //! \param[in] ncomp Number of scalar components in this PDE system
  //! \param[in] ul Left (domain-internal) state
  //! \param[in] fn Unit face normal
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \details The inlet boundary condition specifies a velocity at a
  //!   sideset and assumes a zero bulk pressure and density gradient
  //! \note The function signature must follow tk::StateFn
  static tk::StateFn::result_type
  inlet( ncomp_t ncomp,
            const std::vector< EOS >& mat_blk,
            const std::vector< tk::real >& ul,
            tk::real, tk::real, tk::real, tk::real,
            const std::array< tk::real, 3 >& fn )
  {
    auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
    const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();
    auto& inbc = g_inputdeck.get< tag::bc >()[0].get< tag::inlet >();

    // inlet velocity and material
    auto u_in = inbc[0].get< tag::velocity >();
    auto mat_in = inbc[0].get< tag::materialid >() - 1;
    auto p_in = inbc[0].get< tag::pressure >();
    auto t_in = inbc[0].get< tag::temperature >();

    [[maybe_unused]] auto nsld = numSolids(nmat, solidx);

    Assert( ul.size() == ncomp+nmat+3+nsld*6, "Incorrect size for appended "
            "internal state vector" );

    auto ur = ul;

    // External cell velocity, such that velocity = v_in at face
    auto v1r = u_in[0];
    auto v2r = u_in[1];
    auto v3r = u_in[2];

    // Normal inlet velocity
    auto vn = u_in[0]*fn[0] + u_in[1]*fn[1] + u_in[2]*fn[2];

    // Acoustic speed
    tk::real a(0.0);
    for (std::size_t k=0; k<nmat; ++k)
      if (ul[volfracIdx(nmat, k)] > 1.0e-04)
        a = std::max( a, mat_blk[k].compute< EOS::soundspeed >(
          ul[densityIdx(nmat, k)], ul[ncomp+pressureIdx(nmat, k)],
          ul[volfracIdx(nmat, k)], k ) );

    // Mach number
    auto Ma = vn / a;

    tk::real alphamin = 1e-12;
    tk::real pk(0.0);
    tk::real rho(0.0);
    for (std::size_t k=0; k<nmat; ++k) {
      if (k == mat_in)
        ur[volfracIdx(nmat,k)] = 1.0 -
          (static_cast< tk::real >(nmat-1))*alphamin;
      else
        ur[volfracIdx(nmat,k)] = alphamin;

      // Material pressure, which, for supersonic inflow, is the exterior
      // pressure and the interior pressure for subsonic
      if(Ma <= -1)
        pk = p_in;
      else
        pk = ul[ncomp+pressureIdx(nmat,k)]/ul[volfracIdx(nmat,k)];
      auto rhok = mat_blk[k].compute< EOS::density >(pk, t_in);

      ur[ncomp+pressureIdx(nmat, k)] = ur[volfracIdx(nmat,k)] * pk;
      ur[densityIdx(nmat,k)] = ur[volfracIdx(nmat,k)] * rhok;
      ur[energyIdx(nmat,k)] = ur[volfracIdx(nmat,k)] *
        mat_blk[k].compute< EOS::totalenergy >(rhok, v1r, v2r, v3r, pk);

      // bulk density
      rho += ur[densityIdx(nmat,k)];
    }

    ur[momentumIdx(nmat, 0)] = rho * v1r;
    ur[momentumIdx(nmat, 1)] = rho * v2r;
    ur[momentumIdx(nmat, 2)] = rho * v3r;

    // velocity
    ur[ncomp+velocityIdx(nmat, 0)] = v1r;
    ur[ncomp+velocityIdx(nmat, 1)] = v2r;
    ur[ncomp+velocityIdx(nmat, 2)] = v3r;

    Assert( ur.size() == ncomp+nmat+3+nsld*6, "Incorrect size for appended "
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
  farfield( ncomp_t ncomp,
            const std::vector< EOS >& mat_blk,
            const std::vector< tk::real >& ul,
            tk::real, tk::real, tk::real, tk::real,
            const std::array< tk::real, 3 >& fn )
  {
    auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
    const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

    // Farfield primitive quantities
    auto fp =
      g_inputdeck.get< tag::bc >()[0].get< tag::pressure >();
    auto ft =
      g_inputdeck.get< tag::bc >()[0].get< tag::temperature >();
    auto fu =
      g_inputdeck.get< tag::bc >()[0].get< tag::velocity >();
    auto fmat =
      g_inputdeck.get< tag::bc >()[0].get< tag::materialid >() - 1;

    [[maybe_unused]] auto nsld = numSolids(nmat, solidx);

    Assert( ul.size() == ncomp+nmat+3+nsld*6, "Incorrect size for appended "
            "internal state vector" );

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
        a = std::max( a, mat_blk[k].compute< EOS::soundspeed >(
          ul[densityIdx(nmat, k)], ul[ncomp+pressureIdx(nmat, k)],
          ul[volfracIdx(nmat, k)], k ) );

    // Mach number
    auto Ma = vn / a;

    tk::real alphamin = 1e-12;

    if (Ma <= -1) {  // Supersonic inflow
      // For supersonic inflow, all the characteristics are from outside.
      // Therefore, we calculate the ghost cell state using the primitive
      // variables from outside.
      tk::real rho(0.0);
      for (std::size_t k=0; k<nmat; ++k) {
        if (k == fmat)
          ur[volfracIdx(nmat,k)] = 1.0 -
            (static_cast< tk::real >(nmat-1))*alphamin;
        else
          ur[volfracIdx(nmat,k)] = alphamin;
        auto rhok = mat_blk[k].compute< EOS::density >(fp, ft);
        ur[densityIdx(nmat,k)] = ur[volfracIdx(nmat,k)] * rhok;
        ur[energyIdx(nmat,k)] = ur[volfracIdx(nmat,k)] *
          mat_blk[k].compute< EOS::totalenergy >(rhok, fu[0], fu[1], fu[2], fp);

        // material pressures
        ur[ncomp+pressureIdx(nmat, k)] = ur[volfracIdx(nmat, k)] * fp;

        rho += ur[densityIdx(nmat,k)];
      }
      for (std::size_t i=0; i<3; ++i) {
        ur[momentumIdx(nmat,i)] = rho * fu[i];
        ur[ncomp+velocityIdx(nmat, i)] = fu[i];
      }

    } else if (Ma > -1 && Ma < 0) {  // Subsonic inflow
      // For subsonic inflow, there is 1 outgoing characteristic and 4
      // incoming characteristics. Therefore, we calculate the ghost cell state
      // by taking pressure from the internal cell and other quantities from
      // the outside.
      tk::real rho(0.0);
      for (std::size_t k=0; k<nmat; ++k) {
        if (k == fmat)
          ur[volfracIdx(nmat,k)] = 1.0 -
            (static_cast< tk::real >(nmat-1))*alphamin;
        else
          ur[volfracIdx(nmat,k)] = alphamin;
        auto p = ul[ncomp+pressureIdx(nmat,k)] / ul[volfracIdx(nmat,k)];
        auto rhok = mat_blk[k].compute< EOS::density >(p, ft);
        ur[densityIdx(nmat,k)] = ur[volfracIdx(nmat,k)] * rhok;
        ur[energyIdx(nmat,k)] = ur[volfracIdx(nmat,k)] *
          mat_blk[k].compute< EOS::totalenergy >(rhok, fu[0], fu[1], fu[2], p);

        // material pressures
        ur[ncomp+pressureIdx(nmat, k)] = ur[volfracIdx(nmat, k)] * p;

        rho += ur[densityIdx(nmat,k)];
      }
      for (std::size_t i=0; i<3; ++i) {
        ur[momentumIdx(nmat,i)] = rho * fu[i];
        ur[ncomp+velocityIdx(nmat, i)] = fu[i];
      }

    } else if (Ma >= 0 && Ma < 1) {  // Subsonic outflow
      // For subsonic outflow, there is 1 incoming characteristic and 4
      // outgoing characteristics. Therefore, we calculate the ghost cell state
      // by taking pressure from the outside and other quantities from the
      // internal cell.
      for (std::size_t k=0; k<nmat; ++k) {
        ur[energyIdx(nmat, k)] = ul[volfracIdx(nmat, k)] *
        mat_blk[k].compute< EOS::totalenergy >(
          ur[densityIdx(nmat, k)]/ul[volfracIdx(nmat, k)], v1l, v2l, v3l, fp );

        // material pressures
        ur[ncomp+pressureIdx(nmat, k)] = ul[volfracIdx(nmat, k)] * fp;
      }
    }
    // Otherwise, for supersonic outflow, all the characteristics are from
    // internal cell. Therefore, we calculate the ghost cell state using the
    // conservative variables from internal cell (which is what ur is
    // initialized to).

    Assert( ur.size() == ncomp+nmat+3+nsld*6, "Incorrect size for appended "
            "boundary state vector" );

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
  //! \param[in] fn Unit face normal
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \note The function signature must follow tk::StateFn. For multimat, the
  //!   left or right state is the vector of conserved quantities, followed by
  //!   the vector of primitive quantities appended to it.
  static tk::StateFn::result_type
  noslipwall( ncomp_t ncomp,
              const std::vector< EOS >&,
              const std::vector< tk::real >& ul,
              tk::real, tk::real, tk::real, tk::real,
              const std::array< tk::real, 3 >& fn )
  {
    auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
    const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

    [[maybe_unused]] auto nsld = numSolids(nmat, solidx);

    Assert( ul.size() == ncomp+nmat+3+nsld*6, "Incorrect size for appended "
            "internal state vector" );

    tk::real rho(0.0);
    for (std::size_t k=0; k<nmat; ++k)
      rho += ul[densityIdx(nmat, k)];

    auto ur = ul;

    // Internal cell velocity components
    auto v1l = ul[ncomp+velocityIdx(nmat, 0)];
    auto v2l = ul[ncomp+velocityIdx(nmat, 1)];
    auto v3l = ul[ncomp+velocityIdx(nmat, 2)];
    // Ghost state velocity components
    auto v1r = -v1l;
    auto v2r = -v2l;
    auto v3r = -v3l;
    // Boundary condition
    for (std::size_t k=0; k<nmat; ++k)
    {
      ur[volfracIdx(nmat, k)] = ul[volfracIdx(nmat, k)];
      ur[densityIdx(nmat, k)] = ul[densityIdx(nmat, k)];
      ur[energyIdx(nmat, k)] = ul[energyIdx(nmat, k)];
      if (solidx[k] > 0) {
        // Internal inverse deformation tensor
        std::array< std::array< tk::real, 3 >, 3 > g;
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            g[i][j] = ul[deformIdx(nmat,solidx[k],i,j)];
        // Internal Cauchy stress tensor
        std::array< std::array< tk::real, 3 >, 3 > s;
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            s[i][j] = ul[ncomp+stressIdx(nmat,solidx[k],stressCmp[i][j])];
        // Make reflection matrix
        std::array< std::array< tk::real, 3 >, 3 >
        reflectionMat{{{1,0,0}, {0,1,0}, {0,0,1}}};
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            reflectionMat[i][j] -= 2*fn[i]*fn[j];
        // Reflect g
        g = tk::reflectTensor(g, reflectionMat);
        // Reflect s
        s = tk::reflectTensor(s, reflectionMat);
        // Copy g and s into ur
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j) {
            ur[deformIdx(nmat,solidx[k],i,j)] = g[i][j];
            ur[ncomp+stressIdx(nmat,solidx[k],stressCmp[i][j])] = s[i][j];
          }
      }
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

    Assert( ur.size() == ncomp+nmat+3+nsld*6, "Incorrect size for appended "
            "boundary state vector" );

    return {{ std::move(ul), std::move(ur) }};
  }

  //! \brief Boundary state function providing the left and right state of a
  //!   face at back pressure boundaries
  //! \param[in] ncomp Number of scalar components in this PDE system
  //! \param[in] ul Left (domain-internal) state
  //! \param[in] fn Unit face normal
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \details The back pressure boundary calculation, implemented here, is
  //!   based on the characteristic theory of hyperbolic systems.
  //! \note The function signature must follow tk::StateFn
  static tk::StateFn::result_type
  back_pressure( ncomp_t ncomp,
                 const std::vector< EOS >& mat_blk,
                 const std::vector< tk::real >& ul,
                 tk::real, tk::real, tk::real, tk::real,
                 const std::array< tk::real, 3 >& fn )
  {
    auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
    const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

    // Back pressure
    auto fbp = g_inputdeck.get< tag::bc >()[0].get< tag::back_pressure >().get<
      tag::pressure >();

    [[maybe_unused]] auto nsld = numSolids(nmat, solidx);

    Assert( ul.size() == ncomp+nmat+3+nsld*6, "Incorrect size for appended "
            "internal state vector" );

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
        a = std::max( a, mat_blk[k].compute< EOS::soundspeed >(
          ul[densityIdx(nmat, k)], ul[ncomp+pressureIdx(nmat, k)],
          ul[volfracIdx(nmat, k)], k ) );

    // Mach number
    auto Ma = vn / a;

    if (Ma < 1) {  // Subsonic outflow
      // For subsonic outflow, there is 1 incoming characteristic and 4
      // outgoing characteristics. Therefore, we calculate the ghost cell state
      // by taking pressure from the outside (i.e. the back pressure) and other
      // quantities from the internal cell.
      for (std::size_t k=0; k<nmat; ++k) {
        ur[energyIdx(nmat, k)] = ul[volfracIdx(nmat, k)] *
        mat_blk[k].compute< EOS::totalenergy >(
          ur[densityIdx(nmat, k)]/ul[volfracIdx(nmat, k)], v1l, v2l, v3l, fbp );

        // material pressures
        ur[ncomp+pressureIdx(nmat, k)] = ul[volfracIdx(nmat, k)] * fbp;
      }
    }
    // Otherwise, for supersonic outflow, all the characteristics are from
    // internal cell. Therefore, we calculate the ghost cell state using the
    // conservative variables from internal cell (which is what ur is
    // initialized to).

    Assert( ur.size() == ncomp+nmat+3+nsld*6, "Incorrect size for appended "
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
  //! \note The function signature must follow tk::StateFn. For multimat, the
  //!   left or right state is the vector of gradients of primitive quantities.
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
  //! \note The function signature must follow tk::StateFn. For multimat, the
  //!   left or right state is the vector of gradients of primitive quantities.
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

  //! \brief Boundary gradient function for zero gradient cells
  //! \param[in] ncomp Number of variables whos gradients are needed
  //! \param[in] dul Left (domain-internal) gradients
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \note The function signature must follow tk::StateFn. For multimat, the
  //!   left or right state is the vector of gradients of primitive quantities.
  static tk::StateFn::result_type
  zeroGrad( ncomp_t ncomp,
                const std::vector< EOS >&,
                const std::vector< tk::real >& dul,
                tk::real, tk::real, tk::real, tk::real,
                const std::array< tk::real, 3 >& )
  {
    Assert(dul.size() == 3*ncomp, "Incorrect size of boundary gradient vector");

    std::vector< tk::real > dur(3*ncomp, 0.0);

    return {{ std::move(dul), std::move(dur) }};
  }

} // inciter::

#endif // BCFunctions_h
