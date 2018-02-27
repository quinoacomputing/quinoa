// *****************************************************************************
/*!
  \file      src/PDE/Transport/Physics/CGAdvDiff.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Physics policy for advection-diffusion using continuous Galerkin
  \details   This file defines a Physics policy class for scalar transport using
     continuous Galerkin discretization, defined in PDE/Transport/CGTransport.h.
     The class defined here is used to configure the behavior of CGTransport.
     See PDE/Transport/Physics/CG.h for general requirements on Physics policy
     classes for CGTransport.
*/
// *****************************************************************************
#ifndef TransportPhysicsCGAdvDiff_h
#define TransportPhysicsCGAdvDiff_h

#include <limits>

#include "Types.h"
#include "Inciter/Options/Physics.h"

namespace inciter {
namespace cg {

//! Physics policy for advection-diffusion using continuous Galerkin
class TransportPhysicsAdvDiff {

  public:
    //! Add diffusion contribution to rhs at 2nd time step stage
    //! \param[in] e Equation system index, i.e., which transport equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this PDE
    //! \param[in] deltat Size of time step
    //! \param[in] J Element Jacobi determinant
    //! \param[in] grad Shape function derivatives, nnode*ndim [4][3]
    //! \param[in] N Element node indices
    //! \param[in] u Solution at element nodes at recent time step
    //! \param[in] r Pointers to right hand side at component and offset
    //! \param[in,out] R Right-hand side vector contributing to
    static void
    diffusionRhs( tk::ctr::ncomp_type e,
                  tk::ctr::ncomp_type ncomp,
                  tk::real deltat,
                  tk::real J,
                  const std::array< std::array< tk::real, 3 >, 4 >& grad,
                  const std::array< std::size_t, 4 >& N,
                  const std::vector< std::array< tk::real, 4 > >& u,
                  const std::vector< const tk::real* >& r,
                  tk::Fields& R )
    {
      // diffusivities for all components
      const auto& diff =
        g_inputdeck.get< tag::param, tag::transport, tag::diffusivity >().at(e);
      // add diffusion contribution to right hand side
      const auto d = deltat * J/6.0;
      for (ncomp_t c=0; c<ncomp; ++c)
        for (std::size_t k=0; k<3; ++k) {
          const auto D = diff[ 3*c+k ];
          for (std::size_t a=0; a<4; ++a)
            for (std::size_t b=0; b<4; ++b)
              R.var(r[c],N[a]) -= d * D * grad[a][k] * grad[b][k] * u[c][b];
        }
    }

    //! Compute the minimum time step size based on the diffusion
    //! \param[in] e Equation system index, i.e., which transport equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this PDE
    //! \param[in] L Characteristic length scale
    //! \return Minimum time step size based on diffusion
    static tk::real
    diffusion_dt( tk::ctr::ncomp_type e,
                  tk::ctr::ncomp_type ncomp,
                  tk::real L,
                  const std::vector< std::array< tk::real, 4 > >& )
    {
      // diffusivities for all components
      const auto& df =
        g_inputdeck.get< tag::param, tag::transport, tag::diffusivity >().at(e);
      // compute the minimum diffusion time step size across the four nodes
      tk::real mindt = std::numeric_limits< tk::real >::max();
      for (ncomp_t c=0; c<ncomp; ++c) {
        const auto di = 3*c;
        const auto d = std::max( df[di+2], std::max( df[di+0], df[di+1] ) );
        const auto dt = L * L / (2.0*d);  // dt ~ dx^2/(2D)
        if (dt < mindt) mindt = dt;
      }
      return mindt;
    }

    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::ADVDIFF; }
};

} // cg::
} // inciter::

#endif // TransportPhysicsCGAdvDiff_h
