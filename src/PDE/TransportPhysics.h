// *****************************************************************************
/*!
  \file      src/PDE/TransportPhysics.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Physics configurations for a system of transport equations
  \details   This file defines policy classes for transport equations,
    defined in PDE/Transport.h.

    General requirements on transport equation physics policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::PhysicsType type() noexcept {
          return ctr::PhysicsType::Advection;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef TransportPhysics_h
#define TransportPhysics_h

#include <limits>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Physics.h"

namespace inciter {

//! Transport equation system of PDEs problem: advection
class TransportPhysicsAdvection {

  public:
    //! Add diffusion contribution to rhs (no-op)
    static void
    diffusionRhs( tk::ctr::ncomp_type,
                  tk::ctr::ncomp_type,
                  tk::real,
                  tk::real,
                  tk::real,
                  const std::array< std::size_t, 4 >&,
                  const std::array< std::array< tk::real, 3 >, 4 >&,
                  const std::vector< std::array< tk::real, 4 > >&,
                  const std::vector< const tk::real* >&,
                  tk::Fields& ) {}

    //! Compute the minimum time step size based on the diffusion
    //! \return A large time step size, i.e., ignore
    static tk::real
    diffusion_dt( tk::ctr::ncomp_type,
                  tk::ctr::ncomp_type,
                  tk::real,
                  const std::vector< std::array< tk::real, 4 > >& )
    { return std::numeric_limits< tk::real >::max(); }

    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::ADVECTION; }
};

//! Transport equation system of PDEs problem: advection + diffusion
class TransportPhysicsAdvDiff {

  public:
    //! Add diffusion contribution to rhs
    //! \param[in] e Equation system index, i.e., which transport equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this PDE
    //! \param[in] mult Multiplier differentiating the different stages in
    //!    multi-stage time stepping
    //! \param[in] dt Size of time step
    //! \param[in] J Element Jacobi determinant
    //! \param[in] N Element node indices
    //! \param[in] grad Shape function derivatives, nnode*ndim [4][3]
    //! \param[in] u Solution at element nodes at recent time step stage
    //! \param[in] r Pointers to right hand side at component and offset
    //! \param[in,out] R Right-hand side vector contributing to
    static void
    diffusionRhs( tk::ctr::ncomp_type e,
                  tk::ctr::ncomp_type ncomp,
                  tk::real mult,
                  tk::real dt,
                  tk::real J,
                  const std::array< std::size_t, 4 >& N,
                  const std::array< std::array< tk::real, 3 >, 4 >& grad,
                  const std::vector< std::array< tk::real, 4 > >& u,
                  const std::vector< const tk::real* >& r,
                  tk::Fields& R )
    {
      // diffusivities for all components
      const auto& diff =
        g_inputdeck.get< tag::param, tag::transport, tag::diffusivity >().at(e);
      // add diffusion contribution to right hand side
      for (ncomp_t c=0; c<ncomp; ++c) {
        tk::real a = mult * dt * diff[c] * J;
        for (std::size_t i=0; i<4; ++i)
          for (std::size_t j=0; j<4; ++j)
            for (std::size_t k=0; k<3; ++k)
              R.var(r[c],N[j]) -= a * grad[j][k] * grad[i][k] * u[c][i];
      }
    }

    //! Compute the minimum time step size based on the diffusion
    //! \param[in] e Equation system index, i.e., which transport equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this PDE
    //! \param[in] L Characteristic length scale
    //! \param[in] u Solution at element nodes at recent time step stage
    //! \return Minimum time step size based on diffusion
    static tk::real
    diffusion_dt( tk::ctr::ncomp_type e,
                  tk::ctr::ncomp_type ncomp,
                  tk::real L,
                  const std::vector< std::array< tk::real, 4 > >& u )
    {
      // diffusivities for all components
      const auto& diff =
        g_inputdeck.get< tag::param, tag::transport, tag::diffusivity >().at(e);
      // compute the minimum diffusion time step size across the four nodes
      tk::real mindt = std::numeric_limits< tk::real >::max();
      for (ncomp_t c=0; c<ncomp; ++c)
        for (std::size_t j=0; j<4; ++j) {
          auto dt = L * L * u[c][j] / diff[c];  // dt ~ dx^2/D
          if (dt < mindt) mindt = dt;
        }
      return mindt;
    }

    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::ADVDIFF; }
};
//! List of all Transport equation problem policies
using TransportPhysics = boost::mpl::vector< TransportPhysicsAdvection
                                           , TransportPhysicsAdvDiff >;

} // inciter::

#endif // AdvDifTransport_h
