// *****************************************************************************
/*!
  \file      src/PDE/CompFlowPhysics.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Physics configurations for the compressible flow equations
  \details   This file defines policy classes for the compressible flow
    equations, defined in PDE/CompFlow.h.

    General requirements on flow equations problem policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::PhysicsType type() noexcept {
          return ctr::PhysicsType::NAVIERSTOKES;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef CompFlowPhysics_h
#define CompFlowPhysics_h

#include <vector>
#include <array>
#include <numeric>
#include <limits>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Physics.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! CompFlow system of PDEs problem: Navier-Stokes (viscous flow)
//! \details This class adds the viscous force contributions to the momentum and
//!    energy conservation equations governing compressible flow.
class CompFlowPhysicsNavierStokes {

  public:
    //! Add viscous stress contribution to momentum and energy rhs
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
    viscousRhs( tk::real mult,
                tk::real dt,
                tk::real J,
                const std::array< std::size_t, 4 >& N,
                const std::array< std::array< tk::real, 3 >, 4 >& grad,
                const std::array< std::array< tk::real, 4 >, 5 >& u,
                const std::array< const tk::real*, 5 >& r,
                tk::Fields& R )
    {
      // dynamic viscosity
      auto mu = g_inputdeck.get< tag::param, tag::compflow, tag::mu >()[0];
      // add deviatoric viscous stress contribution to momentum rhs
      auto c = mult * dt * J/6.0 * mu;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t k=0; k<4; ++k) {
            R.var(r[i+1],N[0]) -= c * grad[0][j]*(grad[k][j]*u[i+1][k] +
                                                  grad[k][i]*u[j+1][k])/u[0][k];
            R.var(r[i+1],N[1]) -= c * grad[1][j]*(grad[k][j]*u[i+1][k] +
                                                  grad[k][i]*u[j+1][k])/u[0][k];
            R.var(r[i+1],N[2]) -= c * grad[2][j]*(grad[k][j]*u[i+1][k] +
                                                  grad[k][i]*u[j+1][k])/u[0][k];
            R.var(r[i+1],N[3]) -= c * grad[3][j]*(grad[k][j]*u[i+1][k] +
                                                  grad[k][i]*u[j+1][k])/u[0][k];
          }
      // add isotropic viscous stress contribution to momentum rhs
      c = mult * dt * J/6.0 * mu * 2.0/3.0;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t k=0; k<4; ++k) {
            R.var(r[i+1],N[0]) += c * grad[0][i]*grad[k][j]*u[j+1][k]/u[0][k];
            R.var(r[i+1],N[1]) += c * grad[1][i]*grad[k][j]*u[j+1][k]/u[0][k];
            R.var(r[i+1],N[2]) += c * grad[2][i]*grad[k][j]*u[j+1][k]/u[0][k];
            R.var(r[i+1],N[3]) += c * grad[3][i]*grad[k][j]*u[j+1][k]/u[0][k];
          }
      // add deviatoric viscous stress contribution to energy rhs
      c = mult * dt * J/24.0 * mu;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t k=0; k<4; ++k) {
            R.var(r[4],N[0]) -= c * u[i+1][k]/u[0][k] *
                                grad[0][j]*(grad[k][j]*u[i+1][k] +
                                            grad[k][i]*u[j+1][k])/u[0][k];
            R.var(r[4],N[1]) -= c * u[i+1][k]/u[0][k] *
                                grad[1][j]*(grad[k][j]*u[i+1][k] +
                                            grad[k][i]*u[j+1][k])/u[0][k];
            R.var(r[4],N[2]) -= c * u[i+1][k]/u[0][k] *
                                grad[2][j]*(grad[k][j]*u[i+1][k] +
                                            grad[k][i]*u[j+1][k])/u[0][k];
            R.var(r[4],N[3]) -= c * u[i+1][k]/u[0][k] *
                                grad[3][j]*(grad[k][j]*u[i+1][k] +
                                            grad[k][i]*u[j+1][k])/u[0][k];
          }
      // add isotropic viscous stress contribution to energy rhs
      c = mult * dt * J/24.0 * mu * 2.0/3.0;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t k=0; k<4; ++k) {
            R.var(r[4],N[0]) += c * u[i+1][k]/u[0][k] *
                                grad[0][i]*grad[k][j]*u[j+1][k]/u[0][k];
            R.var(r[4],N[1]) += c * u[i+1][k]/u[0][k] *
                                grad[1][i]*grad[k][j]*u[j+1][k]/u[0][k];
            R.var(r[4],N[2]) += c * u[i+1][k]/u[0][k] *
                                grad[2][i]*grad[k][j]*u[j+1][k]/u[0][k];
            R.var(r[4],N[3]) += c * u[i+1][k]/u[0][k] *
                                grad[3][i]*grad[k][j]*u[j+1][k]/u[0][k];
          }
    }

    //! Compute the minimum time step size based on the viscous force
    //! \param[in] L Characteristic length scale
    //! \param[in] u Solution at element nodes at recent time step stage
    //! \return Minimum time step size based on viscous force
    static tk::real
    viscous_dt( tk::real L,
                const std::array< std::array< tk::real, 4 >, 5 >& u )
    {
      // dynamic viscosity
      auto mu = g_inputdeck.get< tag::param, tag::compflow, tag::mu >()[0];
      // compute the minimum viscous time step size across the four nodes
      tk::real mindt = std::numeric_limits< tk::real >::max();
      for (std::size_t j=0; j<4; ++j) {
        auto& r = u[0][j];              // rho
        auto dt = L * L * r / mu;       // dt ~ dx^2/nu = dx^2*rho/mu
        if (dt < mindt) mindt = dt;
      }
      return mindt;
    }

    //! Add heat conduction contribution to energy rhs
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
    conductRhs( tk::real mult,
                tk::real dt,
                tk::real J,
                const std::array< std::size_t, 4 >& N,
                const std::array< std::array< tk::real, 3 >, 4 >& grad,
                const std::array< std::array< tk::real, 4 >, 5 >& u,
                const std::array< const tk::real*, 5 >& r,
                tk::Fields& R )
    {
      // specific heat at constant volume
      auto cv = g_inputdeck.get< tag::param, tag::compflow, tag::cv >()[0];
      // thermal conductivity
      auto kc = g_inputdeck.get< tag::param, tag::compflow, tag::k >()[0];
      // compute temperature
      std::array< tk::real, 4 > T;
      for (std::size_t i=0; i<4; ++i)
        T[i] = cv*(u[4][i] - (u[1][i]*u[1][i] +
                              u[2][i]*u[2][i] +
                              u[3][i]*u[3][i])/2.0/u[0][i]) / u[0][i];
      // add heat conduction contribution to energy rhs
      auto c = mult * dt * J/24.0 * kc;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t k=0; k<4; ++k) {
          R.var(r[4],N[0]) += c * grad[k][i] * T[k];
          R.var(r[4],N[1]) += c * grad[k][i] * T[k];
          R.var(r[4],N[2]) += c * grad[k][i] * T[k];
          R.var(r[4],N[3]) += c * grad[k][i] * T[k];
        }
    }

    //! Compute the minimum time step size based on thermal diffusion
    //! \param[in] L Characteristic length scale
    //! \param[in] u Solution at element nodes at recent time step stage
    //! \return Minimum time step size based on thermal diffusion
    static tk::real
    heat_diffusion_dt( tk::real L,
                       const std::array< std::array< tk::real, 4 >, 5 >& u )
    {
      // ratio of specific heats
      auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0];
      // specific heat at constant volume
      auto cv = g_inputdeck.get< tag::param, tag::compflow, tag::cv >()[0];
      // thermal conductivity
      auto k = g_inputdeck.get< tag::param, tag::compflow, tag::k >()[0];
      // specific heat at constant pressure
      auto cp = g * cv;
      // compute the minimum conduction time step size across the four nodes
      tk::real mindt = std::numeric_limits< tk::real >::max();
      for (std::size_t j=0; j<4; ++j) {
        auto& r = u[0][j];              // rho
        auto dt = L * L * r * cp / k;   // dt ~ dx^2/alpha = dx^2*rho*C_p/k
        if (dt < mindt) mindt = dt;
      }
      return mindt;
    }

    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::NAVIERSTOKES; }
};

//! CompFlow system of PDEs problem: Euler (inviscid flow)
//! \details This class is a no-op, consistent with no additional physics needed
//!   to make the basic implementation in CompFlow the Euler equations governing
//!   compressible flow.
class CompFlowPhysicsEuler {

  public:
    //! Add viscous stress contribution to momentum and energy rhs (no-op)
    static void
    viscousRhs( tk::real,
                tk::real,
                tk::real,
                const std::array< std::size_t, 4 >&,
                const std::array< std::array< tk::real, 3 >, 4 >&,
                const std::array< std::array< tk::real, 4 >, 5 >&,
                const std::array< const tk::real*, 5 >&,
                tk::Fields& ) {}

    //! Compute the minimum time step size based on the viscous force
    //! \return A large time step size, i.e., ignore
    static tk::real
    viscous_dt( tk::real, const std::array< std::array< tk::real, 4 >, 5 >& )
    { return std::numeric_limits< tk::real >::max(); }

    //! Add heat conduction contribution to energy rhs (no-op)
    static void
    conductRhs( tk::real,
                tk::real,
                tk::real,
                const std::array< std::size_t, 4 >&,
                const std::array< std::array< tk::real, 3 >, 4 >&,
                const std::array< std::array< tk::real, 4 >, 5 >&,
                const std::array< const tk::real*, 5 >&,
                tk::Fields& ) {}

    //! Compute the minimum time step size based on thermal diffusion
    //! \return A large time step size, i.e., ignore
    static tk::real
    heat_diffusion_dt( tk::real,
                       const std::array< std::array< tk::real, 4 >, 5 >& )
    { return std::numeric_limits< tk::real >::max(); }

    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::EULER; }
};

//! List of all CompFlow problem policies
using CompFlowPhysics = boost::mpl::vector< CompFlowPhysicsNavierStokes,
                                            CompFlowPhysicsEuler >;

} // inciter::

#endif // CompFlowPhysics_h
