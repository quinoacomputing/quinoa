// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Physics/DGNavierStokes.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Physics policy for the Navier-Stokes equation using discontinuous
    Galerkin
  \details   This file defines a Physics policy class for the compressible
    single-material viscous flow equations using discontinuous Galerkin
    discretization, defined in PDE/CompFlow/DGCompFlow.h. The class defined here
    is used to configure the behavior of DGCompFlow. See
    PDE/CompFlow/Physics/DG.h for general requirements on Physics policy classes
    for DGCompFlow.
*/
// *****************************************************************************
#ifndef CompFlowPhysicsDGNavierStokes_h
#define CompFlowPhysicsDGNavierStokes_h

#include <array>
#include <limits>

#include "Types.h"
#include "Inciter/Options/Physics.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! CompFlow system of PDEs problem: Navier-Stokes (viscous flow)
//! \details This class adds the viscous force contributions to the momentum and
//!    energy conservation equations governing compressible flow.
class CompFlowPhysicsNavierStokes {

  public:
    //! Add viscous stress contribution to momentum and energy rhs
    //! \param[in] dt Size of time step
    //! \param[in] J Element Jacobi determinant
    //! \param[in] N Element node indices
    //! \param[in] grad Shape function derivatives, nnode*ndim [4][3]
    //! \param[in] u Solution at element nodes at recent time step
    //! \param[in] r Pointers to right hand side at component and offset
    //! \param[in,out] R Right-hand side vector contributing to
    static void
    viscousRhs( tk::real dt,
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
      auto c = dt * J/6.0 * mu;
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
      c = dt * J/6.0 * mu * 2.0/3.0;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t k=0; k<4; ++k) {
            R.var(r[i+1],N[0]) += c * grad[0][i]*grad[k][j]*u[j+1][k]/u[0][k];
            R.var(r[i+1],N[1]) += c * grad[1][i]*grad[k][j]*u[j+1][k]/u[0][k];
            R.var(r[i+1],N[2]) += c * grad[2][i]*grad[k][j]*u[j+1][k]/u[0][k];
            R.var(r[i+1],N[3]) += c * grad[3][i]*grad[k][j]*u[j+1][k]/u[0][k];
          }
      // add deviatoric viscous stress contribution to energy rhs
      c = dt * J/24.0 * mu;
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
      c = dt * J/24.0 * mu * 2.0/3.0;
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
    //! \param[in] u Solution at element nodes at recent time step
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
        auto dt = L * L * r / (2.0*mu); // dt ~ dx^2/nu = dx^2*rho/(2mu)
        if (dt < mindt) mindt = dt;
      }
      return mindt;
    }

    //! Add heat conduction contribution to energy rhs
    //! \param[in] dt Size of time step
    //! \param[in] J Element Jacobi determinant
    //! \param[in] N Element node indices
    //! \param[in] grad Shape function derivatives, nnode*ndim [4][3]
    //! \param[in] u Solution at element nodes at recent time step
    //! \param[in] r Pointers to right hand side at component and offset
    //! \param[in,out] R Right-hand side vector contributing to
    static void
    conductRhs( tk::real dt,
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
      auto c = dt * J/24.0 * kc;
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
    //! \param[in] u Solution at element nodes at recent time step
    //! \return Minimum time step size based on thermal diffusion
    static tk::real
    conduct_dt( tk::real L,
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

} // dg::

} // inciter::

#endif // CompFlowPhysicsDGNavierStokes_h
