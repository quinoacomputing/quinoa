// *****************************************************************************
/*!
  \file      src/PDE/CompFlowPhysics.h
  \author    J. Bakosi
  \date      Mon 29 Aug 2016 11:47:33 AM MDT
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

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Physics.h"

namespace inciter {

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
    //! \param[in] s Solution at element nodes at recent time step stage
    //! \param[in] r Pointers to right hand side at component and offset
    //! \param[in,out] R Right-hand side vector contributing to
    static void
    viscousRhs( tk::real mult,
                tk::real dt,
                tk::real J,
                const std::array< std::size_t, 4 >& N,
                const std::array< std::array< tk::real, 3 >, 4 >& grad,
                const std::array< std::array< tk::real, 4 >, 5 >& s,
                const std::array< const tk::real*, 5 >& r,
                tk::Fields& R )
    {
      // dynamic viscosity
      tk::real mu = g_inputdeck.get< tag::param, tag::compflow, tag::mu >()[0];

      // add deviatoric viscous stress contribution to momentum rhs
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t k=0; k<4; ++k) {
            R.var(r[i+1],N[0]) -= mult * dt * J/6.0 * mu *
                                  grad[0][j]*(grad[k][j]*s[i+1][k] +
                                              grad[k][i]*s[j+1][k])/s[0][k];
            R.var(r[i+1],N[1]) -= mult * dt * J/6.0 * mu *
                                  grad[1][j]*(grad[k][j]*s[i+1][k] +
                                              grad[k][i]*s[j+1][k])/s[0][k];
            R.var(r[i+1],N[2]) -= mult * dt * J/6.0 * mu *
                                  grad[2][j]*(grad[k][j]*s[i+1][k] +
                                              grad[k][i]*s[j+1][k])/s[0][k];
            R.var(r[i+1],N[3]) -= mult * dt * J/6.0 * mu *
                                  grad[3][j]*(grad[k][j]*s[i+1][k] +
                                              grad[k][i]*s[j+1][k])/s[0][k];
          }

      // add isotropic viscous stress contribution to momentum rhs
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t k=0; k<4; ++k) {
            R.var(r[i+1],N[0]) += mult * dt * J/6.0 * 2.0/3.0 * mu *
                                  grad[0][i]*grad[k][j]*s[j+1][k]/s[0][k];
            R.var(r[i+1],N[1]) += mult * dt * J/6.0 * 2.0/3.0 * mu *
                                  grad[1][i]*grad[k][j]*s[j+1][k]/s[0][k];
            R.var(r[i+1],N[2]) += mult * dt * J/6.0 * 2.0/3.0 * mu *
                                  grad[2][i]*grad[k][j]*s[j+1][k]/s[0][k];
            R.var(r[i+1],N[3]) += mult * dt * J/6.0 * 2.0/3.0 * mu *
                                  grad[3][i]*grad[k][j]*s[j+1][k]/s[0][k];
          }

      // add deviatoric viscous stress contribution to energy rhs
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t k=0; k<4; ++k) {
            R.var(r[4],N[0]) -= mult * dt * J/24.0 * mu * s[i+1][k]/s[0][k] *
                                grad[0][j]*(grad[k][j]*s[i+1][k] +
                                            grad[k][i]*s[j+1][k])/s[0][k];
            R.var(r[4],N[1]) -= mult * dt * J/24.0 * mu * s[i+1][k]/s[0][k] *
                                grad[1][j]*(grad[k][j]*s[i+1][k] +
                                            grad[k][i]*s[j+1][k])/s[0][k];
            R.var(r[4],N[2]) -= mult * dt * J/24.0 * mu * s[i+1][k]/s[0][k] *
                                grad[2][j]*(grad[k][j]*s[i+1][k] +
                                            grad[k][i]*s[j+1][k])/s[0][k];
            R.var(r[4],N[3]) -= mult * dt * J/24.0 * mu * s[i+1][k]/s[0][k] *
                                grad[3][j]*(grad[k][j]*s[i+1][k] +
                                            grad[k][i]*s[j+1][k])/s[0][k];
          }

      // add isotropic viscous stress contribution to energy rhs
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t k=0; k<4; ++k) {
            R.var(r[4],N[0]) += mult * dt * J/24.0 * s[i+1][k]/s[0][k] *
                                grad[0][i]*grad[k][j]*s[j+1][k]/s[0][k] *
                                2.0/3.0 * mu;
            R.var(r[4],N[1]) += mult * dt * J/24.0 * s[i+1][k]/s[0][k] *
                                grad[1][i]*grad[k][j]*s[j+1][k]/s[0][k] *
                                2.0/3.0 * mu;
            R.var(r[4],N[2]) += mult * dt * J/24.0 * s[i+1][k]/s[0][k] *
                                grad[2][i]*grad[k][j]*s[j+1][k]/s[0][k] *
                                2.0/3.0 * mu;
            R.var(r[4],N[3]) += mult * dt * J/24.0 * s[i+1][k]/s[0][k] *
                                grad[3][i]*grad[k][j]*s[j+1][k]/s[0][k] *
                                2.0/3.0 * mu;
          }
    }

    //! Add heat conduction contribution to energy rhs
    //! \param[in] mult Multiplier differentiating the different stages in
    //!    multi-stage time stepping
    //! \param[in] dt Size of time step
    //! \param[in] J Element Jacobi determinant
    //! \param[in] N Element node indices
    //! \param[in] grad Shape function derivatives, nnode*ndim [4][3]
    //! \param[in] s Solution at element nodes at recent time step stage
    //! \param[in] r Pointers to right hand side at component and offset
    //! \param[in,out] R Right-hand side vector contributing to
    static void
    conductRhs( tk::real mult,
                tk::real dt,
                tk::real J,
                const std::array< std::size_t, 4 >& N,
                const std::array< std::array< tk::real, 3 >, 4 >& grad,
                const std::array< std::array< tk::real, 4 >, 5 >& s,
                const std::array< const tk::real*, 5 >& r,
                tk::Fields& R )
    {
      // specific heat at constant volume
      tk::real cv = g_inputdeck.get< tag::param, tag::compflow, tag::cv >()[0];
      // thermal conductivity
      tk::real kc = g_inputdeck.get< tag::param, tag::compflow, tag::k >()[0];

      // compute temperature
      std::array< tk::real, 4 > T;
      for (std::size_t i=0; i<4; ++i)
        T[i] = cv*(s[4][i] - (s[1][i]*s[1][i] +
                              s[2][i]*s[2][i] +
                              s[3][i]*s[3][i])/2.0/s[0][i]) / s[0][i];

      // add heat conduction contribution to energy rhs
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t k=0; k<4; ++k) {
          R.var(r[4],N[0]) += mult * dt * J/24.0 * grad[k][i] * T[k] * kc;
          R.var(r[4],N[1]) += mult * dt * J/24.0 * grad[k][i] * T[k] * kc;
          R.var(r[4],N[2]) += mult * dt * J/24.0 * grad[k][i] * T[k] * kc;
          R.var(r[4],N[3]) += mult * dt * J/24.0 * grad[k][i] * T[k] * kc;
        }
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

    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::EULER; }
};

//! List of all CompFlow problem policies
using CompFlowPhysics = boost::mpl::vector< CompFlowPhysicsNavierStokes,
                                            CompFlowPhysicsEuler >;

} // inciter::

#endif // CompFlowPhysics_h
