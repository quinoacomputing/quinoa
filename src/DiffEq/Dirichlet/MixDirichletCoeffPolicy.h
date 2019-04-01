// *****************************************************************************
/*!
  \file      src/DiffEq/Dirichlet/MixDirichletCoeffPolicy.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     MixDirichlet coefficients policies
  \details   This file defines coefficients policy classes for the MixDirichlet
    SDE, defined in DiffEq/MixDirichlet.h.

    General requirements on Dirichlet SDE coefficients policy classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, b, S, and kappa. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_type ncomp,
          const std::vector< kw::sde_b::info::expect::type >& b_,
          const std::vector< kw::sde_S::info::expect::type >& S_,
          const std::vector< kw::sde_kappa::info::expect::type >& k_,
          const std::vector< kw::sde_r::info::expect::type >& rho_,
          const std::vector< kw::sde_r::info::expect::type >& r_,
          std::vector< kw::sde_b::info::expect::type  >& b,
          std::vector< kw::sde_S::info::expect::type >& S,
          std::vector< kw::sde_kappa::info::expect::type >& k,
          std::vector< kw::sde_r::info::expect::type >& r_ );
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of
        MixDirichlet SDEs.
      - Constant references to b_, S_, k_, rho_, r_, which denote three vectors
        of real values used to initialize the parameter vectors of the
        MixDirichlet SDEs. The length of the vectors must be equal to the number
        of components given by ncomp.
      - References to b, S, k, rho, r, which denote the parameter vectors to be
        initialized based on b_, S_, k_, rho_.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::CoeffPolicyType type() noexcept {
          return ctr::CoeffPolicyType::CONST_COEFF;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.

    - Must define the function _update()_, called from
      MixDirichlet::advance(), updating the model coefficients.
      Required signature:
      \code{.cpp}
        void update(
               char depvar,
               ncomp_t ncomp,
               const std::map< tk::ctr::Product, tk::real >& moments,
               const std::vector< kw::sde_rho::info::expect::type >& rho,
               const std::vector< kw::sde_r::info::expect::type >& r,
               const std::vector< kw::sde_kappa::info::expect::type >& kprime,
               std::vector< kw::sde_kappa::info::expect::type >& k,
               std::vector< kw::sde_kappa::info::expect::type >& S ) const {}
      \endcode
      where _depvar_ is the dependent variable associated with the mix
      Dirichlet SDE, specified in the control file by the user, _ncomp_
      is the number of components in the system, _moments_ is the map
      associating moment IDs (tk::ctr::vector< tk::ctr::Term >) to values of
      statistical moments, _rho_, _r_, and _kprime_ are user-defined
      parameters, and _k_ and _S_ are the SDE parameters computed, see
      DiffEq/DiffEq/MixDirichlet.h.
*/
// *****************************************************************************
#ifndef MixDirichletCoeffPolicy_h
#define MixDirichletCoeffPolicy_h

#include <brigand/sequences/list.hpp>

#include "Types.h"
#include "Walker/Options/CoeffPolicy.h"
#include "SystemComponents.h"

namespace walker {

//! Compute parameter vector r based on r_i = rho_N/rho_i - 1
std::vector< kw::sde_r::info::expect::type >
MixDir_r( const std::vector< kw::sde_rho::info::expect::type >& rho );

//! MixDirichlet coefficients policity: constants in time + mean(rho) = const
//! \details User-defined parameters b and kappaprime are constant vectors in
//!   time and, S is constrained to make \f$\mathrm{d}<rho>/\mathrm{d}t = 0\f$.
class MixDirichletHomCoeffConst {

  private:
    using ncomp_t = tk::ctr::ncomp_type;

  public:
    //! Constructor: initialize coefficients
    MixDirichletHomCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_b::info::expect::type >& b_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappa::info::expect::type >& kprime_,
      const std::vector< kw::sde_rho::info::expect::type >& rho_,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_kappa::info::expect::type >& kprime,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_rho::info::expect::type >& rho,
      std::vector< kw::sde_r::info::expect::type >& r,
      std::vector< kw::sde_kappa::info::expect::type >& k );

    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::HOMOGENEOUS; }

    //! Update coefficients
    void update(
      char depvar,
      ncomp_t ncomp,
      const std::map< tk::ctr::Product, tk::real >& moments,
      const std::vector< kw::sde_rho::info::expect::type >& rho,
      const std::vector< kw::sde_r::info::expect::type >& r,
      const std::vector< kw::sde_kappa::info::expect::type >& kprime,
      std::vector< kw::sde_kappa::info::expect::type >& k,
      std::vector< kw::sde_kappa::info::expect::type >& S ) const;
};

//! List of all MixDirichlet's coefficients policies
using MixDirichletCoeffPolicies = brigand::list< MixDirichletHomCoeffConst >;

} // walker::

#endif // MixDirichletCoeffPolicy_h
