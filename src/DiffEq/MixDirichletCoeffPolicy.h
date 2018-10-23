// *****************************************************************************
/*!
  \file      src/DiffEq/MixDirichletCoeffPolicy.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
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
          std::vector< kw::sde_b::info::expect::type  >& b,
          std::vector< kw::sde_S::info::expect::type >& S,
          std::vector< kw::sde_kappa::info::expect::type >& k )
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of
        MixDirichlet SDEs.
      - Constant references to b_, S_, and k_, which denote three vectors of
        real values used to initialize the parameter vectors of the MixDirichlet
        SDEs. The length of the vectors must be equal to the number of
        components given by ncomp.
      - References to b, S, and k, which denote the parameter vectors to be
        initialized based on b_, S_, and k_.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::CoeffPolicyType type() noexcept {
          return ctr::CoeffPolicyType::CONSTANT;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef MixDirichletCoeffPolicy_h
#define MixDirichletCoeffPolicy_h

#include <brigand/sequences/list.hpp>

#include "Types.h"
#include "Walker/Options/CoeffPolicy.h"
#include "SystemComponents.h"

namespace walker {

//! MixDirichlet coefficients policity: constants in time + <rho> = const
//! \details User-defined parameters b and kappa are constant vectors in time
//!   and, S is constrained to make \f$\mathrm{d}\<\rho\>/\mathrm{d}t = 0\f$.
class MixDirichletHomCoeffConst {

  private:
    using ncomp_t = tk::ctr::ncomp_type;

  public:
    //! Constructor: initialize coefficients
    MixDirichletHomCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_b::info::expect::type >& b_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappa::info::expect::type >& k_,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_kappa::info::expect::type >& k,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_rho2::info::expect::type >& rho2 )
    {
      ErrChk( b_.size() == ncomp,
              "Wrong number of MixDirichlet SDE parameters 'b'");
      ErrChk( S_.size() == ncomp,
              "Wrong number of MixDirichlet SDE parameters 'S'");
      ErrChk( k_.size() == ncomp,
              "Wrong number of MixDirichlet SDE parameters 'kappa'");

      b = b_;
      S = S_;
      k = k_;
      rho2 = rho2_;
    }

    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::HOMOGENEOUS; }

    //! \brief Update coefficients vector S so <rho> = const
    void update(
      char depvar,
      ncomp_t ncomp,
      const std::map< tk::ctr::Product, tk::real >& moments,
      const std::vector< kw::sde_rho2::info::expect::type >& /*rho2*/,
      std::vector< kw::sde_kappa::info::expect::type >& S ) const
    {
      using tk::ctr::lookup;
      using tk::ctr::mean;
      // statistics nomenclature:
      //   Y = instantaneous mass fraction,
      //   R = instantaneous density,
      //   y = Y - <Y>, mass fraction fluctuation about its mean,
      //   r = R - <R>, density fluctuation about its mean,
      // <Y> = mean mass fraction,
      // <R> = mean density,


      tk::real d = lookup( mean(depvar,ncomp), moments );       // <R>
      if (d < 1.0e-8) { std::cout << "d:" << d << " "; d = 0.5; }

      for (ncomp_t c=0; c<ncomp; ++c) {
        tk::real m = lookup( mean(depvar,c), moments );         // <Y>
        if (m<1.0e-8 || m>1.0-1.0e-8) m = 0.5;

        //S[c] = 0.5;

        if (S[c] < 0.0 || S[c] > 1.0) {
          std::cout << S[c] << " ";
          //S[c] = 0.5;
        }
      }
    }
};

//! List of all MixDirichlet's coefficients policies
using MixDirichletCoeffPolicies = brigand::list< MixDirichletHomCoeffConst >;

} // walker::

#endif // MixDirichletCoeffPolicy_h
