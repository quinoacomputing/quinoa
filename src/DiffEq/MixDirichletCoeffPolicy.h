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
//!   and, S is constrained to make \f$\mathrm{d}<rho>/\mathrm{d}t = 0\f$.
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
      using tk::ctr::Term;
      using tk::ctr::Moment;
      using tk::ctr::Product;

      // note: ncomp = K = N-1

      // statistics nomenclature:
      //   Y = instantaneous mass fraction
      //   R = instantaneous density
      //   y = Y - <Y>, mass fraction fluctuation about its mean
      //   r = R - <R>, density fluctuation about its mean
      // <Y> = mean mass fraction
      // <R> = mean density

      // <R>
      tk::real R = lookup( mean(depvar,ncomp), moments );
      if (R < 1.0e-8) R = 1.0;

      // <RY^c>
      std::vector< tk::real > RY( ncomp, 0.0 );
      for (ncomp_t c=0; c<ncomp; ++c) {
        Term tR( 'Y', ncomp, Moment::ORDINARY );
        Term tY( 'Y', c, Moment::ORDINARY );
        RY[c] = lookup( Product({tR,tY}), moments );
      }

      // Reynolds means

      // Reynolds means, Yc
      std::vector< tk::real > Y( ncomp, 0.0 );
      for (ncomp_t c=0; c<ncomp; ++c) {
        Y[c] = lookup( mean(depvar,c), moments );
        std::cout << "Y: " << Y[c] << ' ';
      }
      std::cout << std::endl;

      // sum of Yc
      tk::real sumY = 0.0;
      for (ncomp_t c=0; c<ncomp; ++c) sumY += Y[c];

      // Y|Kc
      std::vector< tk::real > YK( ncomp, 0.0 );
      for (ncomp_t c=0; c<ncomp; ++c) {
        YK[c] = sumY - lookup( mean(depvar,c), moments );
        std::cout << "YK: " << YK[c] << ' ';
      }
      std::cout << std::endl;

      // Favre means

      // Ytc
      std::vector< tk::real > Yt( ncomp, 0.0 );
      for (ncomp_t c=0; c<ncomp; ++c) {
        Yt[c] = RY[c] / R;
        std::cout << "Yt: " << Yt[c] << ' ';
      }
      std::cout << std::endl;

      // sum of Ytc
      tk::real sumYt = 0.0;
      for (ncomp_t c=0; c<ncomp; ++c) sumYt += Yt[c];

      // Yt|Kc
      std::vector< tk::real > YtK( ncomp, 0.0 );
      for (ncomp_t c=0; c<ncomp; ++c) {
        YtK[c] = sumYt - Yt[c];
        std::cout << "YtK: " << YtK[c] << ' ';
      }
      std::cout << std::endl;

      // Sc
      for (ncomp_t c=0; c<ncomp; ++c) {
        // 1st attempt at forcing <rho> = const
        //S[c] = 1.0/(1.0-YK[c]) - (1.0-Yt[c])/(1.0-YtK[c]);
        // 2nd attempt at forcing <rho> = const
        S[c] = YK[c]/(1.0-YK[c]) - (1.0-Yt[c])*YtK[c]/(1.0-YtK[c]) + Yt[c];
        //std::cout << "S: " << S[c] << ", YKc: " << YK[c]
        //          << ", Ytc: " << Yt[c] << ", YtKc: " << YtK[c] << ' ';
      }
      //std::cout << std::endl;

      for (ncomp_t c=0; c<ncomp; ++c) {
        if (S[c] < 0.0 || S[c] > 1.0) {
          std::cout << "S[" << c << "] bounds violated: " << S[c] << ' ';
          S[c] = 0.5;
        }
      }
      std::cout << std::endl;
    }
};

//! List of all MixDirichlet's coefficients policies
using MixDirichletCoeffPolicies = brigand::list< MixDirichletHomCoeffConst >;

} // walker::

#endif // MixDirichletCoeffPolicy_h
