// *****************************************************************************
/*!
  \file      src/DiffEq/MixNumberFractionBetaCoeffPolicy.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Mix number-fraction beta SDE coefficients policies
  \details   This file defines coefficients policy classes for the mix
    number-fraction beta SDE, defined in DiffEq/MixNumberFractionBeta.h.

    General requirements on mix number-fraction beta SDE coefficients policy
    classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, b, S, kappa, rho2, and rcomma. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_type ncomp,
          const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
          const std::vector< kw::sde_S::info::expect::type >& S_,
          const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
          const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
          const std::vector< kw::sde_rcomma::info::expect::type >& rcomma_,
          std::vector< kw::sde_bprime::info::expect::type  >& bprime,
          std::vector< kw::sde_S::info::expect::type >& S,
          std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
          std::vector< kw::sde_rho2::info::expect::type >& rho2_,
          std::vector< kw::sde_rcomma::info::expect::type >& rcomma_ );
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of
        mix number-fraction beta SDEs.
      - Constant references to bprime_, S_, kprime_, rho2_, and rcomma_, which
        denote five vectors of real values used to initialize the parameter
        vectors of the system of mix number-fraction beta SDEs. The length of
        the vectors must be equal to the number of components given by ncomp.
      - References to bprime, S, kprime, rho2_, and rcomma, which denote the
        parameter vectors to be initialized based on bprime_, S_, kprime_,
        rho2_, and rcomma_.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::CoeffPolicyType type() noexcept {
          return ctr::CoeffPolicyType::DECAY;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.

    - Must define the function _update()_, called from
      MixNumberFractionBeta::advance(), updating the model coefficients.
      Required signature:
      \code{.cpp}
        void update(
          char depvar,
          ncomp_t ncomp,
          const std::map< tk::ctr::Product, tk::real >& moments,
          const std::vector< kw::sde_bprime::info::expect::type  >& bprime,
          const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
          std::vector< kw::sde_b::info::expect::type  >& b,
          std::vector< kw::sde_kappa::info::expect::type >& k ) const {}
      \endcode
      where _depvar_ is the dependent variable associated with the mix
      number-fraction beta SDE, specified in the control file by the user, _ncomp_ is
      the number of components in the system, _moments_ is the map associating
      moment IDs (tk::ctr::vector< tk::ctr::Term >) to values of statistical
      moments, _bprime_, _kprime_ are user-defined parameters, and _b_, _k_ are
      the SDE parameters computed, see DiffEq/MixNumberFractionBeta.h.
*/
// *****************************************************************************
#ifndef MixNumberFractionBetaCoeffPolicy_h
#define MixNumberFractionBetaCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Walker/Options/CoeffPolicy.h"
#include "SystemComponents.h"

namespace walker {

//! \brief Mix number-fraction beta SDE decay coefficients policity.
//! \details User-defined parameters b' and kappa' are constants in time and
//!   ensure decay in the evolution of <x^2>.
//! \author J. Bakosi
class MixNumFracBetaCoeffDecay {

    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor: initialize coefficients
    MixNumFracBetaCoeffDecay(
      ncomp_t ncomp,
      const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
      const std::vector< kw::sde_rcomma::info::expect::type >& rcomma_,
      std::vector< kw::sde_bprime::info::expect::type  >& bprime,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
      std::vector< kw::sde_rho2::info::expect::type >& rho2,
      std::vector< kw::sde_rcomma::info::expect::type >& rcomma,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_kappa::info::expect::type >& k )
    {
      ErrChk( bprime_.size() == ncomp,
        "Wrong number of mix number-fraction beta SDE parameters 'b''");
      ErrChk( S_.size() == ncomp,
        "Wrong number of mix number-fraction beta SDE parameters 'S'");
      ErrChk( kprime_.size() == ncomp,
        "Wrong number of mix number-fraction beta SDE parameters 'kappa''");
      ErrChk( rho2_.size() == ncomp,
        "Wrong number of mix number-fraction beta SDE parameters 'rho2'");
      ErrChk( rcomma_.size() == ncomp,
        "Wrong number of mix number-fraction beta SDE parameters 'rcomma'");

      bprime = bprime_;
      S = S_;
      kprime = kprime_;
      rho2 = rho2_;
      rcomma = rcomma_;

      b.resize( bprime.size() );
      k.resize( kprime.size() );
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::DECAY; }

    //! \brief Update coefficients using constant coefficients for b' and kappa'
    //! \details This where the mix number-fraction beta SDE is made consistent
    //!   with the no-mix and fully mixed limits by specifying the SDE
    //!   coefficients, b and kappa as functions of b' and kappa'.
    void update(
      char depvar,
      ncomp_t ncomp,
      const std::map< tk::ctr::Product, tk::real >& moments,
      const std::vector< kw::sde_bprime::info::expect::type  >& bprime,
      const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_kappa::info::expect::type >& k ) const
    {
      for (ncomp_t c=0; c<ncomp; ++c) {
        tk::real m = tk::ctr::lookup( tk::ctr::mean(depvar,c), moments );
        tk::real v = tk::ctr::lookup( tk::ctr::variance(depvar,c), moments );

        if (m<1.0e-8 || m>1.0-1.0e-8) m = 0.5;
        if (v<1.0e-8 && v>1.0-1.0e-8) v = 0.5;

        b[c] = bprime[c] * (1.0 - v / m / ( 1.0 - m ));
        k[c] = kprime[c] * v;
      }
    }
};

//! List of all mix numberf-fraction beta's coefficients policies
using MixNumFracBetaCoeffPolicies =
  boost::mpl::vector< MixNumFracBetaCoeffDecay >;

} // walker::

#endif // MixNumberFractionBetaCoeffPolicy_h
