//******************************************************************************
/*!
  \file      src/DiffEq/MixBetaCoeffPolicy.h
  \author    J. Bakosi
  \date      Thu 19 Mar 2015 09:14:02 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Mix beta SDE coefficients policies
  \details   This file defines coefficients policy classes for the mix beta SDE,
    defined in DiffEq/MixBeta.h.

    General requirements on mix beta SDE coefficients policy classes:

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
        mix beta SDEs.
      - Constant references to bprime_, S_, kprime_, rho2_, and rcomma_, which
        denote five vectors of real values used to initialize the parameter
        vectors of the system of mix beta SDEs. The length of the vectors must
        be equal to the number of components given by ncomp.
      - References to bprime, S, kprime, rho2_, and rcomma, which denote the
        parameter vectors to be initialized based on bprime_, S_, kprime_,
        rho2_, and rcomma_.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static tk::ctr::CoeffPolicyType type() noexcept {
          return tk::ctr::CoeffPolicyType::CONSTANT;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.

    - Must define the function _update()_, called from MixBeta::advance(),
      updating the model coefficients. Required signature:
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
      where _depvar_ is the dependent variable associated with the mix beta SDE,
      specified in the control file by the user, _ncomp_ is the number of
      components in the system, _moments_ is the map associating moment IDs
      (tk::ctr::vector< tk::ctr::Term >) to values of statistical moments,
      _bprime_, _kprime_ are user-defined parameters, and _b_, _k_ are the SDE
      parameters computed, see DiffEq/MixBeta.h.
*/
//******************************************************************************
#ifndef MixBetaCoeffPolicy_h
#define MixBetaCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! \brief Number-fraction beta SDE constant coefficients policity: b' and
//!    kappa' constants in time
class MixBetaCoeffConst {

    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor: initialize coefficients
    MixBetaCoeffConst(
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
              "Wrong number of mix beta SDE parameters 'b''");
      ErrChk( S_.size() == ncomp,
              "Wrong number of mix beta SDE parameters 'S'");
      ErrChk( kprime_.size() == ncomp,
              "Wrong number of mix beta SDE parameters 'kappa''");
      ErrChk( rho2_.size() == ncomp,
              "Wrong number of mix beta SDE parameters 'rho2'");
      ErrChk( rcomma_.size() == ncomp,
              "Wrong number of mix beta SDE parameters 'rcomma'");

      bprime = bprime_;
      S = S_;
      kprime = kprime_;
      rho2 = rho2_;
      rcomma = rcomma_;

      b.resize( bprime.size() );
      k.resize( kprime.size() );
    }

    //! Coefficients policy type accessor
    static tk::ctr::CoeffPolicyType type() noexcept
    { return tk::ctr::CoeffPolicyType::CONSTANT; }

    //! Update coefficients using constant coefficients for b' and kappa'
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

//! \brief Number-fraction beta SDE JRRJ coefficients policity
class MixBetaCoeffJRRJ {

  public:
    //! Constructor: initialize coefficients
    MixBetaCoeffJRRJ(
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
              "Wrong number of mix beta SDE parameters 'b''");
      ErrChk( S_.size() == ncomp,
              "Wrong number of mix beta SDE parameters 'S'");
      ErrChk( kprime_.size() == ncomp,
              "Wrong number of mix beta SDE parameters 'k''");
      ErrChk( rho2_.size() == ncomp,
              "Wrong number of mix beta SDE parameters 'rho2'");
      ErrChk( rcomma_.size() == ncomp,
              "Wrong number of mix beta SDE parameters 'rcomma'");

      bprime = bprime_;
      S = S_;
      kprime = kprime_;
      rho2 = rho2_;
      rcomma = rcomma_;

      b.resize( bprime.size() );
      k.resize( kprime.size() );
    }

    //! Coefficients policy type accessor
    static tk::ctr::CoeffPolicyType type() noexcept
    { return tk::ctr::CoeffPolicyType::JRRJ; }

    //! Update coefficients using the JRRJ model
    void update(
      char depvar,
      ncomp_t ncomp,
      const std::map< tk::ctr::Product, tk::real >& moments,
      const std::vector< kw::sde_bprime::info::expect::type  >& bprime,
      const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_kappa::info::expect::type >& k ) const {}
};

//! List of all beta's coefficients policies
using MixBetaCoeffPolicies =
  boost::mpl::vector< MixBetaCoeffConst
                    , MixBetaCoeffJRRJ
                    >;

} // walker::

#endif // MixBetaCoeffPolicy_h
