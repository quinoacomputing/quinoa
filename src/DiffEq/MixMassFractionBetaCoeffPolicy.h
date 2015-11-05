//******************************************************************************
/*!
  \file      src/DiffEq/MixMassFractionBetaCoeffPolicy.h
  \author    J. Bakosi
  \date      Thu 22 Oct 2015 02:11:07 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Mix mass-fraction beta SDE coefficients policies
  \details   This file defines coefficients policy classes for the mix
    mass-fraction beta SDE, defined in DiffEq/MixMassFractionBeta.h.

    General requirements on mix mass-fraction beta SDE coefficients policy
    classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, b, S, kappa, rho2, and r. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_type ncomp,
          const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
          const std::vector< kw::sde_S::info::expect::type >& S_,
          const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
          const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
          const std::vector< kw::sde_r::info::expect::type >& r_,
          std::vector< kw::sde_bprime::info::expect::type  >& bprime,
          std::vector< kw::sde_S::info::expect::type >& S,
          std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
          std::vector< kw::sde_rho2::info::expect::type >& rho2_,
          std::vector< kw::sde_r::info::expect::type >& r_ );
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of
        mix mass-fraction beta SDEs.
      - Constant references to bprime_, S_, kprime_, rho2_, and r_, which
        denote five vectors of real values used to initialize the parameter
        vectors of the system of mix mass-fraction beta SDEs. The length of
        the vectors must be equal to the number of components given by ncomp.
      - References to bprime, S, kprime, rho2_, and r, which denote the
        parameter vectors to be initialized based on bprime_, S_, kprime_,
        rho2_, and r_.

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
      MixMassFractionBeta::advance(), updating the model coefficients.
      Required signature:
      \code{.cpp}
        void update(
          char depvar,
          ncomp_t ncomp,
          const std::map< tk::ctr::Product, tk::real >& moments,
          const std::vector< kw::sde_bprime::info::expect::type  >& bprime,
          const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
          const std::vector< kw::sde_rho2::info::expect::type >& rho2,
          const std::vector< kw::sde_r::info::expect::type >& r,
          std::vector< kw::sde_b::info::expect::type  >& b,
          std::vector< kw::sde_kappa::info::expect::type >& k,
          std::vector< kw::sde_S::info::expect::type >& S ) const {}
      \endcode
      where _depvar_ is the dependent variable associated with the mix
      mass-fraction beta SDE, specified in the control file by the user, _ncomp_
      is the number of components in the system, _moments_ is the map
      associating moment IDs (tk::ctr::vector< tk::ctr::Term >) to values of
      statistical moments, _bprime_, _kprime_, rho2, r, are user-defined
      parameters, and _b_, _k_, _S_, are the SDE parameters computed, see
      DiffEq/MixMassFractionBeta.h.
*/
//******************************************************************************
#ifndef MixMassFractionBetaCoeffPolicy_h
#define MixMassFractionBetaCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Walker/Options/CoeffPolicy.h"

namespace walker {

//! \brief Mix mass-fraction beta SDE decay coefficients policity.
//! \details User-defined parameters b' and kappa' are constants in time and
//!   ensure decay in the evolution of <y^2>.
//! \author J. Bakosi
class MixMassFracBetaCoeffDecay {

    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor: initialize coefficients
    MixMassFracBetaCoeffDecay(
      ncomp_t ncomp,
      const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
      const std::vector< kw::sde_r::info::expect::type >& r_,
      std::vector< kw::sde_bprime::info::expect::type  >& bprime,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
      std::vector< kw::sde_rho2::info::expect::type >& rho2,
      std::vector< kw::sde_r::info::expect::type >& r,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_kappa::info::expect::type >& k )
    {
      ErrChk( bprime_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'b''");
      ErrChk( S_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'S'");
      ErrChk( kprime_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'kappa''");
      ErrChk( rho2_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'rho2'");
      ErrChk( r_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'r'");

      bprime = bprime_;
      S = S_;
      kprime = kprime_;
      rho2 = rho2_;
      r = r_;

      b.resize( bprime.size() );
      k.resize( kprime.size() );
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::DECAY; }

    //! \brief Update coefficients using constant coefficients for b' and kappa'
    //! \details This where the mix mass-fraction beta SDE is made consistent
    //!   with the no-mix and fully mixed limits by specifying the SDE
    //!   coefficients, b and kappa as functions of b' and kappa'. We leave S
    //!   unchanged.
    void update(
      char depvar,
      ncomp_t ncomp,
      const std::map< tk::ctr::Product, tk::real >& moments,
      const std::vector< kw::sde_bprime::info::expect::type  >& bprime,
      const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2,
      const std::vector< kw::sde_r::info::expect::type >& r,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_kappa::info::expect::type >& k,
      std::vector< kw::sde_S::info::expect::type >& S ) const
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

//! \brief Mix mass-fraction beta SDE homogneous decay coefficients policity.
//! \details User-defined parameters b' and kappa' are constants in time and
//!   ensure decay in the evolution of <y^2>. Additionally, S is constrained to
//!   make d<rho>/dt = 0, where <rho> = rho_2/(1+rY).
//! \author J. Bakosi
class MixMassFracBetaCoeffHomDecay {

    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor: initialize coefficients
    MixMassFracBetaCoeffHomDecay(
      ncomp_t ncomp,
      const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
      const std::vector< kw::sde_r::info::expect::type >& r_,
      std::vector< kw::sde_bprime::info::expect::type  >& bprime,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
      std::vector< kw::sde_rho2::info::expect::type >& rho2,
      std::vector< kw::sde_r::info::expect::type >& r,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_kappa::info::expect::type >& k )
    {
      ErrChk( bprime_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'b''");
      ErrChk( S_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'S'");
      ErrChk( kprime_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'kappa''");
      ErrChk( rho2_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'rho2'");
      ErrChk( r_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'r'");

      bprime = bprime_;
      S = S_;
      kprime = kprime_;
      rho2 = rho2_;
      r = r_;

      b.resize( bprime.size() );
      k.resize( kprime.size() );
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::HOMOGENEOUS_DECAY; }

    //! \brief Update coefficients b', kappa', and S
    //! \details This where the mix mass-fraction beta SDE is made consistent
    //!   with the no-mix and fully mixed limits by specifying the SDE
    //!   coefficients, b and kappa as functions of b' and kappa'. We also
    //!   specify S to force d<rho>/dt = 0, where <rho> = rho_2/(1+rY).
    void update(
      char depvar,
      ncomp_t ncomp,
      const std::map< tk::ctr::Product, tk::real >& moments,
      const std::vector< kw::sde_bprime::info::expect::type  >& bprime,
      const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2,
      const std::vector< kw::sde_r::info::expect::type >& r,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_kappa::info::expect::type >& k,
      std::vector< kw::sde_S::info::expect::type >& S ) const
    {
      using tk::ctr::lookup;
      using tk::ctr::mean;
      using tk::ctr::variance;
      using tk::ctr::cen3;
      // statistics nomenclature:
      //   Y = instantaneous mass fraction,
      //   R = instantaneous density,
      //   y = Y - <Y>, mass fraction fluctuation about its mean,
      //   r = R - <R>, density fluctuation about its mean,
      // <Y> = mean mass fraction,
      // <R> = mean density,
      //std::vector< tk::real > M{ 0.5, 0.012, 0.98, 0.37, 0.9 };
      for (ncomp_t c=0; c<ncomp; ++c) {
        tk::real m = lookup( mean(depvar,c), moments );            // <Y>
        tk::real v = lookup( variance(depvar,c), moments );        // <y^2>
        tk::real d = lookup( mean(depvar,c+ncomp), moments );      // <R>
        tk::real d2 = lookup( variance(depvar,c+ncomp), moments ); // <r^2>
        tk::real d3 = lookup( cen3(depvar,c+ncomp), moments );     // <r^3>

        if (m<1.0e-8 || m>1.0-1.0e-8) m = 0.5;
        if (v<1.0e-8 && v>1.0-1.0e-8) v = 0.5;
        b[c] = bprime[c] * (1.0 - v/m/(1.0-m));
        //b[c] = bprime[c] * (1.0 - v/M[c]/(1.0-M[c]));
        k[c] = kprime[c] * v;
        //b[c] = 1.0;
        //k[c] = 0.5*v/(m*(1.0-m));

        if (d < 1.0e-8) {
          std::cout << "d:" << d << " ";
          d = 0.5;
        }
        tk::real R = 1.0 + d2/d/d;
        tk::real B = -1.0/r[c]/r[c];
        tk::real C = (2.0+r[c])/r[c]/r[c];
        tk::real D = -(1.0+r[c])/r[c]/r[c];
        tk::real diff =
          B*d/rho2[c] +
          C*d*d*R/rho2[c]/rho2[c] +
          D*d*d*d*(1.0 + 3.0*d2/d/d + d3/d/d/d)/rho2[c]/rho2[c]/rho2[c];
        S[c] = (rho2[c]/d/R +
                2.0*k[c]/b[c]*rho2[c]*rho2[c]/d/d*r[c]*r[c]/R*diff - 1.0) / r[c];
        if (S[c] < 0.0 || S[c] > 1.0) {
          std::cout << S[c] << " ";
          S[c] = 0.5;
        }
      }
    }
};

//! \brief Mix mass-fraction beta SDE Monte Carlo homogneous decay coefficients
//!   policity.
//! \details User-defined parameters b' and kappa' are constants in time and
//!   ensure decay in the evolution of <y^2>. Additionally, S is constrained to
//!   make d<rho>/dt = 0, where <rho> = rho_2/(1+rY).
//! \author J. Bakosi
class MixMassFracBetaCoeffMonteCarloHomDecay {

    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor: initialize coefficients
    MixMassFracBetaCoeffMonteCarloHomDecay(
      ncomp_t ncomp,
      const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
      const std::vector< kw::sde_r::info::expect::type >& r_,
      std::vector< kw::sde_bprime::info::expect::type  >& bprime,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
      std::vector< kw::sde_rho2::info::expect::type >& rho2,
      std::vector< kw::sde_r::info::expect::type >& r,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_kappa::info::expect::type >& k )
    {
      ErrChk( bprime_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'b''");
      ErrChk( S_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'S'");
      ErrChk( kprime_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'kappa''");
      ErrChk( rho2_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'rho2'");
      ErrChk( r_.size() == ncomp,
        "Wrong number of mix mass-fraction beta SDE parameters 'r'");

      bprime = bprime_;
      S = S_;
      kprime = kprime_;
      rho2 = rho2_;
      r = r_;

      b.resize( bprime.size() );
      k.resize( kprime.size() );
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::MONTE_CARLO_HOMOGENEOUS_DECAY; }

    //! \brief Update coefficients b', kappa', and S
    //! \details This where the mix mass-fraction beta SDE is made consistent
    //!   with the no-mix and fully mixed limits by specifying the SDE
    //!   coefficients, b and kappa as functions of b' and kappa'. We also
    //!   specify S to force d<rho>/dt = 0, where <rho> = rho_2/(1+rY).
    void update(
      char depvar,
      ncomp_t ncomp,
      const std::map< tk::ctr::Product, tk::real >& moments,
      const std::vector< kw::sde_bprime::info::expect::type  >& bprime,
      const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2,
      const std::vector< kw::sde_r::info::expect::type >& r,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_kappa::info::expect::type >& k,
      std::vector< kw::sde_S::info::expect::type >& S ) const
    {
      using tk::ctr::lookup;
      using tk::ctr::mean;
      using tk::ctr::variance;
      using tk::ctr::ord2;
      // statistics nomenclature:
      //   Y = instantaneous mass fraction,
      //   R = instantaneous density,
      //   y = Y - <Y>, mass fraction fluctuation about its mean,
      //   r = R - <R>, density fluctuation about its mean,
      // <Y> = mean mass fraction,
      // <R> = mean density,
      for (ncomp_t c=0; c<ncomp; ++c) {
        tk::real m = lookup( mean(depvar,c), moments );            // <Y>
        tk::real v = lookup( variance(depvar,c), moments );        // <y^2>
        tk::real r2 = lookup( ord2(depvar,c+ncomp), moments );     // <R^2>

        const tk::ctr::Term Y( static_cast<char>(std::toupper(depvar)),
                               c,
                               tk::ctr::Moment::ORDINARY );
        const tk::ctr::Term R( static_cast<char>(std::toupper(depvar)),
                               c+ncomp,
                               tk::ctr::Moment::ORDINARY );
        const tk::ctr::Term OneMinusY( static_cast<char>(std::toupper(depvar)),
                                       c+3*ncomp,
                                       tk::ctr::Moment::ORDINARY );

        const auto YR2 = tk::ctr::Product( { Y, R, R } );
        const auto Y1MYR3 = tk::ctr::Product( { Y, OneMinusY, R, R, R } );

        tk::real yr2 = lookup( YR2, moments );          // <RY^2>
        tk::real y1myr3 = lookup( Y1MYR3, moments );    // <Y(1-Y)R^3>

        if (m<1.0e-8 || m>1.0-1.0e-8) m = 0.5;
        if (v<1.0e-8 || v>1.0-1.0e-8) v = 0.5;
        b[c] = bprime[c] * (1.0 - v/m/(1.0-m));
        k[c] = kprime[c] * v;
        //b[c] = 1.0;
        //k[c] = 0.5*v/(m*(1.0-m));

        if (r2 < 1.0e-8) {
          std::cout << "r2:" << r2 << " ";
          r2 = 0.5;
        }
        S[c] = (yr2 + 2.0*k[c]/b[c]*r[c]/rho2[c]*y1myr3) / r2;
        if (S[c] < 0.0 || S[c] > 1.0) {
          std::cout << "S:" << S[c] << " ";
          S[c] = 0.5;
        }
      }
    }
};

//! List of all mix mass-fraction beta's coefficients policies
using MixMassFracBetaCoeffPolicies =
  boost::mpl::vector< MixMassFracBetaCoeffDecay
                    , MixMassFracBetaCoeffHomDecay
                    , MixMassFracBetaCoeffMonteCarloHomDecay >;

} // walker::

#endif // MixMassFractionBetaCoeffPolicy_h
