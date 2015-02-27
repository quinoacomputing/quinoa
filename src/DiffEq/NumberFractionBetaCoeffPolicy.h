//******************************************************************************
/*!
  \file      src/DiffEq/NumberFractionBetaCoeffPolicy.h
  \author    J. Bakosi
  \date      Fri 27 Feb 2015 07:18:49 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Number-fraction beta SDE coefficients policies
  \details   This file defines coefficients policy classes for the
    number-fraction beta SDE, defined in DiffEq/NumberFractionBeta.h.

    General requirements on number-fraction beta SDE coefficients policy
    classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, b, S, kappa, rho2, and rcomma. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_type ncomp,
          const std::vector< kw::sde_b::info::expect::type >& b_,
          const std::vector< kw::sde_S::info::expect::type >& S_,
          const std::vector< kw::sde_kappa::info::expect::type >& k_,
          const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
          const std::vector< kw::sde_rcomma::info::expect::type >& rcomma_,
          std::vector< kw::sde_b::info::expect::type  >& b,
          std::vector< kw::sde_S::info::expect::type >& S,
          std::vector< kw::sde_kappa::info::expect::type >& k,
          std::vector< kw::sde_rho2::info::expect::type >& rho2_,
          std::vector< kw::sde_rcomma::info::expect::type >& rcomma_ );
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of
        number-fraction beta SDEs.
      - Constant references to b_, S_, k_, rho2_, and rcomma_, which denote five
        vectors of real values used to initialize the parameter vectors of the
        system of number-fraction beta SDEs. The length of the vectors must be
        equal to the number of components given by ncomp.
      - References to b, S, k, rho2_, and rcomma, which denote the parameter
        vectors to be initialized based on b_, S_, k_, rho2_, and rcomma_.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static tk::ctr::CoeffPolicyType type() noexcept {
          return tk::ctr::CoeffPolicyType::CONSTANT;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.

    - Must define the function _lookup()_, called from
      NumberFractionBeta::initialize(), performing pre-lookup of the locations of
      the statistical moments required by the given model. Required signature:
      \code{.cpp}
        void lookup( const tk::Statistics& stat,
                     char depvar,
                     tk::ctr::ncomp_type ncomp ) {}
      \endcode
      where _stat_ is the Statistics object, allowing access to the location of
      the various moments in memory, and _depvar_ is the dependent variable
      associated with the number-fraction beta SDE, given in the control file by
      the user, and _ncomp_ is the number of components.
*/
//******************************************************************************
#ifndef NumberFractionBetaCoeffPolicy_h
#define NumberFractionBetaCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! \brief Number-fraction beta SDE constant coefficients policity: constants in
//!   time
class NumberFractionBetaCoeffConst {

  public:
    //! Constructor: initialize coefficients
    NumberFractionBetaCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_b::info::expect::type >& b_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappa::info::expect::type >& k_,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
      const std::vector< kw::sde_rcomma::info::expect::type >& rcomma_,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappa::info::expect::type >& k,
      std::vector< kw::sde_rho2::info::expect::type >& rho2,
      std::vector< kw::sde_rcomma::info::expect::type >& rcomma )
    {
      ErrChk( b_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'b'");
      ErrChk( S_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'S'");
      ErrChk( k_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'kappa'");
      ErrChk( rho2_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'rho2'");
      ErrChk( rcomma_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'rcomma'");

      b = b_;
      S = S_;
      k = k_;
      rho2 = rho2_;
      rcomma = rcomma_;
    }

    //! Coefficients policy type accessor
    static tk::ctr::CoeffPolicyType type() noexcept
    { return tk::ctr::CoeffPolicyType::CONSTANT; }

    //! Lookup statistical moments required: no-op for constant coefficients
    void lookup( const tk::Statistics&, char, tk::ctr::ncomp_type ) {}

    //! Update coefficients: no-op for constant coefficients
    void update( std::vector< kw::sde_b::info::expect::type  >& b,
                 std::vector< kw::sde_S::info::expect::type >& S,
                 std::vector< kw::sde_kappa::info::expect::type >& k,
                 std::vector< kw::sde_rho2::info::expect::type >& rho2,
                 std::vector< kw::sde_rcomma::info::expect::type >& rcomma )
      const {}
};

//! \brief Number-fraction beta SDE JRRJ coefficients policity
class NumberFractionBetaCoeffJRRJ {

  public:
    //! Constructor: initialize coefficients
    NumberFractionBetaCoeffJRRJ(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_b::info::expect::type >& b_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappa::info::expect::type >& k_,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
      const std::vector< kw::sde_rcomma::info::expect::type >& rcomma_,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappa::info::expect::type >& k,
      std::vector< kw::sde_rho2::info::expect::type >& rho2,
      std::vector< kw::sde_rcomma::info::expect::type >& rcomma )
    {
      ErrChk( b_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'b'");
      ErrChk( S_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'S'");
      ErrChk( k_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'k'");
      ErrChk( rho2_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'rho2'");
      ErrChk( rcomma_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'rcomma'");

      b = b_;
      S = S_;
      k = k_;
      rho2 = rho2_;
      rcomma = rcomma_;
    }

    //! Coefficients policy type accessor
    static tk::ctr::CoeffPolicyType type() noexcept
    { return tk::ctr::CoeffPolicyType::JRRJ; }

    //! Lookup statistical moments required by the JRRJ model
    void lookup( const tk::Statistics& stat,
                 char depvar,
                 tk::ctr::ncomp_type ncomp )
    {
      // Generate tk::ctr::Terms for looking up the locations of required stats
      std::vector< tk::ctr::Product > means;
      std::vector< tk::ctr::Product > vars;
      for (tk::ctr::ncomp_type c=0; c<ncomp; ++c) {
        tk::ctr::Term m( toupper(depvar), c, tk::ctr::Moment::ORDINARY );
        means.push_back( tk::ctr::Product( { m } ) );
        tk::ctr::Term f( tolower(depvar), c, tk::ctr::Moment::CENTRAL );
        vars.push_back( tk::ctr::Product( { f, f } ) );
      }

      // Perform the lookups and store the locations 
      for (const auto& product : means) m_mean.push_back( stat.pos( product ) );
      for (const auto& product : vars) m_var.push_back( stat.pos( product ) );
    }

    //! Update coefficients using the JRRJ model
    void update( std::vector< kw::sde_b::info::expect::type  >& b,
                 std::vector< kw::sde_S::info::expect::type >& S,
                 std::vector< kw::sde_kappa::info::expect::type >& k,
                 std::vector< kw::sde_rho2::info::expect::type >& rho2,
                 std::vector< kw::sde_rcomma::info::expect::type >& rcomma )
      const
    {
      Assert( S.size() == m_mean.size() && S.size() == m_var.size(),
             "Vector lengths must match in NumberFractionBetaCoeffJRRJ::update()" );

      for (std::size_t i=0; i<S.size(); ++i) {
        const auto m = *m_mean[i];
        const auto v = *m_var[i];
        tk::real theta = 1.0 - v / (v + m * (1.0-m));
        S[i] = m + (1.0 - 2.0*m) * (1.0-theta) / (2.0-theta) *
                   (2.0/3.0 - 4.0/3.0*v*v) * (k[i]/b[i]*m*(1.0-m) - 1.0);
        std::cout << i << ": " << S[i] << std::endl;
        Assert( S[i] > 0.0 && S[i] < 1.0,
                "S out of bounds in NumberFractionBetaCoeffJRRJ::update()" );
      }
    }

  private:
    std::vector< const tk::real* > m_mean;      //!< Location of means
    std::vector< const tk::real* > m_var;       //!< Location of variances
};

//! List of all beta's coefficients policies
using NumberFractionBetaCoeffPolicies =
  boost::mpl::vector< NumberFractionBetaCoeffConst
                    , NumberFractionBetaCoeffJRRJ
                    >;

} // walker::

#endif // NumberFractionBetaCoeffPolicy_h
