//******************************************************************************
/*!
  \file      src/DiffEq/GammaCoeffPolicy.h
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 11:34:00 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Gamma coefficients policies
  \details   This file defines coefficients policy classes for the gamma SDE,
    defined in DiffEq/Gamma.h.

    General requirements on gamma SDE coefficients policy classes:

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
      - ncomp denotes the number of scalar components of the system of gamma
        SDEs.
      - Constant references to b_, S_, and k_, which denote three vectors of
        real values used to initialize the parameter vectors of the system of
        gamma SDEs. The length of the vectors must be equal to the number of
        components given by ncomp.
      - References to b, S, and k, which denote the parameter vectors to be
        initialized based on b_, S_, and k_.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static tk::ctr::CoeffPolicyType type() noexcept {
          return tk::ctr::CoeffPolicyType::CONSTANT;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.

    - Must define the function _lookup()_, called from Gamma::initialize(),
      performing pre-lookup of the locations of the statistical moments required
      by the given model. Required signature:
      \code{.cpp}
        void lookup( const tk::Statistics& stat, char depvar ) {}
      \endcode
      where _stat_ is the Statistics object, allowing access to the location of
      the various moments in memory, and _depvar_ is the dependent variable
      associated with the gamma SDE, given in the control file by the user.
*/
//******************************************************************************
#ifndef GammaCoeffPolicy_h
#define GammaCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! Gamma constant coefficients policity: constants in time
class GammaCoeffConst {

  public:
    //! Constructor: initialize coefficients
    GammaCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_b::info::expect::type >& b_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappa::info::expect::type >& k_,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappa::info::expect::type >& k )
    {
      b = b_;
      S = S_;
      k = k_;
      ErrChk( b.size() == ncomp, "Wrong number of gamma SDE parameters 'b'");
      ErrChk( S.size() == ncomp, "Wrong number of gamma SDE parameters 'S'");
      ErrChk( k.size() == ncomp, "Wrong number of gamma SDE parameters 'kappa'");
    }

    static tk::ctr::CoeffPolicyType type() noexcept
    { return tk::ctr::CoeffPolicyType::CONSTANT; }

    //! Lookup statistical moments required: no-op for constant coefficients
    void lookup( const tk::Statistics& stat, char ) {}

    //! Function call: no-op for constant coefficients
    void operator()( const tk::real&,
                     std::vector< tk::real >&,
                     std::vector< tk::real >&,
                     std::vector< tk::real >& ) {}
};

//! List of all gamma's coefficients policies
using GammaCoeffPolicies = boost::mpl::vector< GammaCoeffConst >;

} // walker::

#endif // GammaCoeffPolicy_h
