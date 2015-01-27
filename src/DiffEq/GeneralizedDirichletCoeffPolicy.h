//******************************************************************************
/*!
  \file      src/DiffEq/GeneralizedDirichletCoeffPolicy.h
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 11:40:48 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Lochner's generalized Dirichlet coefficients policies
  \details   This file defines coefficients policy classes for the generalized
    Dirichlet SDE, defined in DiffEq/GeneralizedDirichlet.h.

    General requirements on generalized Dirichlet SDE coefficients policy
    classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, b, S, kappa, and c. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_type ncomp,
          const std::vector< kw::sde_b::info::expect::type >& b_,
          const std::vector< kw::sde_S::info::expect::type >& S_,
          const std::vector< kw::sde_kappa::info::expect::type >& k_,
          const std::vector< kw::sde_c::info::expect::type >& c_,
          std::vector< kw::sde_b::info::expect::type  >& b,
          std::vector< kw::sde_S::info::expect::type >& S,
          std::vector< kw::sde_kappa::info::expect::type >& k,
          std::vector< kw::sde_c::info::expect::type >& c )
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of beta
        SDEs.
      - Constant references to b_, S_, k_, and c_, which denote four vectors of
        real values used to initialize the parameter vectors of the generalized
        Dirichlet SDEs. The length of the vectors b_, S_, and kappa_, must be
        equal to the number of components given by ncomp, while the length of
        vector c_ must be ncomp*(ncomp-1)/2.
      - References to b, S, k, and c, which denote the parameter vectors to be
        initialized based on b_, S_, k_, and c_.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static tk::ctr::CoeffPolicyType type() noexcept {
          return tk::ctr::CoeffPolicyType::CONSTANT;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.

    - Must define the function _lookup()_, called from GeneralizedDirichlet::
      initialize(), performing pre-lookup of the locations of the statistical
      moments required by the given model. Required signature:
      \code{.cpp}
        void lookup( const tk::Statistics& stat, char depvar ) {}
      \endcode
      where _stat_ is the Statistics object, allowing access to the location of
      the various moments in memory, and _depvar_ is the dependent variable
      associated with the generalized Dirichlet SDE, given in the control file
      by the user.
*/
//******************************************************************************
#ifndef GeneralizedDirichletCoeffPolicy_h
#define GeneralizedDirichletCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! Generalized Dirichlet constant coefficients policity: constants in time
class GeneralizedDirichletCoeffConst {

  public:
    //! Constructor: initialize coefficients
    GeneralizedDirichletCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_b::info::expect::type >& b_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappa::info::expect::type >& k_,
      const std::vector< kw::sde_c::info::expect::type >& c_,
      std::vector< kw::sde_b::info::expect::type >& b,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappa::info::expect::type >& k,
      std::vector< kw::sde_c::info::expect::type >& c )
    {
      b = b_;
      S = S_;
      k = k_;
      c = c_;
      ErrChk( b.size() == ncomp,
              "Wrong number of generalized Dirichlet SDE parameters 'b'");
      ErrChk( S.size() == ncomp,
              "Wrong number of generalized Dirichlet SDE parameters 'S'");
      ErrChk( k.size() == ncomp,
              "Wrong number of generalized Dirichlet SDE parameters 'kappa'");
      ErrChk( c.size() == ncomp*(ncomp-1)/2,
              "Wrong number of generalized Dirichlet SDE parameters 'c'");
    }

    static tk::ctr::CoeffPolicyType type() noexcept
    { return tk::ctr::CoeffPolicyType::CONSTANT; }

    //! Lookup statistical moments required: no-op for constant coefficients
    void lookup( const tk::Statistics& stat, char ) {}

    //! Function call: no-op for constant coefficients
    void operator()( const tk::real&,
                     std::vector< tk::real >&,
                     std::vector< tk::real >&,
                     std::vector< tk::real >&,
                     std::vector< tk::real >& ) {}
};

//! List of all generalized Dirichlet's coefficients policies
using GeneralizedDirichletCoeffPolicies =
  boost::mpl::vector< GeneralizedDirichletCoeffConst >;

} // walker::

#endif // GeneralizedDirichletCoeffPolicy_h
