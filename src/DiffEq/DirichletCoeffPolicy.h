//******************************************************************************
/*!
  \file      src/DiffEq/DirichletCoeffPolicy.h
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 11:31:10 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Dirichlet coefficients policies
  \details   This file defines coefficients policy classes for the Dirichlet
    SDE, defined in DiffEq/Dirichlet.h.

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
      - ncomp denotes the number of scalar components of the system of Dirichlet
        SDEs.
      - Constant references to b_, S_, and k_, which denote three vectors of
        real values used to initialize the parameter vectors of the Dirichlet
        SDEs. The length of the vectors must be equal to the number of
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

    - Must define the function _lookup()_, called from Dirichlet::initialize(),
      performing pre-lookup of the locations of the statistical moments required
      by the given model. Required signature:
      \code{.cpp}
        void lookup( const tk::Statistics& stat, char depvar )
      \endcode
      where _stat_ is the Statistics object, allowing access to the location of
      the various moments in memory, and _depvar_ is the dependent variable
      associated with the Dirichlet SDE, given in the control file by the user.

*/
//******************************************************************************
#ifndef DirichletCoeffPolicy_h
#define DirichletCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! Dirichlet constant coefficients policity: constants in time
class DirichletCoeffConst {

  public:
    //! Constructor: initialize coefficients
    DirichletCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_b::info::expect::type >& b_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappa::info::expect::type >& k_,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappa::info::expect::type >& k )
    {
      ErrChk( b_.size() == ncomp,
              "Wrong number of Dirichlet SDE parameters 'b'");
      ErrChk( S_.size() == ncomp,
              "Wrong number of Dirichlet SDE parameters 'S'");
      ErrChk( k_.size() == ncomp,
              "Wrong number of Dirichlet SDE parameters 'kappa'");

      b = b_;
      S = S_;
      k = k_;
    }

    static tk::ctr::CoeffPolicyType type() noexcept
    { return tk::ctr::CoeffPolicyType::CONSTANT; }

    //! Lookup statistical moments required: no-op for constant coefficients
    void lookup( const tk::Statistics&, char ) {}

    //! Function call: no-op for constant coefficients
    void operator()( const tk::real&,
                     std::vector< tk::real >&,
                     std::vector< tk::real >&,
                     std::vector< tk::real >& ) {}
};

//! List of all Dirichlet's coefficients policies
using DirichletCoeffPolicies = boost::mpl::vector< DirichletCoeffConst >;

} // walker::

#endif // DirichletCoeffPolicy_h
