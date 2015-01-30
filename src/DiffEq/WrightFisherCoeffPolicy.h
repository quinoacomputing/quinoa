//******************************************************************************
/*!
  \file      src/DiffEq/WrightFisherCoeffPolicy.h
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 12:02:06 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Wright-Fisher coefficients policies
  \details   This file defines coefficients policy classes for the Wright-Fisher
    SDE, defined in DiffEq/WrightFisher.h.

    General requirements on the Wright-Fisher SDE coefficients policy classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficient vector, omega. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_type ncomp,
          const std::vector< kw::sde_omega::info::expect::type >& omega_,
          std::vector< kw::sde_omega::info::expect::type >& omega )
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of the
        Wright-Fisher SDE.
      - Constant references to omega_ which denote a vector of real values used
        to initialize the parameter vector of the Wright-Fisher SDE. The length
        of the vector must be equal to the number of components given by ncomp.
      - Reference to omega which denote the parameter vector to be initialized
        based on omega_.

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
      WrightFisher::initialize(), performing pre-lookup of the locations of the
      statistical moments required by the given model. Required signature:
      \code{.cpp}
        void lookup( const tk::Statistics& stat, char depvar )
      \endcode
      where _stat_ is the Statistics object, allowing access to the location of
      the various moments in memory, and _depvar_ is the dependent variable
      associated with the Wright-Fisher SDE, given in the control file by the
      user.
*/
//******************************************************************************
#ifndef WrightFisherCoeffPolicy_h
#define WrightFisherCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! Wright-Fisher constant coefficients policity: constants in time
class WrightFisherCoeffConst {

  public:
    //! Constructor: initialize coefficients
    WrightFisherCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_omega::info::expect::type >& omega_,
      std::vector< kw::sde_omega::info::expect::type >& omega )
    {
      omega = omega_;
      ErrChk( omega.size() == ncomp,
              "Wrong number of Wright-Fisher SDE parameters 'omega'");
    }

    static tk::ctr::CoeffPolicyType type() noexcept
    { return tk::ctr::CoeffPolicyType::CONSTANT; }

    //! Lookup statistical moments required: no-op for constant coefficients
    void lookup( const tk::Statistics& stat, char ) {}

    //! Function call: no-op for constant coefficients
    void operator()( const tk::real&, std::vector< tk::real >& ) {}
};

//! List of all Wright-Fisher's coefficients policies
using WrightFisherCoeffPolicies = boost::mpl::vector< WrightFisherCoeffConst >;

} // walker::

#endif // WrightFisherCoeffPolicy_h
