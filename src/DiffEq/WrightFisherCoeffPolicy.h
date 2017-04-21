// *****************************************************************************
/*!
  \file      src/DiffEq/WrightFisherCoeffPolicy.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
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
        static ctr::CoeffPolicyType type() noexcept {
          return ctr::CoeffPolicyType::CONSTANT;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef WrightFisherCoeffPolicy_h
#define WrightFisherCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Walker/Options/CoeffPolicy.h"

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
      ErrChk( omega_.size() == ncomp,
              "Wrong number of Wright-Fisher SDE parameters 'omega'");

      omega = omega_;
    }

    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONSTANT; }
};

//! List of all Wright-Fisher's coefficients policies
using WrightFisherCoeffPolicies = boost::mpl::vector< WrightFisherCoeffConst >;

} // walker::

#endif // WrightFisherCoeffPolicy_h
