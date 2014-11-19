//******************************************************************************
/*!
  \file      src/DiffEq/WFCoeffPolicy.h
  \author    J. Bakosi
  \date      Mon 17 Nov 2014 07:48:58 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Wright-Fisher coefficients policies
  \details   Wright-Fisher coefficients policies
*/
//******************************************************************************
#ifndef WFCoeffPolicy_h
#define WFCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Quinoa/Options/CoeffPolicy.h>

namespace quinoa {

//! Wright-Fisher constant coefficients policity: constants in time
struct WFCoeffConst {

  //! Constructor: default for accessing policy name, type, etc.
  WFCoeffConst() = default;
  //! Constructor: initialize coefficients
  WFCoeffConst( unsigned int ncomp,
                const std::vector< tk::real >& omega_,
                std::vector< tk::real >& omega )
  {
    omega = omega_;
    ErrChk( omega.size() == ncomp,
            "Wrong number of Wright-Fisher SDE parameters 'omega'");
  }

  std::string policy() const noexcept
  { return ctr::CoeffPolicy().name( ctr::CoeffPolicyType::CONSTANT ); }

  ctr::CoeffPolicyType type() const noexcept
  { return ctr::CoeffPolicyType::CONSTANT; }

  //! Function call: no-op for constant coefficients
  void operator()( const tk::real&, std::vector< tk::real >& ) {}
};

//! List of all Wright-Fisher's coefficients policies
using WFCoeffPolicies = boost::mpl::vector< WFCoeffConst >;

} // quinoa::

#endif // WFCoeffPolicy_h
