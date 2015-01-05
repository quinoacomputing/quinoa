//******************************************************************************
/*!
  \file      src/DiffEq/WFCoeffPolicy.h
  \author    J. Bakosi
  \date      Tue 13 Jan 2015 10:52:52 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Wright-Fisher coefficients policies
  \details   Wright-Fisher coefficients policies
*/
//******************************************************************************
#ifndef WFCoeffPolicy_h
#define WFCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! Wright-Fisher constant coefficients policity: constants in time
struct WFCoeffConst {

  //! Constructor: default for accessing policy name, type, etc.
  WFCoeffConst() = default;
  //! Constructor: initialize coefficients
  WFCoeffConst( tk::ctr::ncomp_type ncomp,
                const std::vector< tk::real >& omega_,
                std::vector< tk::real >& omega )
  {
    omega = omega_;
    ErrChk( omega.size() == ncomp,
            "Wrong number of Wright-Fisher SDE parameters 'omega'");
  }

  std::string policy() const noexcept
  { return tk::ctr::CoeffPolicy().name( tk::ctr::CoeffPolicyType::CONSTANT ); }

  tk::ctr::CoeffPolicyType type() const noexcept
  { return tk::ctr::CoeffPolicyType::CONSTANT; }

  //! Function call: no-op for constant coefficients
  void operator()( const tk::real&, std::vector< tk::real >& ) {}
};

//! List of all Wright-Fisher's coefficients policies
using WFCoeffPolicies = boost::mpl::vector< WFCoeffConst >;

} // walker::

#endif // WFCoeffPolicy_h
