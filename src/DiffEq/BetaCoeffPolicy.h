//******************************************************************************
/*!
  \file      src/DiffEq/BetaCoeffPolicy.h
  \author    J. Bakosi
  \date      Sat 25 Oct 2014 06:20:42 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Beta coefficients policies
  \details   Beta coefficients policies
*/
//******************************************************************************
#ifndef BetaCoeffPolicy_h
#define BetaCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Quinoa/Options/CoeffPolicy.h>

namespace quinoa {

//! Beta constant coefficients policity: constants in time
struct BetaCoeffConst {

  //! Constructor: default for accessing policy name, type, etc.
  BetaCoeffConst() = default;
  //! Constructor: initialize coefficients
  BetaCoeffConst( unsigned int ncomp,
                  const std::vector< tk::real >& b_,
                  const std::vector< tk::real >& S_,
                  const std::vector< tk::real >& k_,
                  std::vector< tk::real >& b,
                  std::vector< tk::real >& S,
                  std::vector< tk::real >& k )
  {
    b = b_;
    S = S_;
    k = k_;
    ErrChk( b.size() == ncomp, "Wrong number of beta SDE parameters 'b'");
    ErrChk( S.size() == ncomp, "Wrong number of beta SDE parameters 'S'");
    ErrChk( k.size() == ncomp, "Wrong number of beta SDE parameters 'k'");
  }

  std::string policy() const noexcept
  { return ctr::CoeffPolicy().name( ctr::CoeffPolicyType::CONSTANT ); }

  ctr::CoeffPolicyType type() const noexcept
  { return ctr::CoeffPolicyType::CONSTANT; }

  //! Function call: no-op for constant coefficients
  void operator()( const tk::real&,
                   std::vector< tk::real >&,
                   std::vector< tk::real >&,
                   std::vector< tk::real >& ) {}
};

//! List of all beta's coefficients policies
using BetaCoeffPolicies = boost::mpl::vector< BetaCoeffConst >;

} // quinoa::

#endif // BetaCoeffPolicy_h
