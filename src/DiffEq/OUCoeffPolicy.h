//******************************************************************************
/*!
  \file      src/SDE/OUCoeffPolicy.h
  \author    J. Bakosi
  \date      Wed 08 Oct 2014 10:29:51 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Ornstein-Uhlenbeck coefficients policies
  \details   Ornstein-Uhlenbeck coefficients policies
*/
//******************************************************************************
#ifndef OUCoeffPolicy_h
#define OUCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Quinoa/Options/CoeffPolicy.h>

namespace quinoa {

//! Ornstein-Uhlenbeck constant coefficients policity: constants in time
struct OUCoeffConst {

  //! Constructor: default for accessing policy name, type, etc.
  OUCoeffConst() = default;
  //! Constructor: initialize coefficients
  OUCoeffConst( unsigned int ncomp,
                 const std::vector< tk::real >& sigma_,
                 const std::vector< tk::real >& timescale_,
                 std::vector< tk::real >& sigma,
                 std::vector< tk::real >& timescale )
  {
    sigma = sigma_;
    timescale = timescale_;
    ErrChk( sigma.size() == ncomp, "Wrong number of OU SDE parameters 'sigma'");
    ErrChk( timescale.size() == ncomp,
            "Wrong number of OU SDE parameters 'timscale'");
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

//! List of all Ornstein-Uhlenbeck's coefficients policies
using OUCoeffPolicies = boost::mpl::vector< OUCoeffConst >;

} // quinoa::

#endif // OUCoeffPolicy_h
