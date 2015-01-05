//******************************************************************************
/*!
  \file      src/DiffEq/SkewNormalCoeffPolicy.h
  \author    J. Bakosi
  \date      Tue 13 Jan 2015 10:53:08 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Skew-normal SDE coefficients policies
  \details   Skew-normal SDE coefficients policies
*/
//******************************************************************************
#ifndef SkewNormalCoeffPolicy_h
#define SkewNormalCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! Skew-normal SDE constant coefficients policity: constants in time
struct SkewNormalCoeffConst {

  //! Constructor: default for accessing policy name, type, etc.
  SkewNormalCoeffConst() = default;
  //! Constructor: initialize coefficients
  SkewNormalCoeffConst( tk::ctr::ncomp_type ncomp,
                        const std::vector< tk::real >& timescale_,
                        const std::vector< tk::real >& sigma_,
                        const std::vector< tk::real >& lambda_,
                        std::vector< tk::real >& timescale,
                        std::vector< tk::real >& sigma,
                        std::vector< tk::real >& lambda )
  {
    timescale = timescale_;
    sigma = sigma_;
    lambda = lambda_;
    ErrChk( timescale.size() == ncomp,
            "Wrong number of diagonal Skew-normal SDE parameters 'timescale'");
    ErrChk( sigma.size() == ncomp,
            "Wrong number of diagonal Skew-normal SDE parameters 'sigma'");
    ErrChk( lambda.size() == ncomp,
            "Wrong number of diagonal Skew-normal SDE parameters 'lambda'");
  }

  std::string policy() const noexcept
  { return tk::ctr::CoeffPolicy().name( tk::ctr::CoeffPolicyType::CONSTANT ); }

  tk::ctr::CoeffPolicyType type() const noexcept
  { return tk::ctr::CoeffPolicyType::CONSTANT; }

  //! Function call: no-op for constant coefficients
  void operator()( const tk::real&,
                   std::vector< tk::real >&,
                   std::vector< tk::real >&,
                   std::vector< tk::real >& ) {}
};

//! List of all Skew-normal SDE's coefficients policies
using SkewNormalCoeffPolicies = boost::mpl::vector< SkewNormalCoeffConst >;

} // walker::

#endif // SkewNormalCoeffPolicy_h
