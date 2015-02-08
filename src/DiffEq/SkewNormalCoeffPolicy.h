//******************************************************************************
/*!
  \file      src/DiffEq/SkewNormalCoeffPolicy.h
  \author    J. Bakosi
  \date      Sat 07 Feb 2015 07:30:37 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
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
class SkewNormalCoeffConst {

  public:
    //! Constructor: initialize coefficients
    SkewNormalCoeffConst( tk::ctr::ncomp_type ncomp,
                          const std::vector< tk::real >& timescale_,
                          const std::vector< tk::real >& sigmasq_,
                          const std::vector< tk::real >& lambda_,
                          std::vector< tk::real >& timescale,
                          std::vector< tk::real >& sigmasq,
                          std::vector< tk::real >& lambda )
    {
      ErrChk( timescale_.size() == ncomp,
        "Wrong number of diagonal Skew-normal SDE parameters 'timescale'");
      ErrChk( sigmasq_.size() == ncomp,
        "Wrong number of diagonal Skew-normal SDE parameters 'sigmasq'");
      ErrChk( lambda_.size() == ncomp,
        "Wrong number of diagonal Skew-normal SDE parameters 'lambda'");

      timescale = timescale_;
      sigmasq = sigmasq_;
      lambda = lambda_;
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

//! List of all Skew-normal SDE's coefficients policies
using SkewNormalCoeffPolicies = boost::mpl::vector< SkewNormalCoeffConst >;

} // walker::

#endif // SkewNormalCoeffPolicy_h
