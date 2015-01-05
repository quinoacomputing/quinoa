//******************************************************************************
/*!
  \file      src/DiffEq/GammaCoeffPolicy.h
  \author    J. Bakosi
  \date      Tue 13 Jan 2015 10:49:04 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Gamma coefficients policies
  \details   Gamma coefficients policies
*/
//******************************************************************************
#ifndef GammaCoeffPolicy_h
#define GammaCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! Gamma constant coefficients policity: constants in time
struct GammaCoeffConst {

  //! Constructor: default for accessing policy name, type, etc.
  GammaCoeffConst() = default;
  //! Constructor: initialize coefficients
  GammaCoeffConst( tk::ctr::ncomp_type ncomp,
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
    ErrChk( b.size() == ncomp, "Wrong number of gamma SDE parameters 'b'");
    ErrChk( S.size() == ncomp, "Wrong number of gamma SDE parameters 'S'");
    ErrChk( k.size() == ncomp, "Wrong number of gamma SDE parameters 'k'");
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

//! List of all gamma's coefficients policies
using GammaCoeffPolicies = boost::mpl::vector< GammaCoeffConst >;

} // walker::

#endif // GammaCoeffPolicy_h
