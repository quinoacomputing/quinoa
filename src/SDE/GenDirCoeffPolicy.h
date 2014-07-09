//******************************************************************************
/*!
  \file      src/SDE/GenDirCoeffPolicy.h
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:08:50 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Lochner's generalized Dirichlet coefficients policies
  \details   Lochner's generalized Dirichlet coefficients policies
*/
//******************************************************************************
#ifndef GenDirCoeffPolicy_h
#define GenDirCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Quinoa/Options/CoeffPolicy.h>

namespace quinoa {

//! Generalized Dirichlet constant coefficients policity: constants in time
struct GenDirCoeffConst {

  GenDirCoeffConst() = default;

  std::string policy() const noexcept {
    return ctr::CoeffPolicy().name( ctr::CoeffPolicyType::CONSTANT );
  }

  ctr::CoeffPolicyType type() const noexcept {
    return ctr::CoeffPolicyType::CONSTANT;
  }

  //! Constructor: initialize coefficients
  GenDirCoeffConst( unsigned int ncomp,
                    const std::vector< tk::real >& b_,
                    const std::vector< tk::real >& S_,
                    const std::vector< tk::real >& k_,
                    const std::vector< tk::real >& c_,
                    std::vector< tk::real >& b,
                    std::vector< tk::real >& S,
                    std::vector< tk::real >& k,
                    std::vector< tk::real >& c )
  {
    b = b_;
    S = S_;
    k = k_;
    c = c_;
    ErrChk( b.size() == ncomp,
            "Wrong number of generalized Dirichlet SDE parameters 'b'");
    ErrChk( S.size() == ncomp,
            "Wrong number of generalized Dirichlet SDE parameters 'S'");
    ErrChk( k.size() == ncomp,
            "Wrong number of generalized Dirichlet SDE parameters 'k'");
    ErrChk( c.size() == ncomp*(ncomp-1)/2,
            "Wrong number of generalized Dirichlet SDE parameters 'c'");
  }

  //! Function call: no-op for constant coefficients
  void operator()( const tk::real& t,
                   std::vector< tk::real >& b,
                   std::vector< tk::real >& S,
                   std::vector< tk::real >& k,
                   std::vector< tk::real >& c )
  {
    IGNORE(t);
    IGNORE(b);
    IGNORE(S);
    IGNORE(k);
    IGNORE(c);
  }

};

//! List of all generalized Dirichlet's coefficients policies
using GenDirCoeffPolicies = boost::mpl::vector< GenDirCoeffConst >;

} // quinoa::

#endif // GenDirCoeffPolicy_h
