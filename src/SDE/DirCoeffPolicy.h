//******************************************************************************
/*!
  \file      src/SDE/DirCoeffPolicy.h
  \author    J. Bakosi
  \date      Wed Mar 19 15:43:10 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet coefficients policies
  \details   Dirichlet coefficients policies
*/
//******************************************************************************
#ifndef DirCoeffPolicy_h
#define DirCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Quinoa/Options/CoeffPolicy.h>

namespace quinoa {

//! Dirichlet constant coefficients policity: constants in time
struct DirCoeffConst {

  DirCoeffConst() = default;

  std::string policy() const noexcept {
    return ctr::CoeffPolicy().name( ctr::CoeffPolicyType::CONSTANT );
  }

  ctr::CoeffPolicyType type() const noexcept {
    return ctr::CoeffPolicyType::CONSTANT;
  }

  //! Constructor: initialize coefficients
  DirCoeffConst( unsigned int ncomp,
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
    ErrChk( b.size() == ncomp, tk::ExceptType::FATAL,
            "Wrong number of Dirichlet model parameters 'b'");
    ErrChk( S.size() == ncomp, tk::ExceptType::FATAL,
             "Wrong number of Dirichlet model parameters 'S'");
    ErrChk( k.size() == ncomp, tk::ExceptType::FATAL,
             "Wrong number of Dirichlet model parameters 'k'");
  }

  //! Function call: no-op for constant coefficients
  void operator()( const tk::real& t,
                   std::vector< tk::real >& b,
                   std::vector< tk::real >& S,
                   std::vector< tk::real >& k )
  {
    IGNORE(t);
    IGNORE(b);
    IGNORE(S);
    IGNORE(k);
  }
};

//! List of all Dirichlet's coefficients policies
using DirCoeffPolicies = boost::mpl::vector< DirCoeffConst >;

} // quinoa::

#endif // DirCoeffPolicy_h
