//******************************************************************************
/*!
  \file      src/SDE/DirCoeffPolicy.h
  \author    J. Bakosi
  \date      Wed 08 Oct 2014 10:27:42 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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

  //! Constructor: default for accessing policy name, type, etc.
  DirCoeffConst() = default;
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
    ErrChk( b.size() == ncomp, "Wrong number of Dirichlet SDE parameters 'b'");
    ErrChk( S.size() == ncomp, "Wrong number of Dirichlet SDE parameters 'S'");
    ErrChk( k.size() == ncomp, "Wrong number of Dirichlet SDE parameters 'k'");
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

//! List of all Dirichlet's coefficients policies
using DirCoeffPolicies = boost::mpl::vector< DirCoeffConst >;

} // quinoa::

#endif // DirCoeffPolicy_h
