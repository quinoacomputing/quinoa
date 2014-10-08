//******************************************************************************
/*!
  \file      src/SDE/GenDirCoeffPolicy.h
  \author    J. Bakosi
  \date      Wed 08 Oct 2014 10:27:53 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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

  //! Constructor: default for accessing policy name, type, etc.
  GenDirCoeffConst() = default;
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

  std::string policy() const noexcept
  { return ctr::CoeffPolicy().name( ctr::CoeffPolicyType::CONSTANT ); }

  ctr::CoeffPolicyType type() const noexcept
  { return ctr::CoeffPolicyType::CONSTANT; }

  //! Function call: no-op for constant coefficients
  void operator()( const tk::real&,
                   std::vector< tk::real >&,
                   std::vector< tk::real >&,
                   std::vector< tk::real >&,
                   std::vector< tk::real >& ) {}
};

//! List of all generalized Dirichlet's coefficients policies
using GenDirCoeffPolicies = boost::mpl::vector< GenDirCoeffConst >;

} // quinoa::

#endif // GenDirCoeffPolicy_h
