//******************************************************************************
/*!
  \file      src/DiffEq/GenDirCoeffPolicy.h
  \author    J. Bakosi
  \date      Tue 13 Jan 2015 10:42:40 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Lochner's generalized Dirichlet coefficients policies
  \details   Lochner's generalized Dirichlet coefficients policies
*/
//******************************************************************************
#ifndef GenDirCoeffPolicy_h
#define GenDirCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! Generalized Dirichlet constant coefficients policity: constants in time
struct GenDirCoeffConst {

  //! Constructor: default for accessing policy name, type, etc.
  GenDirCoeffConst() = default;
  //! Constructor: initialize coefficients
  GenDirCoeffConst( tk::ctr::ncomp_type ncomp,
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
  { return tk::ctr::CoeffPolicy().name( tk::ctr::CoeffPolicyType::CONSTANT ); }

  tk::ctr::CoeffPolicyType type() const noexcept
  { return tk::ctr::CoeffPolicyType::CONSTANT; }

  //! Function call: no-op for constant coefficients
  void operator()( const tk::real&,
                   std::vector< tk::real >&,
                   std::vector< tk::real >&,
                   std::vector< tk::real >&,
                   std::vector< tk::real >& ) {}
};

//! List of all generalized Dirichlet's coefficients policies
using GenDirCoeffPolicies = boost::mpl::vector< GenDirCoeffConst >;

} // walker::

#endif // GenDirCoeffPolicy_h
