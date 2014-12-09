//******************************************************************************
/*!
  \file      src/DiffEq/OUCoeffPolicy.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 06:32:43 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Ornstein-Uhlenbeck coefficients policies
  \details   Ornstein-Uhlenbeck coefficients policies
*/
//******************************************************************************
#ifndef OUCoeffPolicy_h
#define OUCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! Ornstein-Uhlenbeck constant coefficients policity: constants in time
struct OUCoeffConst {

  //! Constructor: default for accessing policy name, type, etc.
  OUCoeffConst() = default;
  //! Constructor: initialize coefficients
  OUCoeffConst( unsigned int ncomp,
                const std::vector< tk::real >& sigma_,
                const std::vector< tk::real >& theta_,
                const std::vector< tk::real >& mu_,
                std::vector< tk::real >& sigma,
                std::vector< tk::real >& theta,
                std::vector< tk::real >& mu )
  {
    ErrChk( sigma_.size() == ncomp*(ncomp+1)/2,
            "Wrong number of OU SDE parameters 'sigma'");
    sigma.resize( ncomp * ncomp );
    std::size_t c = 0;
    for (std::size_t i=0; i<ncomp; ++i)
      for (std::size_t j=0; j<ncomp; ++j)
        if (i<=j) sigma[ i*ncomp+j ] = sigma_[ c++ ];

    theta = theta_;
    mu = mu_;
    ErrChk( theta.size() == ncomp, "Wrong number of OU SDE parameters 'theta'");
    ErrChk( mu.size() == ncomp, "Wrong number of OU SDE parameters 'mu'");
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

//! List of all Ornstein-Uhlenbeck's coefficients policies
using OUCoeffPolicies = boost::mpl::vector< OUCoeffConst >;

} // walker::

#endif // OUCoeffPolicy_h
