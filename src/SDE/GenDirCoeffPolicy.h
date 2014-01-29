//******************************************************************************
/*!
  \file      src/SDE/GenDirCoeffPolicy.h
  \author    J. Bakosi
  \date      Wed Jan 29 14:59:48 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Lochner's generalized Dirichlet coefficients policies
  \details   Lochner's generalized Dirichlet coefficients policies
*/
//******************************************************************************
#ifndef GenDirCoeffPolicy_h
#define GenDirCoeffPolicy_h

namespace quinoa {

//! Generalized Dirichlet constant coefficients policity: constants in time
struct GenDirCoeffConst {

  //! Constructor: initialize coefficients
  GenDirCoeffConst( std::string& policy,
                    unsigned int ncomp,
                    const std::vector< tk::real >& b_,
                    const std::vector< tk::real >& S_,
                    const std::vector< tk::real >& k_,
                    const std::vector< tk::real >& c_,
                    std::vector< tk::real >& b,
                    std::vector< tk::real >& S,
                    std::vector< tk::real >& k,
                    std::vector< tk::real >& c )
  {
    policy = "constant";
    b = b_; S = S_; k = k_; c = c_;
    ErrChk( b.size() == ncomp, tk::ExceptType::FATAL,
            "Wrong number of generalized Dirichlet model parameters 'b'");
    ErrChk( S.size() == ncomp, tk::ExceptType::FATAL,
            "Wrong number of generalized Dirichlet model parameters 'S'");
    ErrChk( k.size() == ncomp, tk::ExceptType::FATAL,
            "Wrong number of generalized Dirichlet model parameters 'k'");
    ErrChk( c.size() == ncomp*(ncomp-1)/2, tk::ExceptType::FATAL,
            "Wrong number of generalized Dirichlet model parameters 'c'");
  }

  //! Function call: no-op for constant coefficients
  void operator()( const tk::real& t,
                   std::vector< tk::real >& b,
                   std::vector< tk::real >& S,
                   std::vector< tk::real >& k,
                   std::vector< tk::real >& c ) {}

};

} // quinoa::

#endif // GenDirCoeffPolicy_h
