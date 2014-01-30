//******************************************************************************
/*!
  \file      src/SDE/DirCoeffPolicy.h
  \author    J. Bakosi
  \date      Wed Jan 29 16:24:34 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet coefficients policies
  \details   Dirichlet coefficients policies
*/
//******************************************************************************
#ifndef DirCoeffPolicy_h
#define DirCoeffPolicy_h

namespace quinoa {

//! Dirichlet constant coefficients policity: constants in time
struct DirCoeffConst {

  //! Constructor: initialize coefficients
  DirCoeffConst( std::string& policy,
                 unsigned int ncomp,
                 const std::vector< tk::real >& b_,
                 const std::vector< tk::real >& S_,
                 const std::vector< tk::real >& k_,
                 std::vector< tk::real >& b,
                 std::vector< tk::real >& S,
                 std::vector< tk::real >& k )
  {
    policy = "constant";
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
                   std::vector< tk::real >& k ) {}

};

} // quinoa::

#endif // DirCoeffPolicy_h
