//******************************************************************************
/*!
  \file      src/SDE/DirCoeffPolicy.h
  \author    J. Bakosi
  \date      Tue 14 Jan 2014 09:27:45 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet coefficients policies
  \details   Dirichlet coefficients policies
*/
//******************************************************************************
#ifndef DirCoeffPolicy_h
#define DirCoeffPolicy_h

namespace quinoa {

//! Dirichlet constant coefficients policity: constants
struct DirCoeffConst {

  //! Constructor: initialize coefficients
  DirCoeffConst( unsigned int ncomp,
                 const std::vector< tk::real >& b_,
                 const std::vector< tk::real >& S_,
                 const std::vector< tk::real >& k_,
                 std::vector< tk::real >& b,
                 std::vector< tk::real >& S,
                 std::vector< tk::real >& k )
  {
    b = b_; S = S_; k = k_;
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
