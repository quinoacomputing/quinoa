//******************************************************************************
/*!
  \file      src/SDE/SkewNormalCoeffPolicy.h
  \author    J. Bakosi
  \date      Wed Jan 29 16:25:29 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Skew-normal coefficients policies
  \details   Skew-normal coefficients policies
*/
//******************************************************************************
#ifndef SkewNormalCoeffPolicy_h
#define SkewNormalCoeffPolicy_h

namespace quinoa {

//! Skew-normal constant coefficients policity: constants in time
struct SkewNormalCoeffConst {

  //! Constructor: initialize coefficients
  SkewNormalCoeffConst( std::string& policy,
                        const tk::real& sigma_,
                        const tk::real& timescale_,
                        const tk::real& lambda_,
                        tk::real& sigma,
                        tk::real& timescale,
                        tk::real& lambda )
  {
    policy = "constant";
    sigma = sigma_;
    timescale = timescale_;
    lambda = lambda_;
  }

  //! Function call: no-op for constant coefficients
  void operator()( const tk::real& t,
                   tk::real& sigma,
                   tk::real& timescale,
                   tk::real& lambda ) {}

};

} // quinoa::

#endif // SkewNormalCoeffPolicy_h
