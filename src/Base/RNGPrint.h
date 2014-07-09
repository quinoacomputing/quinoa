//******************************************************************************
/*!
  \file      src/Base/RNGPrint.h
  \author    J. Bakosi
  \date      Sun 06 Jul 2014 08:31:46 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Printer with RNGs
  \details   Printer with RNGs
*/
//******************************************************************************
#ifndef RNGPrint_h
#define RNGPrint_h

#include <Print.h>
#include <tkTypes.h>
#include <Options/RNG.h>

namespace tk {

//! RNGPrint : Print
class RNGPrint : public Print {

  public:
    #ifdef HAS_MKL
    //! Print all fields of MKL RNG parameters
    void MKLParams( const std::vector< ctr::RNGType >& vec,
                    const ctr::RNGMKLParameters& map ) const;
    #endif

    //! Print all fields of RNGSSE parameters
    void RNGSSEParams( const std::vector< ctr::RNGType >& vec,
                       const ctr::RNGSSEParameters& map ) const;

  private:
    #ifdef HAS_MKL
    //! Echo information on MKL random number generator
    void echoMKLParams( const ctr::RNGMKLParam& p ) const;
    #endif

    //! Echo information on RNGSSE random number generator
    void echoRNGSSEParams( const ctr::RNGSSEParam& p,
                           const ctr::RNG& rng,
                           const ctr::RNGType& r ) const;
};

} // tk::

#endif // RNGPrint_h
