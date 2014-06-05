//******************************************************************************
/*!
  \file      src/Base/RNGPrint.h
  \author    J. Bakosi
  \date      Sun 01 Jun 2014 12:42:34 PM MDT
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
    //! Constructor
    explicit RNGPrint() = default;

    //! Destructor
    ~RNGPrint() override = default;

    #ifdef HAS_MKL
    //! Print all fields of MKL RNG parameters
    void MKLParams( const std::vector< ctr::RNGType >& vec,
                    const ctr::RNGMKLParameters& map ) const;
    #endif

    //! Print all fields of RNGSSE parameters
    void RNGSSEParams( const std::vector< ctr::RNGType >& vec,
                       const ctr::RNGSSEParameters& map ) const;

  private:
    //! Don't permit copy constructor
    RNGPrint(const RNGPrint&) = delete;
    //! Don't permit copy assigment
    RNGPrint& operator=(const RNGPrint&) = delete;
    //! Don't permit move constructor
    RNGPrint(RNGPrint&&) = delete;
    //! Don't permit move assigment
    RNGPrint& operator=(RNGPrint&&) = delete;

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
