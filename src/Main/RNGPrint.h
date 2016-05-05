//******************************************************************************
/*!
  \file      src/Main/RNGPrint.h
  \author    J. Bakosi
  \date      Mon 02 May 2016 07:07:14 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Pretty printer base for pretty printers supporting RNGs
  \details   Pretty printer base for pretty printers supporting RNGs.
*/
//******************************************************************************
#ifndef RNGPrint_h
#define RNGPrint_h

#include <iostream>
#include <vector>

#include "Print.h"
#include "RNGParam.h"
#include "Options/RNG.h"

namespace tk {

//! RNGPrint : Print
class RNGPrint : public Print {

  public:
    //! Constructor
    //! \param[in,out] str Verbose stream
    //! \param[in,out] qstr Quiet stream
    //! \see tk::Print::Print    
    explicit RNGPrint( std::ostream& str = std::clog,
                       std::ostream& qstr = std::cout ) : Print( str, qstr ) {}

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
