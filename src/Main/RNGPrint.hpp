// *****************************************************************************
/*!
  \file      src/Main/RNGPrint.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Pretty printer base for pretty printers supporting RNGs
  \details   Pretty printer base for pretty printers supporting RNGs.
*/
// *****************************************************************************
#ifndef RNGPrint_h
#define RNGPrint_h

#include <iostream>
#include <vector>

#include "QuinoaConfig.hpp"
#include "Print.hpp"
#include "RNGParam.hpp"
#include "Options/RNG.hpp"

namespace tk {

//! RNGPrint : Print
class RNGPrint : public Print {

  public:
    //! Constructor
    //! \param[in] screen Screen output filename
    //! \param[in,out] str Verbose stream
    //! \param[in] mode Open mode for screen output file, see
    //!   http://en.cppreference.com/w/cpp/io/ios_base/openmode
    //! \param[in,out] qstr Quiet stream
    //! \see tk::Print::Print    
    explicit RNGPrint( const std::string& screen,
                       std::ostream& str = std::clog,
                       std::ios_base::openmode mode = std::ios_base::out,
                       std::ostream& qstr = std::cout )
      : Print( screen, str, mode, qstr ) {}

    #ifdef HAS_MKL
    //! Print all fields of MKL RNG parameters
    void MKLParams( const std::vector< ctr::RNGType >& vec,
                    const ctr::RNGMKLParameters& map ) const;
    #endif

    //! Print all fields of RNGSSE parameters
    void RNGSSEParams( const std::vector< ctr::RNGType >& vec,
                       const ctr::RNGSSEParameters& map ) const;

    //! Print all fields of Random123 parameters
    void Random123Params( const std::vector< ctr::RNGType >& vec,
                          const ctr::RNGRandom123Parameters& map ) const;

  private:
    #ifdef HAS_MKL
    //! Echo information on MKL random number generator
    void echoMKLParams( const ctr::RNGMKLParam& p ) const;
    #endif

    //! Echo information on RNGSSE random number generator
    void echoRNGSSEParams( const ctr::RNGSSEParam& p,
                           const ctr::RNG& rng,
                           const ctr::RNGType& r ) const;

    //! Echo information on Random123 random number generator
    void echoRandom123Params( const ctr::RNGRandom123Param& p ) const;
};

} // tk::

#endif // RNGPrint_h
