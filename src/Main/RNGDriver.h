//******************************************************************************
/*!
  \file      src/Main/RNGDriver.h
  \author    J. Bakosi
  \date      Tue 03 Jun 2014 10:15:10 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver with RNGs
  \details   Driver with RNGs
*/
//******************************************************************************
#ifndef RNGDriver_h
#define RNGDriver_h

#include <map>

#include <Driver.h>
#include <tkTypes.h>
#include <RNG.h>
#include <Options/RNG.h>

namespace tk {

//! Random number generator factory type
using RNGFactory = std::map< ctr::RNGType, std::function< RNG() > >;

//! RNGDriver
class RNGDriver : public Driver {

  public:
    //! Constructor
    explicit RNGDriver() = default;

    //! Register random number generators into factory
    void initFactory( RNGFactory& factory,
                      int nthreads,
                      #ifdef HAS_MKL
                      const ctr::RNGMKLParameters& mklparam,
                      #endif
                      const ctr::RNGSSEParameters& rngsseparam );

    //! Instantiate selected RNGs
    std::map< tk::ctr::RNGType, tk::RNG >
    createSelected( const tk::RNGFactory& factory,
                    const std::vector< ctr::RNGType >& selected );

  private:
    //! Don't permit copy constructor
    RNGDriver(const RNGDriver&) = delete;
    //! Don't permit assigment constructor
    RNGDriver& operator=(const RNGDriver&) = delete;
    //! Don't permit move constructor
    RNGDriver(RNGDriver&&) = delete;
    //! Don't permit move assignment
    RNGDriver& operator=(RNGDriver&&) = delete;

   #ifdef HAS_MKL
   //! Register MKL RNGs into factory
   void regMKL( RNGFactory& factory,
                int nthreads,
                const ctr::RNGMKLParameters& mklparam );
   #endif

   //! Register RNGSSE RNGs into factory
   void regRNGSSE( RNGFactory& factory,
                   int nthreads,
                   const ctr::RNGSSEParameters& param );
};

} // namespace tk

#endif // RNGDriver_h
