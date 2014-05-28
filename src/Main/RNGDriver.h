//******************************************************************************
/*!
  \file      src/Main/RNGDriver.h
  \author    J. Bakosi
  \date      Tue 27 May 2014 10:50:04 AM MDT
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
                      const ctr::MKLRNGParameters& mklparam,
                      #endif
                      const ctr::RNGSSEParameters& rngsseparam );

    //! Instantiate all registered RNGs
    std::vector< tk::RNG >
      instantiateAll( const tk::RNGFactory& factory );

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
                const ctr::MKLRNGParameters& mklparam );
   #endif

   //! Register RNGSSE RNGs into factory
   void regRNGSSE( RNGFactory& factory,
                   int nthreads,
                   const ctr::RNGSSEParameters& param );
};

} // namespace tk

#endif // RNGDriver_h
