//******************************************************************************
/*!
  \file      src/RNG/RNGStack.h
  \author    J. Bakosi
  \date      Tue 24 Jun 2014 07:47:20 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Stack of random number generators
  \details   Stack of random number generators
*/
//******************************************************************************
#ifndef RNGStack_h
#define RNGStack_h

#include <map>

#include <tkTypes.h>
#include <RNG.h>
#include <Options/RNG.h>

namespace tk {

//! Random number generator factory type
using RNGFactory = std::map< ctr::RNGType, std::function< RNG() > >;

//! RNGStack
class RNGStack {

  public:
    //! Register random number generators into factory
    void initFactory( RNGFactory& factory,
                      int nthreads,
                      #ifdef HAS_MKL
                      const ctr::RNGMKLParameters& mklparam,
                      #endif
                      const ctr::RNGSSEParameters& rngsseparam );

    //! Instantiate selected RNGs
    std::map< std::underlying_type<tk::ctr::RNGType>::type, tk::RNG >
    createSelected( const tk::RNGFactory& factory,
                    const std::vector< ctr::RNGType >& selected );

  private:
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

#endif // RNGStack_h
