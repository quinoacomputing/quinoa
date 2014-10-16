//******************************************************************************
/*!
  \file      src/RNG/RNGStack.h
  \author    J. Bakosi
  \date      Tue 12 Aug 2014 05:06:40 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
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
    //! Constructor: register random number generators into factory
    explicit RNGStack(
                       #ifdef HAS_MKL
                       const ctr::RNGMKLParameters& mklparam,
                       #endif
                       const ctr::RNGSSEParameters& rngsseparam );

    //! Instantiate selected RNGs
    std::map< std::underlying_type< tk::ctr::RNGType >::type, tk::RNG >
    selected( const std::vector< ctr::RNGType >& sel ) const;

    //! Instantiate a RNG
    tk::RNG create( tk::ctr::RNGType r ) const;

  private:
   #ifdef HAS_MKL
   //! Register MKL RNGs into factory
   void regMKL( int nstream, const ctr::RNGMKLParameters& mklparam );
   #endif

   //! Register RNGSSE RNGs into factory
   void regRNGSSE( int nstream, const ctr::RNGSSEParameters& param );

   RNGFactory m_factory;        //!< Random nunmber generator factory
};

} // namespace tk

#endif // RNGStack_h
