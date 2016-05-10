// *****************************************************************************
/*!
  \file      src/RNG/RNGStack.h
  \author    J. Bakosi
  \date      Sun 03 Apr 2016 10:06:07 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Stack of random number generators
  \details   This file declares class RNGStack, which implements various
    functionality related to registering and instantiating random number
    generators interfacing to multiple RNG libraries. Registration and
    instantiation use a random number generator factory, which is a std::map (an
    associative container), associating unique RNG keys to their constructor
    calls. For more details, see the in-code documentation of the constructor.
*/
// *****************************************************************************
#ifndef RNGStack_h
#define RNGStack_h

#include <map>
#include <vector>
#include <functional>
#include <type_traits>

#include "RNG.h"
#include "RNGParam.h"
#include "Options/RNG.h"

namespace tk {

//! Random number generator factory: keys associated to their constructors
//! \author J. Bakosi
using RNGFactory = std::map< ctr::RNGType, std::function< RNG() > >;

//! \brief Random number generator stack
//! \author J. Bakosi
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
