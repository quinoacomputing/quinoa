//******************************************************************************
/*!
  \file      src/Base/PUPUtil.h
  \author    J. Bakosi
  \date      Thu 15 May 2014 08:03:53 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Charm++ Pack/UnPack utilities
  \details   Charm++ Pack/UnPack utilities
*/
//******************************************************************************
#ifndef PUPUtil_h
#define PUPUtil_h

#include <pup_stl.h>

//! PUP operator for enum class RNGType, see also Charm's pup.h
inline void operator|( PUP::er& p, tk::ctr::RNGType& e ) {
  using underlying_type = std::underlying_type< tk::ctr::RNGType >::type;
  underlying_type v = static_cast< underlying_type >( e );
  p | v;
  e = static_cast< tk::ctr::RNGType >( v );
}

//! PUP_tuple: specialization for empty std::tuple
template< std::size_t I = 0, typename... Tp >
inline typename std::enable_if< I == sizeof...(Tp), void >::type
PUP_tuple( PUP::er& p, std::tuple< Tp... >& t ) {}

//! PUP_tuple: specialization for non-empty std::tuple
template< std::size_t I = 0, typename... Tp >
inline typename std::enable_if< I < sizeof...(Tp), void >::type
PUP_tuple( PUP::er& p, std::tuple< Tp... >& t ) {
  p | std::get< I >( t );
  PUP_tuple< I + 1, Tp... >( p, t );
}

//! PUP operator for std::tuple
template< typename... Ts >
inline void operator|( PUP::er& p, std::tuple< Ts... >& t ) {
  if (p.isUnpacking()) t = std::tuple< Ts... >();
  PUP_tuple( p, t );
}

#endif // PUPUtil_h
