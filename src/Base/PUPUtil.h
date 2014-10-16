//******************************************************************************
/*!
  \file      src/Base/PUPUtil.h
  \author    J. Bakosi
  \date      Tue 16 Sep 2014 08:22:45 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Charm++ Pack/UnPack utilities
  \details   Charm++ Pack/UnPack utilities
*/
//******************************************************************************
#ifndef PUPUtil_h
#define PUPUtil_h

#include <unordered_map>

#include <pup_stl.h>
#include <CharmUtil.h>

namespace PUP {

//! Pack/Unpack enum class
template< typename E,
          typename std::enable_if< std::is_enum< E >::value, int >::type = 0 >
inline void pup( PUP::er& p, E& e ) {
  using underlying_type = typename std::underlying_type< E >::type;
  underlying_type v = static_cast< underlying_type >( e );
  p | v;
  e = static_cast< E >( v );
}
// Cannote define operator| for enum class as it would conflict with Charm++'s
// catch-all, template<class T> inline void operator|(PUP::er &p,T &t) {...}

//! PUP_tuple_impl: specialization for empty std::tuple
template< std::size_t I = 0, typename... Tp >
typename std::enable_if< I == sizeof...(Tp), void >::type
pup_tuple_impl( PUP::er&, std::tuple< Tp... >& ) {}
//! PUP_tuple_impl: specialization for non-empty std::tuple
template< std::size_t I = 0, typename... Tp >
typename std::enable_if< I < sizeof...(Tp), void >::type
pup_tuple_impl( PUP::er& p, std::tuple< Tp... >& t ) {
  p | std::get<I>(t);
  pup_tuple_impl< I + 1, Tp... >( p, t );
}
//! Pack/Unpack std::tuple
template< typename... Ts >
inline void pup( PUP::er& p, std::tuple< Ts... >& t ) {
  if (p.isUnpacking()) t = std::tuple< Ts... >();
  pup_tuple_impl( p, t );
}
template< typename... Ts >
inline void operator|( PUP::er& p, std::tuple< Ts... >& t ) { pup( p, t ); }

//! Pack/Unpack std::array
template< class T, std::size_t N >
inline void pup( PUP::er& p, std::array< T, N >& a ) {
  for (std::size_t s=0; s<N; ++s) p | a[s];
}
template< class T, std::size_t N >
inline void operator|( PUP::er& p, std::array< T, N >& a ) { pup( p, a ); }

//! Pack/Unpack std::unordered_map
template< class Key, typename T, class KeyEqual = std::equal_to< Key > >
inline void pup( PUP::er& p, std::unordered_map< Key, T, KeyEqual >& m ) {
  auto size = PUP_stl_container_size( p, m );
  if (p.isUnpacking()) {
    for (std::size_t s=0; s<size; ++s) {
      std::pair< Key, T > node;
      p | node;
      m.emplace( node );
    }
  } else {
    for (auto& t : m) {
      std::pair< Key, T > node( t );
      p | node;
    }
  }
}
template< class Key, typename T, class KeyEqual = std::equal_to< Key > >
inline void operator|( PUP::er& p, std::unordered_map< Key, T, KeyEqual >& m )
{ pup( p, m ); }

} // PUP::

#endif // PUPUtil_h
