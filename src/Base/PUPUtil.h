//******************************************************************************
/*!
  \file      src/Base/PUPUtil.h
  \author    J. Bakosi
  \date      Thu 03 Jul 2014 04:57:32 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Charm++ Pack/UnPack utilities
  \details   Charm++ Pack/UnPack utilities
*/
//******************************************************************************
#ifndef PUPUtil_h
#define PUPUtil_h

#include <make_unique.h>
#include <pup_stl.h>
#include <CharmUtil.h>

namespace tk {

//! Pack/Unpack enum class
template< typename E >
inline void pup( PUP::er& p, E& e ) {
  using underlying_type = typename std::underlying_type< E >::type;
  underlying_type v = static_cast< underlying_type >( e );
  p | v;
  e = static_cast< E >( v );
}

//! Use the enum class pack/unpack only for strongly typed enums
template< typename T,
          typename std::enable_if< std::is_enum<T>::value, int >::type = 0 >
inline void operator|( PUP::er& p, T& e ) { pup(p,e); }

//! Delegate pack/unpack to parent scope operator| for non-strongly-typed enums
template< typename T,
          typename std::enable_if< !std::is_enum<T>::value, int >::type = 0 >
inline void operator|( PUP::er& p, T& e ) { ::operator|(p,e); }

//! Pack/Unpack std::pair
template< class A, class B >
inline void pup( PUP::er& p, typename std::pair< A, B >& v ) {
  p | v.first;
  p | v.second;
}

// Pack/Unpack the length of an STL container
template< class Container >
inline std::size_t pup_container_size( PUP::er& p, Container &c ) {
  auto size = c.size();
  p | size;
  return size;
}

//! Pack/Unpack std::vector of T
template< class T >
inline void operator|( PUP::er& p, std::vector< T >& v ) {
  auto size = pup_container_size( p, v );
  if (p.isUnpacking()) v.resize( size );
  for (auto& t : v) p | t;
}

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
inline void pup_tuple( PUP::er& p, std::tuple< Ts... >& t ) {
  if (p.isUnpacking()) t = std::tuple< Ts... >();
  pup_tuple_impl( p, t );
}
template< typename... Ts >
inline void operator|( PUP::er& p, std::tuple< Ts... >& t ) { pup_tuple(p,t); }

//! Pack/Unpack std::map
template< typename KeyType, typename ValueType >
inline void operator|( PUP::er& p, std::map< KeyType, ValueType >& m ) {
  auto size = pup_container_size( p, m );
  if (p.isUnpacking()) {
    for (std::size_t s=0; s<size; ++s) {
      std::pair< KeyType, ValueType > node;
      pup( p, node );
      m.emplace( node );
    }
  } else {
    for (auto& t : m) {
      std::pair< KeyType, ValueType > node( t );
      pup( p, node );
    }
  }
}

// //! PUP std::unique_ptr using make_unqiue and constructor arguments
// template< class T, class Deleter, typename... Args >
// void pup_unique_ptr( PUP::er& p,
//                      std::unique_ptr< T, Deleter >& t,
//                      Args&... args ) {
//   if (p.isUnpacking()) {
//     std::unique_ptr< T, Deleter > r =
//       tk::make_unique< T, Args... >( std::move(args)... );
//     p | *r;
//     t = std::move(r);
//   } else {
//     p | *t;
//   }
// }

// //! PUP std::unique_ptr using a std::unique_ptr to a creator function
// template< class T, class Deleter >
// void pup_unique_ptr( PUP::er& p,
//                      std::unique_ptr< T, Deleter >& t,
//                      std::unique_ptr< T, Deleter >&& creator ) {
//   if (p.isUnpacking()) {
//     std::unique_ptr< T, Deleter > r = std::move(creator);
//     p | *r;
//     t = std::move(r);
//   } else {
//     p | *t;
//   }
// }

// //! PUP std::vector of std::unqiue_ptr
// template< class T, typename... Args >
// void pup_vector_unique_ptr( PUP::er& p,
//                             std::vector< std::unique_ptr<T> >& v,
//                             Args&... args )
// {
//   auto size = pup_container_size( p, v );
//   if (p.isUnpacking()) v.resize( size );
//   for (auto& t : v) pup_unique_ptr( p, t, std::move(args)... );
// }


// //! PUP std::map of std::unqiue_ptr
// template< class K, class T, typename... Args >
// void pup_map_unique_ptr( PUP::er& p,
//                          std::map< K, std::unique_ptr<T> >& m,
//                          Args&... args )
// {
//   auto size = pup_container_size( p, m );
//   if ( p.isUnpacking() )
//     for (std::size_t s=0; s<size; ++s)
//       m.insert( std::pair< K, T >() );
//   for (auto& t : m) pup_unique_ptr( p, t, std::move(args)... );
// }

} // tk::

#endif // PUPUtil_h
