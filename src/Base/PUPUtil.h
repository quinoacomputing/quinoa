//******************************************************************************
/*!
  \file      src/Base/PUPUtil.h
  \author    J. Bakosi
  \date      Sat 17 Jan 2015 07:27:41 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Charm++ Pack/UnPack utilities
  \brief     This file contains some extensions to Charm++'s Pack/UnPack
    routines.
*/
//******************************************************************************
#ifndef PUPUtil_h
#define PUPUtil_h

#include <iostream>     // NOT NEEDED!
#include <unordered_map>

#include <boost/optional.hpp>

#include <pup_stl.h>

#include <CharmUtil.h>

//! Extensions to Charm++'s Pack/Unpack routines
namespace PUP {

//////////////////// Serialize enum class ////////////////////

//! Pack/Unpack enum class. In Charm++ usually both the pup() overload and an
//! overload for operator| are defined for all serializable types. However, we
//! cannot define operator| for enum class as it would conflict with Charm++'s
//! catch-all, template< class T > inline void operator|(PUP::er &p,T &t) {...}.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] e Enum class to pack/unpack
//! \author J. Bakosi
template< typename E,
          typename std::enable_if< std::is_enum< E >::value, int >::type = 0 >
inline void pup( PUP::er& p, E& e ) {
  using underlying_type = typename std::underlying_type< E >::type;
  underlying_type v = static_cast< underlying_type >( e );
  p | v;
  e = static_cast< E >( v );
}

//////////////////// Serialize std::tuple ////////////////////

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
//! Pack/Unpack std::tuple.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] t std::tuple to pack/unpack
//! \author J. Bakosi
template< typename... Ts >
inline void pup( PUP::er& p, std::tuple< Ts... >& t ) {
  if (p.isUnpacking()) t = std::tuple< Ts... >();
  pup_tuple_impl( p, t );
}
//! Pack/Unpack std::tuple.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] t std::tuple to pack/unpack
//! \author J. Bakosi
template< typename... Ts >
inline void operator|( PUP::er& p, std::tuple< Ts... >& t ) { pup( p, t ); }

//////////////////// Serialize std::array ////////////////////

//! Pack/Unpack std::array.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] a std::array< T, N > of arbitrary type T to pack/unpack
//! \author J. Bakosi
template< class T, std::size_t N >
inline void pup( PUP::er& p, std::array< T, N >& a ) {
  for (std::size_t s=0; s<N; ++s) p | a[s];
}
//! Pack/Unpack std::array.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] a std::array< T, N > of arbitrary type T to pack/unpack
//! \author J. Bakosi
template< class T, std::size_t N >
inline void operator|( PUP::er& p, std::array< T, N >& a ) { pup( p, a ); }

//////////////////// Serialize std::unordered_map ////////////////////

//! Pack/Unpack std::unordered_map.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] m std::unordered_map< Key, T, KeyEqual > to pack/unpack
//! \author J. Bakosi
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
//! Pack/Unpack std::unordered_map.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] m std::unordered_map< Key, T, KeyEqual > to pack/unpack
//! \author J. Bakosi
template< class Key, typename T, class KeyEqual = std::equal_to< Key > >
inline void operator|( PUP::er& p, std::unordered_map< Key, T, KeyEqual >& m )
{ pup( p, m ); }

//////////////////// Serialize boost::optional ////////////////////

//! Pack/Unpack boost::optional.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] o boost::optional< T > of arbitrary type T to pack/unpack
//! \author J. Bakosi
template< class T >
inline void pup( PUP::er& p, boost::optional< T >& o ) {
  // IMPLEMENT !!!
  std::cout << "pup boost::opt\n";
}
//! Pack/Unpack boost::optional.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] o boost::optional< T > of arbitrary type T to pack/unpack
//! \author J. Bakosi
template< class T >
inline void operator|( PUP::er& p, boost::optional< T >& o ) { pup( p, o ); }

} // PUP::

#endif // PUPUtil_h
