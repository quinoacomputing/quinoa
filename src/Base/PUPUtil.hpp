// *****************************************************************************
/*!
  \file      src/Base/PUPUtil.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ Pack/UnPack utilities
  \details   This file contains some extensions to Charm++'s Pack/UnPack
    routines.
*/
// *****************************************************************************
#ifndef PUPUtil_h
#define PUPUtil_h

#include <type_traits>

#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <optional>

#include "NoWarning/pup_stl.hpp"

//! Extensions to Charm++'s Pack/Unpack routines
namespace PUP {

//////////////////// Serialize enum class ////////////////////

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wuninitialized"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wuninitialized"
  #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

//! \brief Pack/Unpack enum class.
//! \details In Charm++ usually both the pup() overload and an overload for
//!   operator| are defined for all serializable types. However, we cannot
//!   define operator| for enum class as it would conflict with Charm++'s
//!   catch-all, template< class T > inline void operator|(PUP::er &p,T &t)
//!   {...}.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] e Enum class to pack/unpack
template< typename E,
          typename std::enable_if< std::is_enum< E >::value, int >::type = 0 >
inline void pup( PUP::er& p, E& e ) {
  auto v = static_cast< typename std::underlying_type< E >::type >( e );
  p | v;
  e = static_cast< E >( v );
}

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

//////////////////// Serialize std::unordered_map ////////////////////

//! Pack/Unpack std::unordered_map.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] m std::unordered_map< Key, T, Hash, KeyEqual > to pack/unpack
template< class Key,
          class T,
          class Hash = std::hash< Key >,
          class KeyEqual = std::equal_to< Key > >
inline void pup( PUP::er& p, std::unordered_map< Key, T, Hash, KeyEqual >& m ) {
  auto size = PUP_stl_container_size( p, m );
  if (p.isUnpacking()) {
    for (decltype(size) s=0; s<size; ++s) {
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
//! \param[in] m std::unordered_map< Key, T, Hash, KeyEqual > to pack/unpack
template< class Key,
          class T,
          class Hash = std::hash< Key >,
          class KeyEqual = std::equal_to< Key > >
inline void operator|( PUP::er& p,
                       std::unordered_map< Key, T, Hash, KeyEqual >& m )
{ pup( p, m ); }

//////////////////// Serialize std::unordered_set ////////////////////

//! Pack/Unpack std::unordered_set.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] s std::unordered_set< Key, Hash, KeyEqual > to pack/unpack
template< class Key,
          class Hash = std::hash< Key >,
          class KeyEqual = std::equal_to< Key > >
inline void pup( PUP::er& p, std::unordered_set< Key, Hash, KeyEqual >& s ) {
  auto size = PUP_stl_container_size( p, s );
  if (p.isUnpacking()) {
    for (decltype(size) i=0; i<size; ++i) {
      Key node;
      p | node;
      s.emplace( node );
    }
  } else {
    for (auto& t : s) {
      Key node( t );
      p | node;
    }
  }
}
//! Pack/Unpack std::unordered_set.
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] s std::unordered_set< Key, Hash, KeyEqual > to pack/unpack
template< class Key,
          class Hash = std::hash< Key >,
          class KeyEqual = std::equal_to< Key > >
inline void operator|( PUP::er& p,
                       std::unordered_set< Key, Hash, KeyEqual >& s )
{ pup( p, s ); }

//////////////////// Serialize std::optional ////////////////////

//! Pack/Unpack std::optional
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] o std::optional< T > of arbitrary type T to pack/unpack
template< class T >
inline void pup( PUP::er& p, std::optional< T >& o ) {
  T underlying_value = o ? *o : T();
  bool exist = o ? true : false;
  p | exist;
  p | underlying_value;
  o = exist ? std::make_optional(underlying_value) : std::nullopt;
}
//! Pack/Unpack std::optional
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] o std::optional< T > of arbitrary type T to pack/unpack
template< class T >
inline void operator|( PUP::er& p, std::optional< T >& o ) { pup( p, o ); }

//////////////////// Serialize std::variant ////////////////////

// Since std::variant when default-constructed is initialized to hold a value of
// the first alternative of its type list, calling PUP that works based on a
// std::visit with a templated operator() would always incorrectly trigger the
// overload for the first type. Thus when PUPing a variant not only its value
// but its type must also be sent during migration. The pup operator template
// below achieves this by reading out not only the value but also its zero-based
// index of the type alternative that is currently held by the variant passed to
// its initializer constructor. The index and the variant are then PUPed and
// when unpacking, as an additional step, the variant is reconstructed using the
// index and the value in the variant. This latter is done by invoking an
// expansion of an initializer list, guaranteed to happen in order, stepping
// through the typelist in the variant.  Thanks to Nils Deppe for simplifying
// the original version of this operation. See UnitTest/tests/Base/TestPUPUtil.h
// or Inciter::SchemeBase.h for puping a variant in action.

//! Pack/Unpack helper for std::variant
//! \param[in,out] index Counter (location) for type in variant
//! \param[in] send_index Target counter (location) for type in variant
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] var std::variant< Ts... > of arbitrary types to pack/unpack
template <class T, class... Ts>
char pup_helper( std::size_t& index,
                 const std::size_t send_index,
                 PUP::er& p,
                 std::variant<Ts...>& var )
{
  if (index == send_index) {
    if (p.isUnpacking()) {
      T t{};
      p | t;
      var = std::move(t);
    } else {
      p | std::get<T>(var);
    }
  }
  index++;
  return '0';
}

//! Pack/Unpack std::variant
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] var std::variant< Ts... > of arbitrary types to pack/unpack
template <class... Ts>
void pup(PUP::er& p, std::variant<Ts...>& var) {
  std::size_t index = 0;
  auto send_index = var.index();
  p | send_index;
  (void)std::initializer_list<char>{
      pup_helper<Ts>(index, send_index, p, var)...};
}

//! Pack/Unpack std::variant
//! \param[in] p Charm++'s pack/unpack object
//! \param[in] d std::variant< Ts... > of arbitrary types to pack/unpack
template <typename... Ts>
inline void operator|(PUP::er& p, std::variant<Ts...>& d) {
  pup(p, d);
}

} // PUP::

#endif // PUPUtil_h
