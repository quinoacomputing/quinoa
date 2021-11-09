// *****************************************************************************
/*!
  \file      src/Base/Has.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     "Has-a" utilities for detecting class internals
  \details   "Has-a" utilities for detecting class internals
*/
// *****************************************************************************
#ifndef Has_h
#define Has_h

#include <type_traits>

namespace tk {

//! Detect if a type defines type 'alias'
template< typename, typename = std::void_t<> >
struct HasTypedef_alias : std::false_type {};

template< typename T >
struct HasTypedef_alias< T, std::void_t< typename T::alias > >
  : std::true_type {};

template < typename T >
inline constexpr bool HasTypedef_alias_v = HasTypedef_alias< T >::value;


//! Detect if a type defines type 'i_am_tagged_tuple'
template< typename, typename = std::void_t<> >
struct HasTypedef_i_am_tagged_tuple : std::false_type {};

template< typename T >
struct HasTypedef_i_am_tagged_tuple< T,
  std::void_t< typename T::i_am_tagged_tuple > > : std::true_type {};

template < typename T >
inline constexpr bool HasTypedef_i_am_tagged_tuple_v =
  HasTypedef_i_am_tagged_tuple< T >::value;


//! Detect if a type defines function 'expect::description()'
template< typename, typename = std::void_t<> >
struct HasFunction_expect_description : std::false_type {};

template< typename T >
struct HasFunction_expect_description< T,
  std::void_t< decltype(std::declval<typename T::expect>().description()) > >
  : std::true_type {};

template < typename T >
inline constexpr bool HasFunction_expect_description_v =
  HasFunction_expect_description< T >::value;


//! Detect if a type defines variable 'expect::lower'
template< typename, typename = std::void_t<> >
struct HasVar_expect_lower : std::false_type {};

template< typename T >
struct HasVar_expect_lower< T,
  std::void_t< decltype(std::declval<typename T::expect>().lower) > >
  : std::true_type {};

template < typename T >
inline constexpr bool HasVar_expect_lower_v = HasVar_expect_lower< T >::value;


//! Detect if a type defines variable 'expect::upper'
template< typename, typename = std::void_t<> >
struct HasVar_expect_upper : std::false_type {};

template< typename T >
struct HasVar_expect_upper< T,
  std::void_t< decltype(std::declval<typename T::expect>().upper) > >
  : std::true_type {};

template < typename T >
inline constexpr bool HasVar_expect_upper_v = HasVar_expect_upper< T >::value;


//! Detect if a type defines function 'expect::choices()'
template< typename, typename = std::void_t<> >
struct HasFunction_expect_choices : std::false_type {};

template< typename T >
struct HasFunction_expect_choices< T,
  std::void_t< decltype(std::declval<typename T::expect>().choices()) > >
  : std::true_type {};

template < typename T >
inline constexpr bool HasFunction_expect_choices_v =
  HasFunction_expect_choices< T >::value;

} // tk::

#endif // Has_h
