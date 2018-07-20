// *****************************************************************************
/*!
  \file      src/Base/ContainerUtil.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Various STL container utilities
  \details   Various STL container utilities.
*/
// *****************************************************************************
#ifndef ContainerUtil_h
#define ContainerUtil_h

#include <vector>
#include <map>
#include <algorithm>
#include <iterator>

#include "Exception.h"

namespace tk {

template< class Container >
void
unique( Container& c )
// *****************************************************************************
//! Make elements of container unique (in-place, overwriting source container)
//! \param[inout] c Container
// *****************************************************************************
{
  std::sort( begin(c), end(c) );
  auto it = std::unique( begin(c), end(c) );
  auto d = std::distance( begin(c), it );
  Assert( d >= 0, "Distance must be non-negative in tk::unique()" );
  c.resize( static_cast< std::size_t >( d ) );
}

template< class Container >
Container
uniquecopy( const Container& src )
// *****************************************************************************
//! Make elements of container unique (on a copy, leaving the source as is)
//! \param[in] src Container
//! \return Container containing only unique elements compared to src
// *****************************************************************************
{
  auto c = src;
  unique( c );
  return c;
}

template< typename Container >
auto cref_find( const Container& map, const typename Container::key_type& key )
  -> const typename Container::mapped_type&
// *****************************************************************************
//! \brief Find and return a constant reference to value for key in container
//!   that provides a find() member function with error handling
//! \param[in] map Map associating values to keys
//! \param[in] key Key to search for
//! \return A constant reference to the value associated to the key in map
//! \note If key is not found an exception is thrown.
// *****************************************************************************
{
  const auto it = map.find( key );
  if (it != end(map)) return it->second; else Throw( "Can't find key" );
}

template< typename Container >
auto ref_find( const Container& map, const typename Container::key_type& key )
  -> typename Container::mapped_type&
// *****************************************************************************
//! \brief Find and return a reference to value for key in a container that
//!   provides a find() member function with error handling
//! \param[in] map Map associating values to keys
//! \param[in] key Key to search for
//! \return A reference to the value associated to the key in map
//! \note If key is not found an exception is thrown.
// *****************************************************************************
{
  return const_cast< typename Container::mapped_type& >( cref_find(map,key) );
}

template< typename T >
std::array< T, 2 >
extents( const std::vector< T >& vec )
// *****************************************************************************
//! \brief Return minimum and maximum values of a vector
//! \param[in] vec Vector whose extents to compute
//! \return Array of two values with the minimum and maximum values
//! \note This function should not be called with heavy T types, as the a copy
//!   of a std::array< T, 2 > is created and returned.
// *****************************************************************************
{
  auto x = std::minmax_element( begin(vec), end(vec) );
  return {{ *x.first, *x.second }};
}

template< typename Container >
auto extents( const Container& map )
  -> std::array< typename Container::mapped_type, 2 >
// *****************************************************************************
//! \brief Find and return minimum and maximum values in associative container
//! \param[in] map Map whose extents of values to find 
//! \return Array of two values with the minimum and maximum values in the map
//! \Note This function should not be called with heavy Value types, as the a
//!   copy of a std::array< Value, 2 > is created and returned.
// *****************************************************************************
{
  using pair_type = typename Container::value_type;
  auto x = std::minmax_element( begin(map), end(map),
             []( const pair_type& a, const pair_type& b )
             { return a.second < b.second; } );
  return {{ x.first->second, x.second->second }};
}

template< class T, class Allocator >
std::vector< T, Allocator >&
operator+=( std::vector< T, Allocator >& dst,
            const std::vector< T, Allocator >& src )
// *****************************************************************************
//! \brief Add all elements of a vector to another one
//! \param[inout] dst Destination vector, i.e., left-hand side of v1 += v2
//! \param[in] src Source vector, i.e., righ-hand side of v1 += v2
//! \return Destination containing v1[0] += v2[0], v1[1] += v2[1], ...
//! \details If src.size() > dst.size() will grow dst to that of src.size()
//!   padding with zeros.
//! \note Will throw exception in DEBUG if src is empty (to warn on no-op), and
//!   if src.size() < dst.size() (to warn on loosing data).
// *****************************************************************************
{
  Assert( !src.empty(), "src empty in std::vector<T,Allocator>::operator+=()" );
  Assert( src.size() >= dst.size(), "src.size() < dst.size() would loose data "
          "in std::vector<T,Allocator>::operator+=()" );
  dst.resize( src.size() );
  std::transform( src.begin(), src.end(), dst.begin(), dst.begin(),
                  []( const T& s, T& d ){ return d += s; } );
  return dst;
}

// *****************************************************************************
//! Test if all keys of two associative containers are equal
//! \param[in] a 1st container to compare
//! \param[in] b 2nd container to compare
//! \return True if the containers have the same size and all keys (and only the
//!   keys) of the two containers are equal
//! \note It is an error to call this function with unequal-size containers,
//!   triggering an exception in DEBUG mode.
//! \note Operator != is used to compare the container keys.
// *****************************************************************************
template< class C1, class C2 >
bool keyEqual( const C1& a, const C2& b ) {
  Assert( a.size() == b.size(), "Size mismatch comparing containers" );
  auto ia = a.cbegin();
  auto ib = b.cbegin();
  while (ia != a.cend()) {
    if (ia->first != ib->first) return false;
    ++ia;
    ++ib;
  }
  return true;
}

// *****************************************************************************
//! Compute the sum of the sizes of a container of containers
//! \param[in] c Container of containers
//! \return Sum of the sizes of the containers of the container
// *****************************************************************************
template< class Container >
std::size_t sumsize( const Container& c ) {
  std::size_t sum = 0;
  for (const auto& s : c) sum += s.size();
  return sum;
}

// *****************************************************************************
//! Compute the sum of the sizes of the values of an associative container
//! \tparam Map Container of containers type
//! \param[in] c Container of containers
//! \return Sum of the sizes of the values of the associative container
// *****************************************************************************
template< class Map >
std::size_t sumvalsize( const Map& c ) {
  std::size_t sum = 0;
  for (const auto& s : c) sum += s.second.size();
  return sum;
}

// *****************************************************************************
//! Free memory of a container
//! \param[in] c Container defining a swap() member function
//! \details See http://stackoverflow.com/a/10465032 as to why this is done with
//!   the swap() member function of the container.
//! \see Specializations of std::swap are documented at
//!   http://en.cppreference.com/w/cpp/algorithm/swap
// *****************************************************************************
template< class Container >
void destroy( Container& c ) {
  typename std::remove_reference< decltype(c) >::type().swap( c );
}

// *****************************************************************************
//! Remove items from container based on predicate
//! \tparam Container Type of container to remove from
//! \tparam Predicate Type for functor defining the predicate
//! \param items Container object to remove from
//! \param predicate Predicate object instance to use
// *****************************************************************************
template< typename Container, typename Predicate >
void erase_if( Container& items, const Predicate& predicate ) {
  for ( auto it = items.begin(); it != items.end(); ) {
    if ( predicate(*it) ) it = items.erase(it);
    else ++it;
  }
};

} // tk::

#endif // ContainerUtil_h
