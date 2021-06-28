// *****************************************************************************
/*!
  \file      src/Base/ContainerUtil.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Various STL container utilities
  \details   Various STL container utilities.
*/
// *****************************************************************************
#ifndef ContainerUtil_h
#define ContainerUtil_h

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <unordered_set>
#include <unordered_map>
#include <type_traits>
#include <sstream>

#include "Exception.hpp"

namespace tk {

//! Make elements of container unique (in-place, overwriting source container)
//! \param[in,out] c Container
template< class Container >
void
unique( Container& c )
{
  std::sort( begin(c), end(c) );
  auto it = std::unique( begin(c), end(c) );
  auto d = std::distance( begin(c), it );
  Assert( d >= 0, "Distance must be non-negative in tk::unique()" );
  c.resize( static_cast< std::size_t >( d ) );
}

//! Make elements of container unique (on a copy, leaving the source as is)
//! \param[in] src Container
//! \return Container containing only unique elements compared to src
template< class Container >
Container
uniquecopy( const Container& src )
{
  auto c = src;
  unique( c );
  return c;
}

//! \brief Find and return a constant reference to value for key in container
//!   that provides a find() member function with error handling
//! \param[in] map Map associating values to keys
//! \param[in] key Key to search for
//! \return A constant reference to the value associated to the key in map
//! \note If key is not found an exception is thrown.
template< typename Container >
auto cref_find( const Container& map, const typename Container::key_type& key )
noexcept(ndebug)
  -> const typename Container::mapped_type&
{
  const auto it = map.find( key );
  Assert( it != end(map), "Can't find key" );
  return it->second;
}

//! \brief Find and return a reference to value for key in a container that
//!   provides a find() member function with error handling
//! \param[in] map Map associating values to keys
//! \param[in] key Key to search for
//! \return A reference to the value associated to the key in map
//! \note If key is not found an exception is thrown.
template< typename Container >
auto ref_find( const Container& map, const typename Container::key_type& key )
noexcept(ndebug)
  -> typename Container::mapped_type&
{
  return const_cast< typename Container::mapped_type& >( cref_find(map,key) );
}

//! \brief Return minimum and maximum values of a vector
//! \param[in] vec Vector whose extents to compute
//! \return Array of two values with the minimum and maximum values
//! \note This function should not be called with heavy T types, as the a copy
//!   of a std::array< T, 2 > is created and returned.
template< typename T >
std::array< T, 2 >
extents( const std::vector< T >& vec )
{
  auto x = std::minmax_element( begin(vec), end(vec) );
  return {{ *x.first, *x.second }};
}

//! \brief Find and return minimum and maximum values in associative container
//! \param[in] map Map whose extents of values to find 
//! \return Array of two values with the minimum and maximum values in the map
//! \note This function should not be called with heavy Value types, as the a
//!   copy of a std::array< Value, 2 > is created and returned.
template< typename Container >
auto extents( const Container& map )
  -> std::array< typename Container::mapped_type, 2 >
{
  auto x = std::minmax_element( begin(map), end(map),
             []( const auto& a, const auto& b )
             { return a.second < b.second; } );
  return {{ x.first->second, x.second->second }};
}

//! Add all elements of an array to those of another one
//! \param[in,out] dst Destination array, i.e., left-hand side of a1 += a2
//! \param[in] src Source array, i.e., righ-hand side of a1 += a2
//! \return Destination containing a1[0] += a2[0], a1[1] += a2[1], ...
template< class T, std::size_t N >
std::array< T, N >&
operator+=( std::array< T, N >& dst, const std::array< T, N >& src ) {
  std::transform( src.cbegin(), src.cend(), dst.begin(), dst.begin(),
                  []( const T& s, T& d ){ return d += s; } );
  return dst;
}

//! Add all elements of a vector to those of another one
//! \param[in,out] dst Destination vector, i.e., left-hand side of v1 += v2
//! \param[in] src Source vector, i.e., righ-hand side of v1 += v2
//! \return Destination containing v1[0] += v2[0], v1[1] += v2[1], ...
//! \details If src.size() > dst.size() will grow dst to that of src.size()
//!   padding with zeros.
//! \note Will throw exception in DEBUG if src is empty (to warn on no-op), and
//!   if src.size() < dst.size() (to warn on loosing data).
template< class T, class Allocator >
std::vector< T, Allocator >&
operator+=( std::vector< T, Allocator >& dst,
            const std::vector< T, Allocator >& src )
{
  Assert( !src.empty(), "src empty in std::vector<T,Allocator>::operator+=()" );
  Assert( src.size() >= dst.size(), "src.size() < dst.size() would loose data "
          "in std::vector<T,Allocator>::operator+=()" );
  dst.resize( src.size() );
  std::transform( src.cbegin(), src.cend(), dst.begin(), dst.begin(),
                  []( const T& s, T& d ){ return d += s; } );
  return dst;
}

//! Divide all elements of a vector with those of another one
//! \param[in,out] dst Destination vector, i.e., left-hand side of v1 /= v2
//! \param[in] src Source vector, i.e., righ-hand side of v1 /= v2
//! \return Destination containing v1[0] /= v2[0], v1[1] /= v2[1], ...
//! \details If src.size() > dst.size() will grow dst to that of src.size()
//!   padding with zeros.
//! \note Will throw exception in DEBUG if src is empty (to warn on no-op), and
//!   if src.size() < dst.size() (to warn on loosing data).
template< class T, class Allocator >
std::vector< T, Allocator >&
operator/=( std::vector< T, Allocator >& dst,
            const std::vector< T, Allocator >& src )
{
  Assert( !src.empty(), "src empty in std::vector<T,Allocator>::operator/=()" );
  Assert( src.size() >= dst.size(), "src.size() < dst.size() would loose data "
          "in std::vector<T,Allocator>::operator/=()" );
  dst.resize( src.size() );
  std::transform( src.cbegin(), src.cend(), dst.begin(), dst.begin(),
                  []( const T& s, T& d ){ return d /= s; } );
  return dst;
}

//! Test if all keys of two associative containers are equal
//! \param[in] a 1st container to compare
//! \param[in] b 2nd container to compare
//! \return True if the containers have the same size and all keys (and only the
//!   keys) of the two containers are equal
//! \note It is an error to call this function with unequal-size containers,
//!   triggering an exception in DEBUG mode.
//! \note Operator != is used to compare the container keys.
template< class C1, class C2 >
bool keyEqual( const C1& a, const C2& b ) {
  Assert( a.size() == b.size(), "Size mismatch comparing containers" );
  std::set< typename C1::key_type > sorted_keys_of_a;
  for (const auto& c : a) sorted_keys_of_a.insert( c.first );
  std::set< typename C2::key_type > sorted_keys_of_b;
  for (const auto& c : b) sorted_keys_of_b.insert( c.first );
  return sorted_keys_of_a == sorted_keys_of_b;
}

//! Compute the sum of the sizes of a container of containers
//! \param[in] c Container of containers
//! \return Sum of the sizes of the containers of the container
template< class Container >
std::size_t sumsize( const Container& c ) {
  std::size_t sum = 0;
  // cppcheck-suppress useStlAlgorithm
  for (const auto& s : c) sum += s.size();
  return sum;
}

//! Compute the number of unique values in a container of containers
//! \param[in] c Container of containers
//! \return Number of unique values in a container of containers
template< class Container >
std::size_t numunique( const Container& c ) {
  using value_type = typename Container::value_type::value_type;
  static_assert( std::is_integral<value_type>::value,
    "Container::value_type::value_type must be an integral type." );
  std::unordered_set< value_type > u;
  for (const auto& r : c) u.insert( begin(r), end(r) );
  return u.size();
}

//! Compute the sum of the sizes of the values of an associative container
//! \tparam Map Container of containers type
//! \param[in] c Container of containers
//! \return Sum of the sizes of the values of the associative container
template< class Map >
std::size_t sumvalsize( const Map& c ) {
  std::size_t sum = 0;
  // cppcheck-suppress useStlAlgorithm
  for (const auto& s : c) sum += s.second.size();
  return sum;
}

//! Free memory of a container
//! \param[in] c Container defining a swap() member function
//! \details See http://stackoverflow.com/a/10465032 as to why this is done with
//!   the swap() member function of the container.
//! \see Specializations of std::swap are documented at
//!   http://en.cppreference.com/w/cpp/algorithm/swap
template< class Container >
void destroy( Container& c ) {
  typename std::remove_reference< decltype(c) >::type().swap( c );
}

//! Remove items from container based on predicate
//! \tparam Container Type of container to remove from
//! \tparam Predicate Type for functor defining the predicate
//! \param items Container object to remove from
//! \param predicate Predicate object instance to use
template< typename Container, typename Predicate >
void erase_if( Container& items, const Predicate& predicate ) {
  for ( auto it = items.begin(); it != items.end(); ) {
    if ( predicate(*it) ) it = items.erase(it);
    else ++it;
  }
}

//! Concatenate vectors of T
//! \tparam T Vector value type
//! \param[in,out] src Source vector (moved from)
//! \param[in,out] dst Destination vector
template< class T >
void concat( std::vector< T >&& src, std::vector< T >& dst )
{
  if (dst.empty())
    dst = std::move(src);
  else {
    dst.reserve( dst.size() + src.size() );
    std::move( std::begin(src), std::end(src), std::back_inserter(dst) );
    src.clear();
  }
}

//! Overwrite vectors of pair< bool, tk::real >
//! \tparam T Vector value type
//! \param[in,out] src Source vector (moved from)
//! \param[in,out] dst Destination vector
template< class T >
void concat( std::vector< std::pair< bool, T > >&& src,
             std::vector< std::pair< bool, T > >& dst )
{
  dst = std::move(src);
}

//! Concatenate unordered sets
//! \tparam Key Set key
//! \tparam Hash Set hasher
//! \tparam Eq Set equality operator
//! \param[in,out] src Source set (moved from)
//! \param[in,out] dst Destination set
template< class Key,
          class Hash = std::hash< Key >,
          class Eq = std::equal_to< Key > >
void concat( std::unordered_set< Key, Hash,Eq >&& src,
             std::unordered_set< Key, Hash, Eq >& dst )
{
  if (dst.empty())
    dst = std::move(src);
  else {
    dst.reserve( dst.size() + src.size() );
    std::move( std::begin(src), std::end(src), std::inserter(dst,end(dst)) );
    src.clear();
  }
}

//! Operator << for writing value_type of a standard map to output streams
//! \param[in,out] os Output stream to write to
//! \param[in] v value_type entry of a map
//! \return Updated output stream
template< class Key, class Value >
std::ostream&
operator<< ( std::ostream& os, const std::pair< const Key, Value  >& v ) {
  os << v.first << ':' << v.second;
  return os;
}

//! \brief Convert and return value as string
//! \tparam T Value type for input
//! \param[in] v Value for input to return as a string
//! \return String for input value
template< typename T >
std::string parameter( const T& v ) {
  std::stringstream s;
  s << v;
  return s.str();
}

//! \brief Convert and return values from container as string
//! \tparam V Container range for works on
//! \param[in] v Container whose components to return as a string
//! \return Concatenated string of values read from a container
template< typename V >
std::string parameters( const V& v ) {
  std::stringstream s;
  s << "{ ";
  for (auto p : v) s << p << ' ';
  s << "}";
  return s.str();
}

} // tk::

#endif // ContainerUtil_h
