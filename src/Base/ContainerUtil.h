//******************************************************************************
/*!
  \file      src/Base/ContainerUtil.h
  \author    J. Bakosi
  \date      Wed 17 Feb 2016 03:52:47 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Various STL container utilities
  \details   Various STL container utilities.
*/
//******************************************************************************
#ifndef ContainerUtil_h
#define ContainerUtil_h

#include <vector>
#include <map>
#include <algorithm>

#include "Exception.h"

namespace tk {

template< class Container >
void
unique( Container& c )
//******************************************************************************
//! Make elements of container unique
//! \param[inout] c Container
//! \author  J. Bakosi
//******************************************************************************
{
  std::sort( begin(c), end(c) );
  auto it = std::unique( begin(c), end(c) );
  auto d = std::distance( begin(c), it );
  Assert( d >= 0, "Distance must be non-negative in tk::unique()" );
  c.resize( static_cast< std::size_t >( d ) );
}

template< typename Container >
auto cref_find( const Container& map, const typename Container::key_type& key )
  -> const typename Container::mapped_type&
//******************************************************************************
//! \brief Find and return a constant reference to value for key in container
//!   that provides a find() member function with error handling
//! \param[in] map Map associating values to keys
//! \param[in] key Key to search for
//! \return A constant reference to the value associated to the key in map
//! \author J. Bakosi
//******************************************************************************
{
  const auto it = map.find( key );
  if (it != end(map)) return it->second; else Throw( "Can't find key" );
}

template< typename Container >
auto ref_find( const Container& map, const typename Container::key_type& key )
  -> typename Container::mapped_type&
//******************************************************************************
//! \brief Find and return a reference to value for key in a container that
//!   provides a find() member function with error handling
//! \param[in] map Map associating values to keys
//! \param[in] key Key to search for
//! \return A reference to the value associated to the key in map
//! \author J. Bakosi
//******************************************************************************
{
  return const_cast< typename Container::mapped_type& >( cref_find(map,key) );
}

template< typename T >
std::array< T, 2 >
extents( const std::vector< T >& vec )
//******************************************************************************
//! \brief Return minimum and maximum values of a vector
//! \param[in] vec Vector whose extents to compute
//! \return Array of two values with the minimum and maximum values
//! \Note This function should not be called with heavy T types, as the a copy
//!   of a std::array< T, 2 > is created and returned.
//! \author J. Bakosi
//******************************************************************************
{
  auto x = std::minmax_element( begin(vec), end(vec) );
  return {{ *x.first, *x.second }};
}

template< typename Container >
auto extents( const Container& map )
  -> std::array< typename Container::mapped_type, 2 >
//******************************************************************************
//! \brief Find and return minimum and maximum values in associative container
//! \param[in] map Map whose extents of values to find 
//! \return Array of two values with the minimum and maximum values in the map
//! \Note This function should not be called with heavy Value types, as the a
//!   copy of a std::array< Value, 2 > is created and returned.
//! \author J. Bakosi
//******************************************************************************
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
//******************************************************************************
//! \brief Add all elements of a vector to another one
//! \param[inout] dst Destination vector, i.e., left-hand side of v1 += v2
//! \param[in] src Source vector, i.e., righ-hand side of v1 += v2
//! \return Destination containing v1[0] += v2[0], v1[1] += v2[1], ...
//! \details If src.size() > dst.size() will grow dst to that of src.size()
//!   padding with zeros.
//! \note Will throw exception in DEBUG if src is empty (to warn on no-op), and
//!   if src.size() < dst.size() (to warn on loosing data).
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !src.empty(), "src empty in std::vector<T,Allocator>::operator+=()" );
  Assert( src.size() >= dst.size(), "src.size() < dst.size() would loose data "
          "in std::vector<T,Allocator>::operator+=()" );
  dst.resize( src.size() );
  std::transform( cbegin(src), cend(src), begin(dst), begin(dst),
                  []( const T& s, T& d ){ return d += s; } );
  return dst;
}

} // tk::

#endif // ContainerUtil_h
