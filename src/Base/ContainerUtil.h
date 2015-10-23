//******************************************************************************
/*!
  \file      src/Base/ContainerUtil.h
  \author    J. Bakosi
  \date      Fri 23 Oct 2015 06:06:47 AM MDT
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

template< typename T >
std::vector< std::pair< std::string, T > >
average( const std::map< std::string, std::vector< T > >& mapvec,
         const std::string& addendum = "" )
//******************************************************************************
//  Compute average of values of vector in std::map::mapped_type 
//! \param[in] mapvec Map of vectors associated to labels
//! \param[in] addendum Optional string to add to the label
//! \return Vector of pairs of label and average
//! \author J. Bakosi
//******************************************************************************
{
  std::vector< std::pair< std::string, T > > s;
  for (const auto& t : mapvec) {
    T sum = 0.0;
    for (const auto& v : t.second) sum += v;
    s.emplace_back( t.first + addendum, sum/static_cast<T>(t.second.size()) );
  }
  return s;
}

template< typename T >
std::vector< std::pair< std::string, T > >
variance( const std::map< std::string, std::vector< T > >& mapvec,
          const std::vector< std::pair< std::string, T > >& avg,
          const std::string& addendum = "" )
//******************************************************************************
//  Compute variance of values of vector in std::map::mapped_type 
//! \param[in] mapvec Map of vectors associated to labels
//! \param[in] avg Vector of labels and averages (labels unused)
//! \param[in] addendum Optional string to add to the label
//! \return Vector of pairs of label and variance
//! \author J. Bakosi
//******************************************************************************
{
  Assert( mapvec.size() == avg.size(),
          "Map and vector must be equal size for variance calculation" );

  std::vector< std::pair< std::string, T > > s;
  std::size_t i = 0;
  for (const auto& t : mapvec) {
    T sum = 0.0;
    for (const auto& v : t.second) sum += (v-avg[i].second)*(v-avg[i].second);
    s.emplace_back( t.first + addendum, sum/static_cast<T>(t.second.size()) );
    ++i;
  }

  return s;
}

template< typename Key, typename Value >
Value
val( const std::map< Key, Value >& map, Key key )
//******************************************************************************
//! Find and return a copy of value for key in std::map with error handling
//! \param[in] map Map associating values to keys
//! \param[in] key
//! \return A copy of the value associated to the key in map
//! \Note This function should not be called with heavy Key types, as the key is
//!    passed by value.
//! \author J. Bakosi
//******************************************************************************
{
  const auto it = map.find( key );

  if (it != map.end())
    return it->second;
  else
    Throw( "Can't find key " + std::to_string(key) );
}

template< typename Key, typename Value >
const Value&
ref( const std::map< Key, Value >& map, Key key )
//******************************************************************************
//! \brief Find and return a constant reference to value for key in std::map with
//!   error handling
//! \param[in] map Map associating values to keys
//! \param[in] key
//! \return A constant reference to the value associated to the key in map
//! \Note This function should not be called with heavy Key types, as the key is
//!   passed by value.
//! \author J. Bakosi
//******************************************************************************
{
  const auto it = map.find( key );

  if (it != map.end())
    return it->second;
  else
    Throw( "Can't find key " + std::to_string(key) );
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

template< typename Key, typename Value >
std::array< Value, 2 >
extents( const std::map< Key, Value >& map )
//******************************************************************************
//! \brief Find and return minimum and maximum values in std::map
//! \param[in] map Map whose extents of values to find 
//! \return Array of two values with the minimum and maximum values in the map
//! \Note This function should not be called with heavy Value types, as the a
//!   copy of a std::array< Value, 2 > is created and returned.
//! \author J. Bakosi
//******************************************************************************
{
  using pair_type = std::pair< const Key, Value >;
  auto x = std::minmax_element( begin(map), end(map),
             []( const pair_type& a, const pair_type& b )
             { return a.second < b.second; } );
  return {{ x.first->second, x.second->second }};
}

} // tk::

#endif // ContainerUtil_h
