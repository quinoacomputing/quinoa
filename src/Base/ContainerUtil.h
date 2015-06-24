//******************************************************************************
/*!
  \file      src/Base/ContainerUtil.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:10:02 PM MDT
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
typename std::map< Key, Value >::size_type
lid( const std::map< Key, Value >& map, Key key )
//******************************************************************************
//! Find and return value for key in std::map with error handling
//! \param[in] map Map associating values to keys
//! \param[in] key
//! \return Value associated to key in map
//! \Note This functions should not be instantiated with heavy Key types, as it
//!   is passed by value.
//! \author J. Bakosi
//******************************************************************************
{
  const auto it = map.find( key );

  if (it != map.end())
    return it->second;
  else
    Throw( "Can't find key " + std::to_string(key) );
}

} // tk::

#endif // ContainerUtil_h
