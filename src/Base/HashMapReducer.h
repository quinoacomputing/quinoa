// *****************************************************************************
/*!
  \file      src/Base/HashMapReducer.h
  \copyright 2012-2015, J. Bakosi, 2016-2019, Triad National Security, LLC.
  \brief     Custom Charm++ reducer for merging std::unordered_maps across PEs
  \details   Custom Charm++ reducer for merging std::unordered_maps across PEs.
*/
// *****************************************************************************
#ifndef HashMapReducer_h
#define HashMapReducer_h

#include <vector>
#include <unordered_map>

#include "NoWarning/charm++.h"

#include "Make_unique.h"
#include "ContainerUtil.h"

namespace tk {

//! Serialize std::unordered_map to raw memory stream
//! \tparam Key Map key
//! \tparam T Map value
//! \tparam Hash Map hasher
//! \tparam Eq Map equality operator
//! \param[in] m Hash map to serialize
//! \return Pair of the length and the raw stream containing the serialized map
template< class Key,
          class T,
          class Hash = std::hash< Key >,
          class Eq = std::equal_to< Key > >
std::pair< int, std::unique_ptr<char[]> >
serialize( const std::unordered_map< Key, T, Hash, Eq >& m ) {
   // Prepare for serializing map to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | const_cast< std::unordered_map< Key, T, Hash, Eq >& >( m );

  // Create raw character stream to store the serialized map
  std::unique_ptr<char[]> flatData = tk::make_unique<char[]>( sizer.size() );

  // Serialize map, each message will contain a map
  PUP::toMem packer( flatData.get() );
  packer | const_cast< std::unordered_map< Key, T, Hash, Eq >& >( m );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

//! Concatenate vectors of T
//! \tparam T Vector value type
//! \param[in] src Source vector
//! \param[in,out] dst Destination vector
template< class T >
void concat( const std::vector< T >& src, std::vector< T >& dst ) {
  dst.insert( end(dst), begin(src), end(src) );
}

//! Overwrite vectors of pair< bool, tk::real >
//! \tparam T Vector value type
//! \param[in] src Source vector
//! \param[in,out] dst Destination vector
template< class T >
void concat( const std::vector< std::pair< bool, T > >& src,
             std::vector< std::pair< bool, T > >& dst ) {
  dst = src;
}

//! Concatenate unordered sets
//! \tparam Key Set key
//! \tparam Hash Set hasher
//! \tparam Eq Set equality operator
//! \param[in] src Source set
//! \param[in,out] dst Destination set
template< class Key,
          class Hash = std::hash< Key >,
          class Eq = std::equal_to< Key > >
void concat( const std::unordered_set< Key, Hash,Eq >& src,
             std::unordered_set< Key, Hash, Eq >& dst )
{
  dst.insert( begin(src), end(src) );
}

//! \brief Charm++ custom reducer for merging std::unordered_maps during
//!   reduction across PEs
//! \tparam Key Map key
//! \tparam T Map value
//! \tparam Hash Map hasher
//! \tparam Eq Map equality operator
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized maps
//! \return Aggregated std::unordered_maps built for further aggregation if
//!    needed
//! \details During aggregation the map keys are inserted, i.e., the keys remain
//!   unique and the mapped values, assuming containers defining begin() and
//!   end() iterators() are concatenated.
//! \note The mapped type must be a container, i.e., must provide iterators
//!   begin() and end().
template< class Key,
          class T,
          class Hash = std::hash< Key >,
          class Eq = std::equal_to< Key > >
CkReductionMsg*
mergeHashMap( int nmsg, CkReductionMsg **msgs ) {
  // Will store deserialized map
  std::unordered_map< Key, T, Hash, Eq > p;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize map from raw stream
  creator | p;

  for (int m=1; m<nmsg; ++m) {
    // Unpack map
    std::unordered_map< Key, T, Hash, Eq > u;
    PUP::fromMem curCreator( msgs[m]->getData() );
    curCreator | u;
    // Concatenate maps
    for (const auto& c : u) concat( c.second, p[c.first] );
  }

  // Serialize concatenated maps to raw stream
  auto stream = tk::serialize( p );

  // Forward serialized hash map
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // tk::

#endif // HashMapReducer_h
