// *****************************************************************************
/*!
  \file      src/Base/HashMapReducer.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Custom Charm++ reducer for merging std::unordered_maps across PEs
  \details   Custom Charm++ reducer for merging std::unordered_maps across PEs.
*/
// *****************************************************************************
#ifndef HashMapReducer_h
#define HashMapReducer_h

#include <unordered_map>

#include "NoWarning/charm++.h"

#include "Make_unique.h"
#include "ContainerUtil.h"

namespace tk {

//! Serialize std::unordered_map to raw memory stream
//! \param[in] m Hash map to serialize
//! \return Pair of the length and the raw stream containing the serialized map
//! \author J. Bakosi
template< class Key,
          class T,
          class Hash = std::hash< Key >,
          class KeyEqual = std::equal_to< Key > >
std::pair< int, std::unique_ptr<char[]> >
serialize( const std::unordered_map< Key, T, Hash, KeyEqual >& m ) {
   // Prepare for serializing map to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | const_cast< std::unordered_map< Key, T, Hash, KeyEqual >& >( m );

  // Create raw character stream to store the serialized map
  std::unique_ptr<char[]> flatData = tk::make_unique<char[]>( sizer.size() );

  // Serialize map, each message will contain a map
  PUP::toMem packer( flatData.get() );
  packer | const_cast< std::unordered_map< Key, T, Hash, KeyEqual >& >( m );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

//! \brief Charm++ custom reducer for merging std::unordered_maps during
//!   reduction across PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized maps
//! \return Aggregated std::unordered_maps built for further aggregation if
//!    needed
//! \details During aggregation the map keys are inserted, i.e., the keys remain
//!   unique and the mapped values, assuming containers defining begin() and
//!   end() iterators() are concatenated.
//! \note The mapped type must be a container, i.e., must provide iterators
//!   begin() and end().
//! \author J. Bakosi
template< class Key,
          class T,
          class Hash = std::hash< Key >,
          class KeyEqual = std::equal_to< Key > >
CkReductionMsg*
mergeHashMap( int nmsg, CkReductionMsg **msgs ) {
  // Will store deserialized map
  std::unordered_map< Key, T, Hash, KeyEqual > p;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize map from raw stream
  creator | p;

  for (int m=1; m<nmsg; ++m) {
    // Unpack map
    std::unordered_map< Key, T, Hash, KeyEqual > u;
    PUP::fromMem curCreator( msgs[m]->getData() );
    curCreator | u;
    // Concatenate maps
    for (const auto& c : u) {
      auto& b = p[ c.first ];
      b.insert( end(b), begin(c.second), end(c.second) );
    }
  }

  // Serialize concatenated maps to raw stream
  auto stream = tk::serialize( p );

  // Forward serialized hash map
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // tk::

#endif // HashMapReducer_h
