// *****************************************************************************
/*!
  \file      src/Base/HashSetReducer.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for merging std::unordered_sets across PEs
  \details   Custom Charm++ reducer for merging std::unordered_sets across PEs.
*/
// *****************************************************************************
#ifndef HashSetReducer_h
#define HashSetReducer_h

#include <unordered_set>
#include <memory>

#include "NoWarning/charm++.hpp"

namespace tk {

//! Serialize std::unordered_set to raw memory stream
//! \tparam Key Set key/value type
//! \tparam Hash Set hasher
//! \tparam Eq Set equality operator
//! \param[in] m Hash set to serialize
//! \return Pair of the length and the raw stream containing the serialized set
template< class Key,
          class Hash = std::hash< Key >,
          class Eq = std::equal_to< Key > >
std::pair< int, std::unique_ptr<char[]> >
serialize( const std::unordered_set< Key, Hash, Eq >& m ) {
   // Prepare for serializing set to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | const_cast< std::unordered_set< Key, Hash, Eq >& >( m );

  // Create raw character stream to store the serialized set
  std::unique_ptr<char[]> flatData = std::make_unique<char[]>( sizer.size() );

  // Serialize set, each message will contain a set
  PUP::toMem packer( flatData.get() );
  packer | const_cast< std::unordered_set< Key, Hash, Eq >& >( m );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

//! \brief Charm++ custom reducer for merging std::unordered_sets during
//!   reduction across PEs
//! \tparam Key Set key/value type
//! \tparam Hash Set hasher
//! \tparam Eq Set equality operator
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized sets
//! \return Aggregated std::unordered_sets built for further aggregation if
//!    needed
//! \details During aggregation the set keys/values are inserted, i.e., the keys
//!   remain unique.
template< class Key,
          class Hash = std::hash< Key >,
          class Eq = std::equal_to< Key > >
CkReductionMsg*
mergeHashSet( int nmsg, CkReductionMsg **msgs ) {
  // Will store deserialized set
  std::unordered_set< Key, Hash, Eq > p;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize set from raw stream
  creator | p;

  for (int m=1; m<nmsg; ++m) {
    // Unpack set
    std::unordered_set< Key, Hash, Eq > u;
    PUP::fromMem curCreator( msgs[m]->getData() );
    curCreator | u;
    // Concatenate sets
    p.insert( begin(u), end(u) );
  }

  // Serialize concatenated sets to raw stream
  auto stream = tk::serialize( p );

  // Forward serialized hash set
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // tk::

#endif // HashSetReducer_h
