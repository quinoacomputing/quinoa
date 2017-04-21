// *****************************************************************************
/*!
  \file      src/Inciter/FieldsMerger.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Custom Charm++ reducer for merging mesh node indices across PEs
  \details   Custom Charm++ reducer for merging mesh node indices across PEs.
*/
// *****************************************************************************
#ifndef FieldsMerger_h
#define FieldsMerger_h

#include <vector>

#include "NoWarning/charm++.h"

#include "Make_unique.h"
#include "ContainerUtil.h"

namespace inciter {

//! Serialize mesh node indices categorized by chares to raw memory stream
//! \param[in] m Chare mesh node indices to serialize
//! \return Pair of the length and the raw stream containing the serialized data
//! \author J. Bakosi
template< typename T >
std::pair< int, std::unique_ptr<char[]> >
serialize( const std::vector< std::pair< int, T > >& m ) {
   // Prepare for serializing node indices to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer |
   const_cast< std::vector< std::pair< int, T > >& >( m );

  // Create raw character stream to store the serialized node indices
  std::unique_ptr<char[]> flatData = tk::make_unique<char[]>( sizer.size() );

  // Serialize, each message will contain a list of node indices per chare
  PUP::toMem packer( flatData.get() );
  packer |
   const_cast< std::vector< std::pair< int, T > >& >( m );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

//! \brief Charm++ custom reducer for merging mesh node indices categorized by
//!   chares during reduction across PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized mesh
//!   node indices categorized by chares
//! \return Aggregated mesh node indices categorized by chares built for further
//!   aggregation if needed
//! \details During aggregation the outer vector is simply concatenated.
//! \author J. Bakosi
template< typename T >
CkReductionMsg*
mergeFields( int nmsg, CkReductionMsg **msgs ) {
  // Will store deserialized mesh node indices categorized by chares
  std::vector< std::pair< int, T > > p;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize from raw stream
  creator | p;

  for (int m=1; m<nmsg; ++m) {
    // Unpack mesh node indices
    std::vector< std::pair< int, T > > u;
    PUP::fromMem curCreator( msgs[m]->getData() );
    curCreator | u;
    // Concatenate mesh node indices categorized by chares
    p.insert( end(p), begin(u), end(u) );
  }

  // Serialize concatenated node indices categorized by chares to raw stream
  auto stream = serialize( p );

  // Forward serialized mesh node indices categorized by chares
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // inciter::

#endif // FieldsMerger_h
