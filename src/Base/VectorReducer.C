// *****************************************************************************
/*!
  \file      src/Statistics/VectorReducer.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Custom Charm++ reducer for merging std::vectors across PEs
  \details   Custom Charm++ reducer for merging std::vectors across PEs.
*/
// *****************************************************************************

#include "VectorReducer.h"
#include "Make_unique.h"
#include "ContainerUtil.h"

namespace tk {

std::pair< int, std::unique_ptr<char[]> >
serialize( const std::vector< std::size_t >& v )
// *****************************************************************************
// Serialize std::vectors to raw memory stream
//! \param[in] v Vector
//! \return Pair of the length and the raw stream containing the serialized
//!   vectors
//! \author J. Bakosi
// *****************************************************************************
{
  // Prepare for serializing vectors to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | const_cast< std::vector< std::size_t >& >( v );

  // Create raw character stream to store the serialized vectors
  std::unique_ptr<char[]> flatData = tk::make_unique<char[]>( sizer.size() );

  // Serialize vector, each message will contain a vector
  PUP::toMem packer( flatData.get() );
  packer | const_cast< std::vector< std::size_t >& >( v );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

CkReductionMsg*
mergeVector( int nmsg, CkReductionMsg **msgs )
// *****************************************************************************
// Charm++ custom reducer for merging std::vectors during reduction across PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized vectors
//! \return Aggregated std::vectors built for further aggregation if needed
//! \author J. Bakosi
// *****************************************************************************
{
  // Will store deserialized vector
  std::vector< std::size_t > v;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize vector from raw stream
  creator | v;

  for (int m=1; m<nmsg; ++m) {
    // Unpack vector
    std::vector< std::size_t > u;
    PUP::fromMem curCreator( msgs[m]->getData() );
    curCreator | u;
    // Concatenate vector and make it unique
    v.insert( end(v), begin(u), end(u) );
    tk::unique( v );
  }

  // Serialize concatenated vectors to raw stream
  auto stream = tk::serialize( v );

  // Forward serialized vector
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // tk::
