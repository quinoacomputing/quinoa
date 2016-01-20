//******************************************************************************
/*!
  \file      src/Inciter/NodesReducer.C
  \author    J. Bakosi
  \date      Wed 20 Jan 2016 04:53:23 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Custom Charm++ reducer for merging global mesh node IDs across PEs
  \details   Custom Charm++ reducer for merging global mesh node IDs across PEs.
*/
//******************************************************************************

#include <iostream>     // NOT NEEDED!

#include "NodesReducer.h"
#include "Make_unique.h"
#include "PUPUtil.h"

namespace tk {

std::pair< int, std::unique_ptr<char[]> >
serialize( const std::vector< int >& pe,
           const std::vector< std::vector< std::size_t > >& gid )
//******************************************************************************
// Serialize global mesh nodes IDs to raw memory stream
//! \param[in] pe PE global node IDs coming from
//! \param[in] gid Global mesh node indices associated to each PE resulting from
//!   the contributing PE reading its contiguously-numbered mesh elements from
//!   file
//! \return Pair of the length and the raw stream containing the serialized
//!   global node IDs
//! \author J. Bakosi
//******************************************************************************
{
  // Prepare for serializing PE & node IDs to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | const_cast< std::vector< int >& >( pe );
  sizer | const_cast< std::vector< std::vector< std::size_t > >& >( gid );

  // Create raw character stream to store PE and the serialized global node IDs
  std::unique_ptr<char[]> flatData = tk::make_unique<char[]>( sizer.size() );

  // Serialize PE & node IDs
  PUP::toMem packer( flatData.get() );
  packer | const_cast< std::vector< int >& >( pe );
  packer | const_cast< std::vector< std::vector< std::size_t > >& >( gid );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

CkReductionMsg*
mergeNodes( int nmsg, CkReductionMsg **msgs )
//******************************************************************************
// Charm++ custom reducer for merging mesh node IDs during reduction across PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized node IDs
//! \return Aggregated PE + node IDs built for further aggregation if needed
//! \author J. Bakosi
//******************************************************************************
{
  // Will store deserialized PEs and global mesh node IDs
  std::vector< int > pe;
  std::vector< std::vector< std::size_t > > gid;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize global mesh node IDs from raw stream
  creator | pe;
  creator | gid;

  for (int m=1; m<nmsg; ++m) {
    // Unpack global PE mesh node IDs
    std::vector< int > p;
    std::vector< std::vector< std::size_t > > id;
    PUP::fromMem curCreator( msgs[m]->getData() );
    curCreator | p;
    curCreator | id;
    // Merge contributing PEs
    pe.insert( end(pe), std::make_move_iterator(begin(p)),
               std::make_move_iterator(end(p)) );
    // Merge global mesh node IDs
    gid.insert( end(gid), std::make_move_iterator(begin(id)),
                std::make_move_iterator(end(id)) );
  }

  // Serialize merged global mesh node IDs and associated PEs to raw stream
  auto stream = tk::serialize( pe, gid );

  // Forward serialized global mesh node IDs and associated PEs
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // tk::
