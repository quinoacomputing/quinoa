// *****************************************************************************
/*!
  \file      src/Inciter/DiagReducer.C
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Custom Charm++ reducer for merging std::vectors across PEs
  \details   Custom Charm++ reducer for merging std::vectors across PEs.
*/
// *****************************************************************************

#include "DiagReducer.h"
#include "Make_unique.h"
#include "ContainerUtil.h"

namespace tk {

std::pair< int, std::unique_ptr<char[]> >
serialize( const std::vector< std::vector< tk::real > >& d )
// *****************************************************************************
// Serialize std::vectors to raw memory stream
//! \param[in] d Diagnostics vector of vectors (of eq components)
//! \return Pair of the length and the raw stream containing the serialized
//!   vectors
// *****************************************************************************
{
  // Prepare for serializing diagnostics to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | const_cast< std::vector< std::vector< tk::real > >& >( d );

  // Create raw character stream to store the serialized vectors
  std::unique_ptr<char[]> flatData = tk::make_unique<char[]>( sizer.size() );

  // Serialize vector, each message will contain a vector
  PUP::toMem packer( flatData.get() );
  packer | const_cast< std::vector< std::vector< tk::real > >& >( d );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

CkReductionMsg*
mergeDiag( int nmsg, CkReductionMsg **msgs )
// *****************************************************************************
// Charm++ custom reducer for merging diagnostics during reduction across PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized
//!   diagnostics
//! \return Aggregated diagnostics built for further aggregation if needed
// *****************************************************************************
{
  // Will store deserialized diagnostics vector of vectors
  std::vector< std::vector< tk::real > > v;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize vector from raw stream
  creator | v;

  for (int m=1; m<nmsg; ++m) {
    // Unpack vector
    std::vector< std::vector< tk::real > > w;
    PUP::fromMem curCreator( msgs[m]->getData() );
    curCreator | w;
    // Aaggregate diagnostics vector
    Assert( v.size() == w.size(),
            "Size mismatch during diagnostics aggregation" );
    for (std::size_t i=0; i<v.size(); ++i)
      Assert( v[i].size() == w[i].size(),
              "Size mismatch during diagnostics aggregation" );
    // Apply diagnostics aggregation policy, see also LinSysMerger::updateDiag()
    // Sum for L2 normal of the numerical solution for all scalar components
    for (std::size_t i=0; i<v[0].size(); ++i) v[0][i] += w[0][i];
    // Sum for the L2 norm of the numerical - analytical solution for all comps
    for (std::size_t i=0; i<v[1].size(); ++i) v[1][i] += w[1][i];
    // Max for the Linf norm of the numerical - analytical solution for all comp
    for (std::size_t i=0; i<v[2].size(); ++i)
      if (w[2][i] > v[2][i]) v[2][i] = w[2][i];
  }

  // Serialize concatenated diagnostics vector to raw stream
  auto stream = tk::serialize( v );

  // Forward serialized diagnostics
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // tk::
