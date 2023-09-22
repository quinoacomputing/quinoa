// *****************************************************************************
/*!
  \file      src/Inciter/DiagReducer.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for merging std::vectors across PEs
  \details   Custom Charm++ reducer for merging std::vectors across PEs.
*/
// *****************************************************************************

#include <stddef.h>
#include <type_traits>
#include <memory>

#include "DiagReducer.hpp"
#include "Diagnostics.hpp"
#include "Exception.hpp"

namespace inciter {

std::pair< int, std::unique_ptr<char[]> >
serialize( std::size_t meshid,
           std::size_t ncomp,
           const std::vector< std::vector< tk::real > >& d )
// *****************************************************************************
// Serialize std::vectors to raw memory stream
//! \param[in] meshid Mesh ID
//! \param[in] ncomp Number of scalar components being solved
//! \param[in] d Diagnostics vector of vectors (of eq components)
//! \return Pair of the length and the raw stream containing the serialized
//!   vectors
// *****************************************************************************
{
  // Prepare for serializing diagnostics to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | meshid;
  sizer | ncomp;
  sizer | const_cast< std::vector< std::vector< tk::real > >& >( d );

  // Create raw character stream to store the serialized vectors
  std::unique_ptr<char[]> flatData = std::make_unique<char[]>( sizer.size() );

  // Serialize vector, each message will contain a vector
  PUP::toMem packer( flatData.get() );
  packer | meshid;
  packer | ncomp;
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
  std::size_t meshid, ncomp;
  std::vector< std::vector< tk::real > > v;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize vector from raw stream
  creator | meshid;
  creator | ncomp;
  creator | v;

  for (int m=1; m<nmsg; ++m) {
    // Unpack vector
    std::size_t mid, nc;
    std::vector< std::vector< tk::real > > w;
    PUP::fromMem curCreator( msgs[m]->getData() );
    curCreator | mid;
    curCreator | nc;
    curCreator | w;
    // Aggregate diagnostics vector
    meshid = mid;
    ncomp = nc;
    Assert( v.size() == w.size(),
            "Size mismatch during diagnostics aggregation" );
    Assert( v.size() == inciter::NUMDIAG,
            "Size mismatch during diagnostics aggregation" );
    for (std::size_t i=0; i<v.size(); ++i)
      Assert( v[i].size() == w[i].size(),
              "Size mismatch during diagnostics aggregation" );
    // Apply diagnostics aggregation policy
    // Sum for L2 normal of the numerical solution for all scalar components
    for (std::size_t i=0; i<v[L2SOL].size(); ++i) v[L2SOL][i] += w[L2SOL][i];
    // Sum for the L2 norm of the numerical - analytical solution for all comps
    for (std::size_t i=0; i<v[L2ERR].size(); ++i) v[L2ERR][i] += w[L2ERR][i];
    // Sum for the L2 norm of the residual of all components
    for (std::size_t i=0; i<v[L2RES].size(); ++i) v[L2RES][i] += w[L2RES][i];
    // Max for the Linf norm of the numerical - analytical solution for all comp
    for (std::size_t i=0; i<v[LINFERR].size(); ++i)
      if (w[LINFERR][i] > v[LINFERR][i]) v[LINFERR][i] = w[LINFERR][i];
    // Sum of the total energy over the entire domain
    v[TOTALSOL][0] += w[TOTALSOL][0];
    // Copy ITER, TIME, DT
    for (std::size_t j=v.size()-3; j<v.size(); ++j)
      for (std::size_t i=0; i<v[j].size(); ++i)
        v[j][i] = w[j][i];
  }

  // Serialize concatenated diagnostics vector to raw stream
  auto stream = serialize( meshid, ncomp, v );

  // Forward serialized diagnostics
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // inciter::
