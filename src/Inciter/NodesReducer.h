//******************************************************************************
/*!
  \file      src/Inciter/NodesReducer.h
  \author    J. Bakosi
  \date      Wed 20 Jan 2016 04:53:05 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Custom Charm++ reducer for merging global mesh node IDs across PEs
  \details   Custom Charm++ reducer for merging global mesh node IDs across PEs.
*/
//******************************************************************************
#ifndef NodesReducer_h
#define NodesReducer_h

#include <unordered_map>
#include <unordered_set>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <charm++.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace tk {

//! Serialize global mesh nodes IDs to raw memory stream
std::pair< int, std::unique_ptr<char[]> >
serialize( const std::vector< int >& pe,
           const std::vector< std::vector< std::size_t > >& gid );

//! \brief Charm++ custom reducer for merging global node IDs during reduction
//!   across PEs
CkReductionMsg*
mergeNodes( int nmsg, CkReductionMsg **msgs );

} // tk::

#endif // NodesReducer_h
