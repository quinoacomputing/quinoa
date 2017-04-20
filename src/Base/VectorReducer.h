// *****************************************************************************
/*!
  \file      src/Statistics/VectorReducer.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Custom Charm++ reducer for merging std::vectors across PEs
  \details   Custom Charm++ reducer for merging std::vectors across PEs.
*/
// *****************************************************************************
#ifndef VectorReducer_h
#define VectorReducer_h

#include <vector>

#include "NoWarning/charm++.h"

namespace tk {

//! Serialize std::vector to raw memory stream
std::pair< int, std::unique_ptr<char[]> >
serialize( const std::vector< std::size_t >& v );

//! Charm++ custom reducer for merging std::vectors during reduction across PEs
CkReductionMsg*
mergeVector( int nmsg, CkReductionMsg **msgs );

} // tk::

#endif // VectorReducer_h
