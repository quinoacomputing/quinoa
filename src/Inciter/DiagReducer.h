// *****************************************************************************
/*!
  \file      src/Inciter/DiagReducer.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Custom Charm++ reducer for merging diagnostics across PEs
  \details   Custom Charm++ reducer for merging diagnostics across PEs.
*/
// *****************************************************************************
#ifndef DiagReducer_h
#define DiagReducer_h

#include <vector>

#include "NoWarning/charm++.h"

#include "Types.h"

namespace inciter {

//! Serialize std::vector to raw memory stream
std::pair< int, std::unique_ptr<char[]> >
serialize( const std::vector< std::vector< tk::real > >& d );

//! Charm++ custom reducer for merging std::vectors during reduction across PEs
CkReductionMsg*
mergeDiag( int nmsg, CkReductionMsg **msgs );

} // inciter::

#endif // DiagReducer_h
