// *****************************************************************************
/*!
  \file      src/Inciter/HoleReducer.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for aggregating hole data across PEs
  \details   Custom Charm++ reducer for aggregating hole data across PEs.
*/
// *****************************************************************************
#ifndef HoleReducer_h
#define HoleReducer_h

#include <vector>
#include <unordered_map>
#include <memory>
#include <utility>

#include "NoWarning/charm++.hpp"

#include "Types.hpp"

namespace inciter {

//! Serialize std::vector to raw memory stream
std::pair< int, std::unique_ptr<char[]> >
serialize( std::size_t meshid,
           const std::unordered_map< std::size_t,
                                     std::vector< tk::real > >& d );

//! Charm++ custom reducer for merging std::vectors during reduction across PEs
CkReductionMsg*
mergeHole( int nmsg, CkReductionMsg **msgs );

} // inciter::

#endif // HoleReducer_h
