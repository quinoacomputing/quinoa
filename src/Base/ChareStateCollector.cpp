// *****************************************************************************
/*!
  \file      src/Base/ChareStateCollector.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare state collector group
  \details   Charm++ chare state collectory group used for debugging.
*/
// *****************************************************************************

#include <vector>
#include <unordered_map>

#include "ChareStateCollector.hpp"
#include "HashMapReducer.hpp"

namespace tk {

static CkReduction::reducerType stateMerger;

} // tk::

using tk::ChareStateCollector;

void
ChareStateCollector::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//! \details Since this is a [initnode] routine, the runtime system executes the
//!   routine exactly once on every logical node early on in the Charm++ init
//!   sequence. Must be static as it is called without an object. See also:
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  stateMerger =
    CkReduction::addReducer( tk::mergeHashMap< int, decltype(m_state) > );
}

void
ChareStateCollector::insert( const std::string& ch, int id, int pe, uint64_t it,
                             const std::string& fn )
// *****************************************************************************
//  Insert new state entry
//! \param[in] ch Chare name
//! \param[in] id Chare thisIndex
//! \param[in] pe Chare PE happens to reside on
//! \param[in] it Iteration count
//! \param[in] fn Chare member function name
// *****************************************************************************
{
  m_state.push_back( ChareState{{ ch, id, pe, it, fn, m_timer.dsec() }} );
}

void
ChareStateCollector::collect( bool error, CkCallback cb )
// *****************************************************************************
//  Collect chare state
//! \param[in] error If true we are called due to an error, if false, user is
//!   just curious
//! \param[in] cb Callback to use for the reduction
// *****************************************************************************
{
  // Will aggregate states across all PEs categorized by PEs
  std::unordered_map< int, std::vector< ChareState > > s{{ CkMyPe(), m_state }};

  // If we are called due to an error, encode a -1 chare id from which the
  // receiving side will know about the error.
  if (error) s[-1];

  // Aggregate states collected by this PE across all PEs
  auto stream = tk::serialize( s );
  contribute( stream.first, stream.second.get(), stateMerger, cb );
}

#include "NoWarning/charestatecollector.def.h"
