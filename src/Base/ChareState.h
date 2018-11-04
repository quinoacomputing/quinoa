// *****************************************************************************
/*!
  \file      src/Base/ChareState.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ chare state collector group
  \details   Charm++ chare state collectory group used for debugging.
*/
// *****************************************************************************
#ifndef ChareState_h
#define ChareState_h

#include <string>

#include "Tags.h"
#include "Types.h"
#include "TaggedTuple.h"

namespace tk {

//! Chare state
using ChareState = tuple::tagged_tuple<
                     tag::ch,   std::string     // chare name
                   , tag::id,   int             // thisIndex
                   , tag::pe,   int             // PE
                   , tag::it,   uint64_t        // iteration count
                   , tag::fn,   std::string     // member function name
                   , tag::time, tk::real        // time stamp
                   >;

} // tk::

#endif // ChareState_h
