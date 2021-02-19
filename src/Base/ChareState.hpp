// *****************************************************************************
/*!
  \file      src/Base/ChareState.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare state collector group
  \details   Charm++ chare state collectory group used for debugging.
*/
// *****************************************************************************
#ifndef ChareState_h
#define ChareState_h

#include <string>

#include "Tags.hpp"
#include "Types.hpp"
#include "TaggedTuple.hpp"

namespace tk {

//! Chare state
using ChareState = TaggedTuple< brigand::list<
    tag::ch,   std::string   // chare name
  , tag::id,   int           // thisIndex
  , tag::pe,   int           // PE
  , tag::it,   uint64_t      // iteration count
  , tag::fn,   std::string   // member function name
  , tag::time, tk::real      // time stamp
> >;

} // tk::

#endif // ChareState_h
