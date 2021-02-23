// *****************************************************************************
/*!
  \file      src/Base/Callback.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Tagged tuple types used for passing Charm++ callbacks
  \details   Tagged tuple types used for passing Charm++ callbacks.
*/
// *****************************************************************************
#ifndef Callbacks_h
#define Callbacks_h

#include "NoWarning/charm++.hpp"

#include "Tags.hpp"
#include "TaggedTuple.hpp"

namespace tk {

using PartitionerCallback =
  tk::TaggedTuple< brigand::list<
      tag::load,           CkCallback
    , tag::partitioned,    CkCallback
    , tag::distributed,    CkCallback
    , tag::refinserted,    CkCallback
    , tag::refined,        CkCallback
  > >;

using RefinerCallback =
  tk::TaggedTuple< brigand::list<
      tag::queried,        CkCallback
    , tag::responded,      CkCallback
    , tag::compatibility,  CkCallback
    , tag::bndint,         CkCallback
    , tag::matched,        CkCallback
    , tag::refined,        CkCallback
  > >;

using SorterCallback =
  tk::TaggedTuple< brigand::list<
      tag::queried,        CkCallback
    , tag::responded,      CkCallback
    , tag::discinserted,   CkCallback
    , tag::workinserted,   CkCallback
  > >;

} // tk::

#endif // Callbacks_h
