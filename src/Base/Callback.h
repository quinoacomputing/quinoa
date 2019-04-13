// *****************************************************************************
/*!
  \file      src/Base/Callback.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Tagged tuple types used for passing Charm++ callbacks
  \details   Tagged tuple types used for passing Charm++ callbacks.
*/
// *****************************************************************************
#ifndef Callbacks_h
#define Callbacks_h

#include "NoWarning/charm++.h"

#include "Tags.h"
#include "TaggedTuple.h"

namespace tk {

using PartitionerCallback =
  tk::tuple::tagged_tuple< tag::load,           CkCallback
                         , tag::distributed,    CkCallback
                         , tag::refinserted,    CkCallback
                         , tag::refined,        CkCallback
                         >;

using RefinerCallback =
  tk::tuple::tagged_tuple< tag::edges,          CkCallback
                         , tag::compatibility,  CkCallback
                         , tag::bndint,         CkCallback
                         , tag::matched,        CkCallback
                         , tag::refined,        CkCallback
                         >;

using SorterCallback =
  tk::tuple::tagged_tuple< tag::queried,        CkCallback
                         , tag::responded,      CkCallback
                         , tag::discinserted,   CkCallback
                         , tag::workinserted,   CkCallback
                         >;

} // tk::

#endif // Callbacks_h
