// *****************************************************************************
/*!
  \file      src/Base/Callbacks.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
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
                         , tag::created,        CkCallback
                         , tag::flattened,      CkCallback
                         , tag::avecost,        CkCallback
                         , tag::stdcost,        CkCallback
                         , tag::coord,          CkCallback
                         >;

using RefinerCallback =
  tk::tuple::tagged_tuple< tag::matched,        CkCallback
                         , tag::refined,        CkCallback
                         >;

using SolverCallback =
  tk::tuple::tagged_tuple< tag::com,            CkCallback
                         , tag::coord,          CkCallback
                         >;

} // tk::

#endif // Callbacks_h
