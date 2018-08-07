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
  tk::tuple::tagged_tuple< tag::lboff,          CkCallback
                         , tag::stat,           CkCallback
                         , tag::load,           CkCallback
                         , tag::distributed,    CkCallback
                         , tag::refinserted,    CkCallback
                         , tag::refined,        CkCallback
                         , tag::flattened,      CkCallback
                         >;

using RefinerCallback =
  tk::tuple::tagged_tuple< tag::matched,        CkCallback
                         , tag::refined,        CkCallback
                         >;

using SolverCallback =
  tk::tuple::tagged_tuple< tag::nchare,         CkCallback
                         , tag::bounds,         CkCallback
                         , tag::com,            CkCallback
                         , tag::disccreated,    CkCallback
                         >;

using SorterCallback =
  tk::tuple::tagged_tuple< tag::flattened,      CkCallback
                         , tag::bounds,         CkCallback
                         , tag::discinserted,   CkCallback
                         , tag::workinserted,   CkCallback
                         >;

} // tk::

#endif // Callbacks_h
