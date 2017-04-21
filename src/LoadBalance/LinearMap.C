// *****************************************************************************
/*!
  \file      src/LoadBalance/LinearMap.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Advanced Charm++ array creation with a map in a linear fashion
  \details   Advanced Charm++ array creation refers to various ways arrays can
     be created with the Charm++ runtime system. See
     http://charm.cs.illinois.edu/manuals/html/charm++/manual.html, Sec.
     Advanced Array Creation. This class does a simple linear distribution.

     Note that to help with performance, it is not advised to do heavy
     computations in the overridden member functions, procNum() and
     populateInitial(), since they can be potentially called many times. See
     also this note on the inner workings of the Charm++ runtime system
     regarding map objects and array element creation from long-time Charm++
     developer, Eric Bohm, also at:
     http://lists.cs.uiuc.edu/pipermail/charm/2015-May/002047.html

     _Procnum will be called to construct your chare, the first time a message
     is sent to it from a node, and each time subsequent sends do not find a
     precached location record for it. The latter event can occur when many
     sends have pushed it out of cache, or after migration._

     _The potential global memory footprint of location management caching is
     proportional to the total number of objects multiplied by the number of
     nodes. Therefore, the runtime system keeps a finite number on each node.
     At the limit, procnum could be called for nearly every message send,
     therefore procnum should be designed to be inexpensive._

     The heavy portion of array element placement should therefore be done in
     the constructor.
*/
// *****************************************************************************

#include <string>

#include "NoWarning/charm.h"
#include "LinearMap.h"

using tk::LinearMap;

int
LinearMap::procNum( int, const CkArrayIndex& idx )
// *****************************************************************************
//  Return the home processor number for the array element for linear
//  distribution
//! \param[in] idx Charm++ array index object containing the array element index
//!   to assign a PE to
//! \return PE assigned
//! \author J. Bakosi
// *****************************************************************************
{
  int elem = *idx.data();       // array element we assign PE for
  auto pe = elem / m_chunksize;
  if (pe >= CkNumPes()) pe = CkNumPes()-1;

  Assert( pe < CkNumPes(), "Assigned PE (" + std::to_string(pe) +
          ") larger than NumPEs (" + std::to_string(CkNumPes()) + ")" );

  return pe;
}

void
LinearMap::populateInitial( int, CkArrayOptions& opt, void *msg, CkArrMgr *mgr )
// *****************************************************************************
// Create initial set of array elements based on linear distribution
//! \param[in] opt Charm++ array options object containing the number of initial
//!   array elements to be created
//! \param[in] msg Charm++ messsage to use for array element creation
//! \param[in] mgr Array manager to use
//! \author J. Bakosi
// *****************************************************************************
{
  int nelem = *opt.getNumInitial().data(); // number of array elements requested
  if (nelem == 0) return;                  // no initial elements requested

  auto lower = CkMyPe() * m_chunksize;
  auto upper = lower + m_chunksize;
  auto remainder = nelem % CkNumPes();
  if (remainder && CkMyPe() == CkNumPes()-1) upper += remainder;

  for (int e=0; e<nelem; ++e)
    if (e >= lower && e < upper)
      mgr->insertInitial( CkArrayIndex(e), CkCopyMsg(&msg) );

  mgr->doneInserting();
  CkFreeMsg( msg );
}

#include "NoWarning/linearmap.def.h"
