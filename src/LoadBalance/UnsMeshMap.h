// *****************************************************************************
/*!
  \file      src/LoadBalance/UnsMeshMap.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Advanced Charm++ array creation with a map using an unstructured
             grid
  \details   Advanced Charm++ array creation refers to various ways arrays can
     be created with the Charm++ runtime system. See
     http://charm.cs.illinois.edu/manuals/html/charm++/manual.html, Sec.
     Advanced Array Creation. This class does a distribution that is based on
     which portion of a distributed sparse matrix resulting from discretization
     on an unstructured grid residing on a PE should hold a given chare array
     element. (The one that owns most on-PE rows to minimize off-PE
     communication.)

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
#ifndef UnsMeshMap_h
#define UnsMeshMap_h

#include <vector>
#include <cstddef>

#include "NoWarning/unsmeshmap.decl.h"

namespace tk {

//! \brief Charm++ array map for initial placement of array elements using an
//!   unstructured grid
//! \details The map object is used by the Charm++ array manager to determine
//!   the "home" PE of each element. The home PE is the PE upon which the array
//!   element is initially placed, which will retain responsibility for
//!   maintaining the location of the element.
class UnsMeshMap : public CkArrayMap {

  public:
    //! Constructor
    explicit
    UnsMeshMap( std::size_t npoin,
                const std::vector< std::vector< std::size_t > >& point );

    //! \brief Return the home processor number for the array element based on
    //!   unstructured-mesh-aware distribution computed in the constructor
    int procNum( int, const CkArrayIndex& idx ) override;

    //! \brief Create initial set of array elements based on the
    //!   unstructured-mesh-aware distribution computed in the constructor
    void populateInitial( int, CkArrayOptions& opt, void *msg, CkArrMgr *mgr )
    override;

  private:
    std::vector< std::size_t > m_pe;    //!< PE assigned to each array element

    //! Check that all PEs create at least a single array element, fix if not
    void fixPEs();
};

} // tk::

#endif // UnsMeshMap_h
