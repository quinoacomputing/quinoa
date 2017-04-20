// *****************************************************************************
/*!
  \file      src/LoadBalance/LinearMap.h
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
#ifndef LinearMap_h
#define LinearMap_h

#include "NoWarning/linearmap.decl.h"

#include "Exception.h"

namespace tk {

//! Charm++ array map for initial placement of array elements in linear fashion
//! \details The map object is used by the Charm++ array manager to determine
//!   the "home" PE of each element. The home PE is the PE upon which the array
//!   element is initially placed, which will retain responsibility for
//!   maintaining the location of the element.
class LinearMap : public CkArrayMap {

  public:
    //! Constructor
    //! \param[in] nelem Total number of array elements
    explicit LinearMap( int nelem ) :
      m_chunksize( nelem > CkNumPes() ? nelem/CkNumPes() : 1 )
    { Assert( nelem > 0, "Number of array elements must be positive" ); }

    //! \brief Return the home processor number for the array element for linear
    //!   distribution
    int procNum( int, const CkArrayIndex& idx ) override;

    //! Create initial set of array elements based on linear distribution
    void populateInitial( int, CkArrayOptions& opt, void *msg, CkArrMgr *mgr )
    override;

  private:
    int m_chunksize;            //!< Number of array elements per PE
};

} // tk::

#endif // LinearMap_h
