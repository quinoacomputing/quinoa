// *****************************************************************************
/*!
  \file      src/LoadBalance/UnsMeshMap.C
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

#include <algorithm>

#include "NoWarning/charm.h"
#include "Exception.h"
#include "UnsMeshMap.h"

using tk::UnsMeshMap;

UnsMeshMap::UnsMeshMap( std::size_t npoin,
                        const std::vector< std::vector< std::size_t > >& point )
 : m_pe()
// *****************************************************************************
// Constructor
//! \param[in] npoin Total number of points in mesh
//! \param[in] point Global mesh point ids owned by each array element
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( npoin > 0, "Need at least a single mesh point" );
  Assert( point.size() >= static_cast< std::size_t >( CkNumPes() ),
          "UnsMeshMap only works with nchare >= numPEs" );

  // Vector of maps to associate the number of mesh points an array element
  // contributes to a pe to the array element
  std::vector< std::map< std::size_t, std::size_t > > owner( point.size() );

  // Compute number of mesh points per PE (chunksize)
  const auto numpes = static_cast< std::size_t >( CkNumPes() );
  const auto chunksize = npoin>numpes ? npoin/numpes : 1;

  // Count up number of mesh points contributing to a PE for each array element
  for (std::size_t e=0; e<point.size(); ++e)    // for all array elements
    for (auto p : point[e]) {           // for all points array element e owns
      auto pe = p / chunksize;          // PE array element e contributes to
      if (pe == static_cast<std::size_t>(CkNumPes())) --pe;
      ++owner[e][pe];  // count up number of points element e contributing to pe
    }

  // Find and store the PE with the largest number of points owned by an array
  // element. This is the key (PE) of the largest mapped_type of the element's
  // map, owner[e].
  using pr = decltype(owner)::value_type::value_type;
  m_pe.resize( owner.size() );
  for (std::size_t e=0; e<owner.size(); ++e) {
    const auto most = std::max_element( begin(owner[e]), end(owner[e]),
                                        []( const pr& p1, const pr& p2 )
                                          { return p1.second < p2.second; } );
    m_pe[e] = most->first;
  }

  // Check if all PEs have at least one array element to create, fix if not
  fixPEs();
}

void
UnsMeshMap::fixPEs()
// *****************************************************************************
//  Check that all PEs create at least a single array element, fix if not
//! \details This is required because if the array elements, placed using this
//!   map object, send messages to some Charm++ chare group branches and the
//!   group happens to use Charm++'s structured dagger, such as LinSysMerger,
//!   memory problems will occur.
//! \author J. Bakosi
// *****************************************************************************
{
  // Build unique set of PEs the array elements are assigned to
  std::set< std::size_t > nkind;
  for (auto p : m_pe) nkind.insert( p );

  // If not all PEs have at least one array element to create, go through all
  // PEs that have no array elements assigned, and construct a map that
  // associates the number of PEs to array elements. Then as a replacement PE,
  // pick the PE that has the most array elements assigned. Then assign the
  // first of the elements assigned to the PE with the most work to the
  // replacement PE.
  if (nkind.size() != static_cast<std::size_t>(CkNumPes()))
    for (std::size_t p=0; p<static_cast<std::size_t>(CkNumPes()); ++p)
      if (!nkind.count(p)) {    // if PE p has no array elements to create
        // Count up number elements for each PE
        std::map< std::size_t, std::size_t > npe;
        using pr = decltype(npe)::value_type;
        for (auto q : m_pe) ++npe[q];
        // Find the PE with the most array elements assigned
        const auto most = std::max_element( begin(npe), end(npe),
                                            []( const pr& p1, const pr& p2 )
                                            { return p1.second < p2.second; } );
        for (std::size_t q=0; q<m_pe.size(); ++q)
          if (m_pe[q] == most->first) { // first element of the most loaded PE
            //std::cout << CkMyPe() << ": replace " << m_pe[q] << " by " << p
            //          << " for e " << q << std::endl;
            m_pe[q] = p;        // replace PE of m_pe[q] of element q by p
            q = m_pe.size();    // quit (enough to replace one)
          }
      }

  nkind.clear();
  for (auto p : m_pe) nkind.insert( p );
  Assert( nkind.size() == static_cast<std::size_t>(CkNumPes()),
          "There are still PE(s) with no array elements assigned" );
}

int
UnsMeshMap::procNum( int, const CkArrayIndex& idx )
// *****************************************************************************
//  Return the home processor number for the array element based on
//  unstructured-mesh-aware distribution computed in the constructor
//! \param[in] idx Charm++ array index object containing the array element index
//!   to assign a PE to
//! \return PE assigned
//! \author J. Bakosi
// *****************************************************************************
{
  // array element we assign PE for
  auto elem = static_cast<std::size_t>(*idx.data());

  Assert( elem < m_pe.size(),
          "Array index larger than the number of chares for which the "
          "UnsMeshMap object holds indices for" );

  auto pe = static_cast< int >( m_pe[ elem ] );

  Assert( pe < CkNumPes(), "Assigned PE (" + std::to_string(pe) +
          ") larger than NumPEs (" + std::to_string(CkNumPes()) + ")" );

  return pe;
}

void
UnsMeshMap::populateInitial( int, CkArrayOptions& opt,
                             void *msg,
                             CkArrMgr *mgr )
// *****************************************************************************
//  Create initial set of array elements based on the unstructured-mesh-aware
//  distribution computed in the constructor
//! \param[in] opt Charm++ array options object containing the number of initial
//!   array elements to be created
//! \param[in] msg Charm++ messsage to use for array element creation
//! \param[in] mgr Array manager to use
//! \author J. Bakosi
// *****************************************************************************
{
  int nelem = *opt.getNumInitial().data(); // number of array elements requested
  if (nelem == 0) return;                  // no initial elements requested

  Assert( static_cast<std::size_t>(nelem) <= m_pe.size(),
          "Number of initial array elements larger than "
          "the UnsMeshMap object holds PEs for" );

  for (int e=0; e<nelem; ++e)
    if (m_pe[ static_cast<std::size_t>(e) ] ==
        static_cast<std::size_t>(CkMyPe()))
      mgr->insertInitial( CkArrayIndex(e), CkCopyMsg(&msg) );

  mgr->doneInserting();
  CkFreeMsg( msg );
}

#include "NoWarning/unsmeshmap.def.h"
