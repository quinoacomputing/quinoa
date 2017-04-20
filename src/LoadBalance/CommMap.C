// *****************************************************************************
/*!
  \file      src/LoadBalance/CommMap.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Calculation of communication maps for unstructured meshes
  \details   Calculation of communication maps for unstructured meshes.
*/
// *****************************************************************************

#include <vector>
#include <set>
#include <map>
#include <iosfwd>
#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>

#include "CommMap.h"
#include "DerivedData.h"
#include "Exception.h"
#include "UnsMesh.h"

namespace tk {

std::vector< std::map< std::size_t, std::vector< std::size_t > > >
poinCommMaps( std::size_t graphsize,
              const std::vector< std::size_t >& chp,
              const std::vector< std::size_t >& tetinpoel,
              std::size_t nchare,
              std::string&& toofine )
// *****************************************************************************
//  Compute point-based communication maps
//! \param[in] graphsize Size of unstructured mesh graph object
//! \param[in] chp Array of chare ownership IDs mapping graph points to
//!   concurrent async chares
//! \param[in] tetinpoel Tetrahedra element connectivity
//! \param[in] nchare Number of work units (Charm++ chares)
//! \param[in] toofine Error message to print triggered for too large
//!   overdecomposition
//! \return Point-based communication map for all chares
//! \details This is a _point-based_ export map, because it stores the global
//!   ids of the mesh points that chares need to export to fellow chares
//!   computed based on which chare owns a mesh point. This is for algorithms
//!   that work by computing data on the mesh by looping over mesh points and
//!   their surrounding points, e.g., edge-based algorithms. The communication
//!   map computed here is stored in a vector of maps associating a vector of
//!   global mesh point ids to chare ids that the points need to be exported to.
//!   Using the map produced here amounts to each chare taking its own map
//!   indexing the outermost vector with its chare index. The map computed here
//!   is an an export map from which the import map can be computed if needed.
//!
//!   In the MPI paradigm, these chare export maps correspond to the export
//!   lists, i.e., lists of global ids exported by a given rank to a set of
//!   receiver ranks and their associated mesh points sent at which data are to
//!   be sent (exported). This is as opposed to import maps which are lists of
//!   global ids imported by a given rank to a set of sender ranks and their
//!   associated mesh points sent at which data are to be received (imported).
//!   For example, in Zoltan, this export/import roughly corresponds to the
//!   "exported" and "imported" mesh node ids after partitioning. Actually,
//!   Zoltan_LB_Partition() already returns this information. However, if the
//!   partitioning with Zoltan is done using less MPI ranks than the number of
//!   desired mesh partitions, i.e., overdecomposition (as is the case here),
//!   there are more mesh partitions than MPI ranks and thus the arrays returned
//!   by Zoltan are not sufficient to determine the export and import maps for
//!   all the chares. Thus we compute the export mapping here.
//! \note This function is only supposed to operate on MPI rank 0.
//! \note This function does not and should not modify global-scope data.
//! \author J. Bakosi
// *****************************************************************************
{
  // Map to associate a chare id to a map of receiver chare ids associated to
  // unique global point ids sent (export map)
  std::map< std::size_t,
            std::map< std::size_t, std::set< std::size_t > > > comm;

  // Generate points surrounding points
  auto psup = tk::genPsup( tetinpoel, 4, tk::genEsup(tetinpoel,4) );

  // Construct point-based export maps
  for (std::size_t p=0; p<graphsize; ++p)  // for all mesh points
    for (auto i=psup.second[p]+1; i<=psup.second[p+1]; ++i) {
      auto q = psup.first[i];
      if (chp[p] != chp[q])   // if the point-colors differ, store global id
        comm[ chp[p] ][ chp[q] ].insert( p );
    }

  // This check should always be done, as it can result from incorrect user
  // input compared to the mesh size and not due to programmer error.
  ErrChk( comm.size() == nchare, std::move(toofine) );

  #ifndef NDEBUG
  std::size_t c = 0;
  for (const auto& e : comm)
    Assert( e.first == c++,
            "Export/import maps should not be missing for chare id " +
            std::to_string(c-1) );
  #endif

  // Construct final product: a vector of export maps associating receiver
  // chare ids to unique communicated global point ids for all chare ids, and
  // store it in global scope so that the Charm++ chares can access it
  std::vector< std::map< std::size_t, std::vector< std::size_t > > >
    pcomm( nchare );
  for (const auto& e : comm) {
    for (const auto& x : e.second)
      for (auto p : x.second)
        pcomm[ e.first ][ x.first ].push_back( p );
  }

  Assert( pcomm.size() == nchare,
          "Number of export/import maps must equal the number of chares" );

//   std::size_t h = 0;
//   for (const auto& m : pcomm) {
//     std::cout << h++ << " p-> ";
//     for (const auto& x : m) {
//       std::cout << x.first << ": ";
//       for (auto p : x.second)
//         std::cout << p << " ";
//     }
//     std::cout << '\n';
//   }

  return pcomm;
}

std::vector< std::map< std::size_t, std::vector< std::size_t > > >
elemCommMaps(
  const std::vector< std::size_t >& chp,
  const std::vector< std::size_t >& tetinpoel,
  const std::vector< std::vector< std::vector< std::size_t > > >& element,
  std::size_t nchare )
// *****************************************************************************
//! Compute element-based communication maps
//! \param[in] chp Array of chare ownership IDs mapping graph points to
//!   concurrent async chares
//! \param[in] tetinpoel Tetrahedra element connectivity
//! \param[in] element Global mesh element ids owned by each chare distributed
//!   to PEs
//! \param[in] nchare Number of work units (Charm++ chares)
//! \return Element-based communication map for all chares
//! \details This is an _element-based_ export map, because it stores the global
//!   ids of the mesh points that chares need to export to fellow chares
//!   computed based on which chare owns an element the mesh point is a vertex
//!   of. This is for algorithms that work by computing data on the mesh by
//!   looping over mesh elements, e.g., element-based algorithms. The
//!   communication map computed here is stored in a vector of maps associating
//!   a vector of global mesh point ids to chare ids that the points need to be
//!   exported to. Using the map produced here amounts to each chare taking its
//!   own map indexing the outermost vector with its chare index. The map
//!   computed here is an an export map from which the import map can be
//!   computed if needed.
//!
//!   In the MPI paradigm, these chare export maps correspond to the export
//!   lists, i.e., lists of global ids exported by a given rank to a set of
//!   receiver ranks and their associated mesh points sent at which data are to
//!   be sent (exported). This is as opposed to import maps which are lists of
//!   global ids imported by a given rank to a set of sender ranks and their
//!   associated mesh points sent at which data are to be received (imported).
//!   For example, in Zoltan, this export/import roughly corresponds to the
//!   "exported" and "imported" mesh node ids after partitioning. Actually,
//!   Zoltan_LB_Partition() already returns this information. However, if the
//!   partitioning with Zoltan is done using less MPI ranks than the number of
//!   desired mesh partitions, i.e., overdecomposition (as is the case here),
//!   there are more mesh partitions than MPI ranks and thus the arrays returned
//!   by Zoltan are not sufficient to determine the export and import maps for
//!   all the chares. Thus we compute the export mapping here.
//! \note This function is only supposed to operate on MPI rank 0.
//! \note This function does not and should not modify global-scope data.
//! \author J. Bakosi
// *****************************************************************************
{
  // Map to associate a chare id to a map of receiver chare ids associated to
  // unique global point ids sent (export map)
  std::map< std::size_t,
            std::map< std::size_t, std::set< std::size_t > > > comm;

  // Construct element-based export maps
  for (const auto& pel : element)
    for (std::size_t c=0; c<pel.size(); ++c)
      for (auto e : pel[c])
        for (std::size_t n=0; n<4; ++n) {
          auto p = tetinpoel[e*4+n];
          if (chp[p] != c)      // if the point-colors differ, store global id
            comm[ c ][ chp[p] ].insert( p );
        }

//   for (const auto& c : comm) {
//     std::cout << c.first << " -> ";
//     for (const auto& t : c.second) {
//       std::cout << t.first << ": ";
//       for (auto p : t.second)
//         std::cout << p << " ";
//       std::cout << ", ";
//     }
//   std::cout << '\n';
//   }

  // Construct final product: a vector of export maps associating receiver
  // chare ids to unique communicated global point ids for all chare ids, and
  // store it in global scope so that the Charm++ chares can access it
  std::vector< std::map< std::size_t, std::vector< std::size_t > > >
    ecomm( nchare );
  for (const auto& e : comm)
    for (const auto& x : e.second)
      for (auto p : x.second)
        ecomm[ e.first ][ x.first ].push_back( p );

//   std::size_t h = 0;
//   for (const auto& m : ecomm) {
//     std::cout << h++ << " e-> ";
//     for (const auto& x : m) {
//       std::cout << x.first << ": ";
//       for (auto p : x.second)
//         std::cout << p << " ";
//     }
//     std::cout << '\n';
//   }

  return ecomm;
}

} // tk::
