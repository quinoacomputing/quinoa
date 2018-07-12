// *****************************************************************************
/*!
  \file      src/Inciter/Partitioner.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ chare partitioner group used to perform mesh partitioning
  \details   Charm++ chare partitioner group used to perform mesh partitioning.
*/
// *****************************************************************************
#ifndef Partitioner_h
#define Partitioner_h

#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>

#include "ExodusIIMeshReader.h"
#include "ContainerUtil.h"
#include "ZoltanInterOp.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Options/PartitioningAlgorithm.h"
#include "Solver.h"
#include "DerivedData.h"
#include "UnsMesh.h"
#include "FaceData.h"
#include "Sorter.h"
#include "Refiner.h"
#include "Callback.h"

#include "NoWarning/partitioner.decl.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! Partitioner Charm++ chare group class
//! \details Instantiations of Partitioner comprise a processor aware Charm++
//!   chare group. When instantiated, a new object is created on each PE and not
//!   more (as opposed to individual chares or chare array object elements). See
//!   also the Charm++ interface file partitioner.ci.
class Partitioner : public CBase_Partitioner {

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-parameter"
    #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-parameter"
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  #elif defined(__INTEL_COMPILER)
    #pragma warning( push )
    #pragma warning( disable: 1478 )
  #endif
  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  Partitioner_SDAG_CODE
  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #elif defined(__INTEL_COMPILER)
    #pragma warning( pop )
  #endif

  public:
    //! Constructor
    Partitioner( const tk::PartitionerCallback& cbp,
                 const tk::RefinerCallback& cbr,
                 const tk::SorterCallback& cbs,
                 const CProxy_Transporter& host,
                 const tk::CProxy_Solver& solver,
                 const CProxy_Refiner& refiner,
                 const CProxy_Sorter& sorter,
                 const Scheme& scheme,
                 const std::map< int, std::vector< std::size_t > >& bface,
                 const std::vector< std::size_t >& triinpoel,
                 const std::map< int, std::vector< std::size_t > >& bnode );

    //! Partition the computational mesh into a number of chares
    void partition( int nchare );

    //! Receive number of uniquely assigned global mesh node IDs from lower PEs
    void offset( int p, std::size_t u );

    //! Request new global node IDs for old node IDs
    void request( int p, const std::unordered_set< std::size_t >& nd );

    // Request new global node IDs for edges
    void request( int p, const tk::UnsMesh::EdgeSet& ed );

    //! Receive new (reordered) global node IDs and coordinates
    void neworder( const std::unordered_map< std::size_t,
           std::tuple< std::size_t, tk::UnsMesh::Coord > >& nodes );

    //! Receive mesh elements and their node coordinates after partitioning
    void addPeMesh( int frompe,
                    const std::vector< std::size_t >& inpoel,
                    const tk::UnsMesh::CoordMap& cm );

    //! Receive mesh associated to chares we own after refinement
    void addMesh( int frompe,
                  const std::unordered_map< int,
                          std::tuple< std::vector< std::size_t >,
                            tk::UnsMesh::CoordMap > >& chmesh );

    //! Acknowledge received mesh after initial mesh refinement
    void recvMesh();

    //! Optionally start refining the mesh
    void refine();

  private:
    //! Charm++ callbacks associated to compile-time tags for partitioner
    tk::PartitionerCallback m_cbp;
    //! Charm++ callbacks associated to compile-time tags for refiner
    tk::RefinerCallback m_cbr;
    //! Charm++ callbacks associated to compile-time tags for sorter
    tk::SorterCallback m_cbs;
    //! Host proxy
    CProxy_Transporter m_host;
    //! Linear system solver proxy
    tk::CProxy_Solver m_solver;
    //! Mesh refiner proxy
    CProxy_Refiner m_refiner;
    //! Mesh sorter proxy
    CProxy_Sorter m_sorter;
    //! Discretization scheme
    Scheme m_scheme;
    //! Element connectivity of this PE's chunk of the mesh (global ids)
    std::vector< std::size_t > m_ginpoel;
    //! Coordinates of mesh nodes of this PE's mesh chunk
    tk::UnsMesh::Coords m_coord;
    //! Element connectivity with local node IDs of this PE's mesh chunk
    std::vector< std::size_t > m_inpoel;
    //! Global node IDs of elements of this PE's mesh chunk
    std::vector< std::size_t > m_gid;
    //! Global->local node IDs of elements of this PE's mesh chunk
    //! \details Key: global node id, value: local node id
    std::unordered_map< std::size_t, std::size_t > m_lid;
    //! Queue of requested node IDs from PEs
    std::vector< std::pair< int, std::unordered_set<std::size_t> > > m_reqNodes;
    //! \brief Starting global mesh node ID for node reordering on this PE
    //!   during mesh node reordering
    std::size_t m_start;
    //! Counter during mesh distribution
    std::size_t m_ndist;
    //! \brief Counter for number of offsets
    //! \details This counts the to-be-received node IDs received while
    //!   computing global mesh node ID offsets for each PE rquired for node
    //!   reordering later
    std::size_t m_noffset;
    //! \brief Counter for number of queries for global mesh node IDs
    //! \details This counts the number of queries received while
    //!   gathering the node IDs that need to be received (instead of uniquely
    //!   assigned) by each PE
    std::size_t m_nquery;
    //! \brief Counter for number of masks of to-be-received global mesh node
    //!   IDs received
    //! \details This counts the to-be-received node ID masks received while
    //!   gathering the node IDs that need to be received (instead of uniquely
    //!   assigned) by each PE
    std::size_t m_nmask;
    //! Coordinates associated to global node IDs of our mesh chunk
    tk::UnsMesh::CoordMap m_coordmap;
    //! Total number of chares across all PEs
    int m_nchare;
    //! Chare IDs (value) associated to global mesh node IDs (key)
    //! \details Multiple chares can contribute to a single node, hence vector
    //!   for map value.
    std::unordered_map< std::size_t, std::vector< int > > m_nodech;
    //! \brief Map associating new node IDs (as in producing contiguous-row-id
    //!   linear system contributions) as map-values to old node IDs (as in
    //!   file) as map-keys
    std::unordered_map< std::size_t, std::size_t > m_linnodes;
    //! Mesh connectivity using global node IDs associated to chares owned
    std::unordered_map< int, std::vector< std::size_t > > m_chinpoel;
    //! Coordinates associated to global node IDs of our mesh chunk for chares
    std::unordered_map< int, tk::UnsMesh::CoordMap > m_chcoordmap;
    //! \brief Maps associating old node IDs to new node IDs (as in producing
    //!   contiguous-row-id linear system contributions) categorized by chares.
    //! \details Maps associating old node IDs (as in file) as map-values to new
    //!   node IDs (as in producing contiguous-row-id linear system
    //!   contributions) as map-keys, associated to chare IDs (outer keys).
    std::unordered_map< int,
      std::unordered_map< std::size_t, std::size_t > > m_chfilenodes;
    //! Communication cost of linear system merging for our PE
    tk::real m_cost;
    //! \brief Map associating a list of chare IDs to old (as in file) global
    //!   mesh node IDs on the chare boundaries
    //! \details Note that a single global mesh node ID can be associated to
    //!   multiple chare IDs as multiple chares can contribute to a single node.
    std::unordered_map< std::size_t, std::vector< int > > m_bnodechares;
    //! \brief Global mesh node IDs associated to chare IDs bordering the mesh
    //!   chunk held by (and associated to) chare IDs this PE owns
    //! \details msum: (M)esh chunks (S)urrounding (M)esh chunks storing mesh
    //!   nodes. Outer map-key: chare IDs this PE owns whose neighbors are
    //!   stored, inner map-key: chare IDs of those chares that hold mesh chunks
    //!   surrounding the outer-key chare's mesh, map-values: global new
    //!   (reordered, as in producing contiguous-row-id linear system
    //!   contributions) mesh node IDs along the border of chares (at which the
    //!   chares will need to communicate) during time stepping.
    std::unordered_map< int,
      std::unordered_map< int, std::unordered_set< std::size_t > > > m_msum;
    //! List of boundary faces associated to side-set IDs
    std::map< int, std::vector< std::size_t > > m_bface;
    //! Boundary face-node connectivity
    std::vector< std::size_t > m_triinpoel;
    //! List of boundary nodes associated to side-set IDs
    std::map< int, std::vector< std::size_t > > m_bnode;

    //! Compute element centroid coordinates
    std::array< std::vector< tk::real >, 3 >
    centroids( const std::vector< std::size_t >& inpoel,
               const tk::UnsMesh::Coords& coord );

    //!  Categorize mesh elements (given by their gobal node IDs) by target
    std::unordered_map< int, std::vector< std::size_t > >
    categorize( const std::vector< std::size_t >& che,
                const std::vector< std::size_t >& inpoel ) const;

    //! Extract coordinates associated to global nodes of a mesh chunk
    tk::UnsMesh::CoordMap coordmap( const std::vector< std::size_t >& inpoel );

    //! Distribute mesh to target PEs after mesh partitioning
    void distribute(
           std::unordered_map< int, std::vector< std::size_t > >&& elems );

    //! Compute chare (partition) distribution
    std::array< int, 2 > distribution( int npart ) const;

    //! Return processing element for chare id
    int pe( int id ) const;

    //! Compute communication cost of linear system merging for our PE
    tk::real cost( std::size_t l, std::size_t u );
};

} // inciter::

#endif // Partitioner_h
