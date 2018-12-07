// *****************************************************************************
/*!
  \file      src/Inciter/Partitioner.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ chare partitioner nodegroup used to perform mesh
             partitioning
  \details   Charm++ chare partitioner nodegroup used to perform mesh read and
             partitioning, one worker per compute node.
*/
// *****************************************************************************
#ifndef Partitioner_h
#define Partitioner_h

#include <array>
#include <stddef.h>

#include "ContainerUtil.h"
#include "ZoltanInterOp.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Options/PartitioningAlgorithm.h"
#include "DerivedData.h"
#include "UnsMesh.h"
#include "FaceData.h"
#include "Sorter.h"
#include "Refiner.h"
#include "Callback.h"

#include "NoWarning/partitioner.decl.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! Partitioner Charm++ chare nodegroup class
//! \details Instantiations of Partitioner comprise a processor aware Charm++
//!   chare node group. When instantiated, a new object is created on each
//!   compute node and not more (as opposed to individual chares or chare array
//!   object elements). See also the Charm++ interface file partitioner.ci.
class Partitioner : public CBase_Partitioner {

  private:
    //! \brief Mesh data used for categorizing mesh chunks assigned to chares
    //!    after mesh partitioning and before mesh distribution across chares
    using MeshData =
      std::tuple<
        // Tetrahedron (domain element) connectivity
        std::vector< std::size_t >,
        // Boundary face connectivity for each side set
        std::unordered_map< int, std::vector< std::size_t > >,
        // Boundary node lists for each side set
        std::unordered_map< int, std::vector< std::size_t > > >;

  public:
    //! Constructor
    Partitioner( const tk::PartitionerCallback& cbp,
                 const tk::RefinerCallback& cbr,
                 const tk::SorterCallback& cbs,
                 const CProxy_Transporter& host,
                 const CProxy_Refiner& refiner,
                 const CProxy_Sorter& sorter,
                 const Scheme& scheme,
                 const std::map< int, std::vector< std::size_t > >& belem,
                 const std::map< int, std::vector< std::size_t > >& faces,
                 const std::map< int, std::vector< std::size_t > >& bnode );

    //! Turn off automatic load balancing
    void lboff();

    //! Partition the computational mesh into a number of chares
    void partition( int nchare );

    //! Receive mesh associated to chares we own after refinement
    void addMesh( int fromnode,
                  const std::unordered_map< int,
                    std::tuple<
                      std::vector< std::size_t >,
                      tk::UnsMesh::CoordMap,
                      std::unordered_map< int, std::vector< std::size_t > >,
                      std::unordered_map< int, std::vector< std::size_t > >
                    > >& chmesh );

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
    //! Mesh refiner proxy
    CProxy_Refiner m_refiner;
    //! Mesh sorter proxy
    CProxy_Sorter m_sorter;
    //! Discretization scheme
    Scheme m_scheme;
    //! Element connectivity of this compute node's mesh chunk (global ids)
    std::vector< std::size_t > m_ginpoel;
    //! Coordinates of mesh nodes of this compute node's mesh chunk
    tk::UnsMesh::Coords m_coord;
    //! \brief Element connectivity with local node IDs of this compute node's
    //!   mesh chunk
    std::vector< std::size_t > m_inpoel;
    //! Global->local node IDs of elements of this compute node's mesh chunk
    //! \details Key: global node id, value: local node id
    std::unordered_map< std::size_t, std::size_t > m_lid;
    //! Counter during mesh distribution
    std::size_t m_ndist;
    //! Total number of chares across all compute nodes
    int m_nchare;
    //! Counters (for each chare owned) for assigning face ids in parallel
    std::unordered_map< int, std::size_t > m_nface;
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
    //! Side set id + boundary face id for each chare
    std::unordered_map< int,
      std::map< int, std::vector< std::size_t > > > m_chbface;
    //! Boundary face connectivity for each chare
    std::map< int, std::vector< std::size_t > > m_chtriinpoel;
    //! Side set id + boundary nodes for each chare
    std::unordered_map< int,
      std::map< int, std::vector< std::size_t > > > m_chbnode;
    //! \brief Map associating a list of chare IDs to old (as in file) global
    //!   mesh node IDs on the chare boundaries
    //! \details Note that a single global mesh node ID can be associated to
    //!   multiple chare IDs as multiple chares can contribute to a single node.
    std::unordered_map< std::size_t, std::vector< int > > m_bnodechares;
    //! Boundary face IDs associated associated to side set IDs
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
    std::unordered_map< int, MeshData >
    categorize( const std::vector< std::size_t >& che ) const;

    //! Extract coordinates associated to global nodes of a mesh chunk
    tk::UnsMesh::CoordMap coordmap( const std::vector< std::size_t >& inpoel );

    //! Distribute mesh to target compute nodes after mesh partitioning
    void distribute( std::unordered_map< int, MeshData >&& mesh );

    //! Compute chare (partition) distribution across compute nodes
    std::array< int, 2 > distribution( int npart ) const;

    //! Return nodegroup id for chare id
    int node( int id ) const;

    //! Keep only those nodes for side sets that reside on this compute node
    void ownBndNodes(
      const std::unordered_map< std::size_t, std::size_t >& lid,
      std::map< int, std::vector< std::size_t > >& bnode );
};

} // inciter::

#endif // Partitioner_h
