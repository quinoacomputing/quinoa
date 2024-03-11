// *****************************************************************************
/*!
  \file      src/Inciter/Partitioner.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include "ContainerUtil.hpp"
#include "ZoltanInterOp.hpp"
#include "Inciter/InputDeck/New2InputDeck.hpp"
#include "Options/PartitioningAlgorithm.hpp"
#include "DerivedData.hpp"
#include "UnsMesh.hpp"
#include "FaceData.hpp"
#include "Sorter.hpp"
#include "Refiner.hpp"
#include "Callback.hpp"

#include "NoWarning/partitioner.decl.h"

namespace inciter {

extern ctr::New2InputDeck g_newinputdeck;

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
        std::unordered_map< int, std::vector< std::size_t > >,
        // Mesh block ids (value) associated to local tet ids (index)
        std::vector< std::size_t > >;

  public:
    //! Constructor
    Partitioner( std::size_t meshid,
                 const std::string& filename,
                 const tk::PartitionerCallback& cbp,
                 const tk::RefinerCallback& cbr,
                 const tk::SorterCallback& cbs,
                 const CProxy_Transporter& host,
                 const CProxy_Refiner& refiner,
                 const CProxy_Sorter& sorter,
                 const tk::CProxy_MeshWriter& meshwriter,
                 const std::vector< Scheme >& scheme,
                 const std::map< int, std::vector< std::size_t > >& bface,
                 const std::map< int, std::vector< std::size_t > >& faces,
                 const std::map< int, std::vector< std::size_t > >& bnode );

    //! Destructor
    ~Partitioner() override;

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    explicit Partitioner( CkMigrateMessage* m ) : CBase_Partitioner( m ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Partition the computational mesh into a number of chares
    void partition( int nchare );

    //! Receive mesh associated to chares we own after refinement
    void addMesh( int fromnode,
                  const std::unordered_map< int,
                    std::tuple<
                      std::vector< std::size_t >,
                      tk::UnsMesh::CoordMap,
                      std::unordered_map< int, std::vector< std::size_t > >,
                      std::unordered_map< int, std::vector< std::size_t > >,
                      std::vector< std::size_t >
                    > >& chmesh );

    //! Acknowledge received mesh after initial mesh refinement
    void recvMesh();

    //! Optionally start refining the mesh
    void refine();

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \note This is a Charm++ nodegroup, pup() is thus only for
    //!    checkpoint/restart.
    void pup( PUP::er &p ) override {
      p | m_meshid;
      p | m_cbp;
      p | m_cbr;
      p | m_cbs;
      p | m_host;
      p | m_refiner;
      p | m_sorter;
      p | m_meshwriter;
      p | m_scheme;
      p | m_ginpoel;
      p | m_coord;
      p | m_inpoel;
      p | m_lid;
      p | m_elemBlockId;
      p | m_ndist;
      p | m_nchare;
      p | m_nface;
      p | m_nodech;
      p | m_linnodes;
      p | m_chinpoel;
      p | m_chcoordmap;
      p | m_chbface;
      p | m_chtriinpoel;
      p | m_chbnode;
      p | m_chelemblockid;
      p | m_bnodechares;
      p | m_bface;
      p | m_triinpoel;
      p | m_bnode;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i Partitioner object reference
    friend void operator|( PUP::er& p, Partitioner& i ) { i.pup(p); }
    //@}

  private:
    //! Mesh ID
    std::size_t m_meshid;
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
    //! Mesh writer proxy
    tk::CProxy_MeshWriter m_meshwriter;
    //! Discretization schemes (one per mesh)
    std::vector< Scheme > m_scheme;
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
    //! List of elements for each block-id.
    //! \details key: block id, value: set of elements in corresponding block
    std::unordered_map< std::size_t, std::set< std::size_t > > m_elemBlockId;
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
    //! Mesh block ids associated to local tet ids for each chare
    //! \details outer key: chare id, vector index: tet id, value: block id of
    //!   corresponding tet.
    std::unordered_map< int, std::vector< std::size_t > > m_chelemblockid;
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
