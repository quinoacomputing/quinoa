// *****************************************************************************
/*!
  \file      src/Inciter/Refiner.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh refiner for interfacing the mesh refinement library
  \details   Mesh refiner is a Charm++ chare array and is used to interface the
   mesh refinement object which does not know about parallelization and thus the
   distributed nature of the mesh it operates on, i.e., it operates on mesh
   chunks. Thus it does not do parallel communication and also does not know
   about global vs local IDs. Instead this Charm++ chare array is the one that
   does all parallel computing aspects, i.e., communcation, and using the mesh
   refiner object as a library.
*/
// *****************************************************************************
#ifndef Refiner_h
#define Refiner_h

#include <vector>
#include <unordered_map>

#include "PUPAMR.h"
#include "AMR/mesh_adapter.h"
#include "Inciter/Options/AMRInitial.h"
#include "TaggedTuple.h"
#include "Tags.h"
#include "Callback.h"
#include "UnsMesh.h"
#include "Base/Fields.h"
#include "SchemeBase.h"
#include "DiagCG.h"
#include "ALECG.h"
#include "DG.h"

#include "NoWarning/transporter.decl.h"
#include "NoWarning/refiner.decl.h"

namespace inciter {

//! Mesh refiner for interfacing the mesh refinement library
class Refiner : public CBase_Refiner {

  public:
    //! Constructor
    explicit Refiner( const CProxy_Transporter& transporter,
                      const CProxy_Sorter& sorter,
                      const tk::CProxy_MeshWriter& meshwriter,
                      const Scheme& scheme,
                      const tk::RefinerCallback& cbr,
                      const tk::SorterCallback& cbs,
                      const std::vector< std::size_t >& ginpoel,
                      const tk::UnsMesh::CoordMap& coordmap,
                      const std::map< int, std::vector< std::size_t > >& bface,
                      const std::vector< std::size_t >& triinpoel,
                      const std::map< int, std::vector< std::size_t > >& bnode,
                      int nchare );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit Refiner( CkMigrateMessage* ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ reduction types
    static void registerReducers();

    //! Query Sorter and update local mesh with the reordered one
    void reorder();

    //! Start new step of initial mesh refinement
    void start();

    //! Continue after finishing a refinement step
    void next();

    //! Start mesh refinement (during time stepping, t>0)
    void dtref( const std::map< int, std::vector< std::size_t > >& bface,
                const std::map< int, std::vector< std::size_t > >& bnode,
                const std::vector< std::size_t >& triinpoel );

    //! Receive boundary edges from all PEs (including this one)
    void addBndEdges( CkReductionMsg* msg );

    //! Receive newly added mesh edges and locks on our chare boundary
    void addRefBndEdges( int fromch,
                         const AMR::EdgeData& ed,
                         const std::unordered_set<size_t>& intermediates );

    //! Correct refinement to arrive at conforming mesh across chare boundaries
    void correctref();

    //! Communicate refined edges after a refinement step
    void comExtra();

    //! Decide what to do after a mesh refinement step
    void eval();

    //! Send Refiner proxy to Discretization objects
    void sendProxy();

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_host;
      p | m_sorter;
      p | m_meshwriter;
      p | m_scheme;
      p | m_cbr;
      p | m_cbs;
      p | m_ginpoel;
      p | m_el;
      if (p.isUnpacking()) {
        m_inpoel = std::get< 0 >( m_el );
        m_gid = std::get< 1 >( m_el );
        m_lid = std::get< 2 >( m_el );
      }
      p | m_coordmap;
      p | m_coord;
      p | m_bface;
      p | m_bnode;
      p | m_triinpoel;
      p | m_nchare;
      p | m_initial;
      p | m_initref;
      p | m_refiner;
      p | m_nref;
      p | m_extra;
      p | m_ch;
      p | m_edgedata;
      p | m_edgedataCh;
      p | m_intermediates;
      p | m_bndEdges;
      p | m_msumset;
      p | m_oldTetIdMap;
      p | m_addedNodes;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] r Refiner object reference
    friend void operator|( PUP::er& p, Refiner& r ) { r.pup(p); }
    //@}

  private:
    //! \brief Used to associate a pair of side set id and adjacent tet id to a
    //!   boundary triangle face
    using BndFaces = std::unordered_map< tk::UnsMesh::Face,
                                         std::pair< int, std::size_t >,
                                         tk::UnsMesh::Hash<3>,
                                         tk::UnsMesh::Eq<3> >;

    //! Unique set of edges
    using EdgeSet = std::unordered_set< tk::UnsMesh::Edge,
                                        tk::UnsMesh::Hash<2>,
                                        tk::UnsMesh::Eq<2> >;

    //! Host proxy
    CProxy_Transporter m_host;
    //! Mesh sorter proxy
    CProxy_Sorter m_sorter;
    //! Mesh writer proxy
    tk::CProxy_MeshWriter m_meshwriter;
    //! Discretization scheme
    Scheme m_scheme;
    //! Charm++ callbacks associated to compile-time tags for refiner
    tk::RefinerCallback m_cbr;
    //! Charm++ callbacks associated to compile-time tags for sorter
    tk::SorterCallback m_cbs;
    //! Tetrtahedron element connectivity of our chunk of the mesh (global ids)
    std::vector< std::size_t > m_ginpoel;
    //! Elements of the mesh chunk we operate on
    //! \details The first vector is the element connectivity (local IDs), the
    //!   second vector is the global node IDs of owned elements, while the
    //!   third one is a map of global->local node IDs.
    tk::UnsMesh::Chunk m_el;
    //! Alias to element connectivity with local node IDs in m_el
    std::vector< std::size_t >& m_inpoel = std::get<0>( m_el );
    //! Alias to global node IDs of owned elements in m_el
    std::vector< std::size_t >& m_gid = std::get<1>( m_el );
    //! \brief Alias to local node IDs associated to the global ones of owned
    //!    elements in m_el
    std::unordered_map< std::size_t, std::size_t >& m_lid = std::get<2>( m_el );
    //! Coordinates associated to global node IDs of our mesh chunk
    tk::UnsMesh::CoordMap m_coordmap;
    //! Coordinates of mesh nodes of our chunk of the mesh
    tk::UnsMesh::Coords m_coord;
    //! List of boundary faces associated to side-set IDs
    std::map< int, std::vector< std::size_t > > m_bface;
    //! List of boundary nodes associated to side-set IDs
    std::map< int, std::vector< std::size_t > > m_bnode;
    //! Boundary face-node connectivity
    std::vector< std::size_t > m_triinpoel;
    //! Total number of refiner chares
    int m_nchare;
    //! True if initial AMR, false if during time stepping
    bool m_initial;
    //! Initial mesh refinement type list (in reverse order)
    std::vector< ctr::AMRInitialType > m_initref;
    //! Number of initial mesh refinement steps
    std::size_t m_ninitref;
    //! Mesh refiner (library) object
    AMR::mesh_adapter_t m_refiner;
    //! Counter during distribution of newly added nodes to chare-boundary edges
    std::size_t m_nref;
    //! Number of chare-boundary newly added nodes that need correction
    std::size_t m_extra;
    //! Chares we share at least a single edge with
    std::unordered_set< int > m_ch;
    //! Refinement data associated to edges
    AMR::EdgeData m_edgedata;
    //! Refinement data associated to edges shared with other chares
    std::unordered_map< int, AMR::EdgeData > m_edgedataCh;
    //! Intermediate nodes
    std::unordered_set< size_t> m_intermediates;
    //! Boundary edges associated to chares we share these edges with
    std::unordered_map< int, EdgeSet > m_bndEdges;
    //! \brief Global mesh node IDs bordering the mesh chunk held by fellow
    //!    worker chares associated to their chare IDs
    //! \details msum: mesh chunks surrounding mesh chunks and their neighbor
    //!   points. This is the same data as in Discretization::m_msum, but the
    //!   nodelist is stored as a set.
    std::unordered_map< int, std::unordered_set< std::size_t > > m_msumset;
    //! ...
    std::vector< std::size_t > m_oldTetIdMap;
    //! Newly added mesh nodes (local id) and their parents (local ids)
    std::unordered_map< std::size_t, tk::UnsMesh::Edge > m_addedNodes;

    //! Generate flat coordinate data from coordinate map
    tk::UnsMesh::Coords flatcoord( const tk::UnsMesh::CoordMap& coordmap );

    //! Output mesh to file before a new step of mesh refinement
    void t0ref();

    //! Generate boundary edges and send them to all chares
    void bndEdges();

    //! Refine mesh
    void refine();

    //! Finish initiel mesh refinement
    void endt0ref();

    //! Do uniform mesh refinement
    void uniformRefine();

    //! Do error-based mesh refinement
    void errorRefine();

    //! Do mesh refinement based on user explicitly tagging edges
    void edgelistRefine();

    //! Do mesh refinement based on tagging edges based on end-point coordinates
    void coordRefine();

    //! Do mesh refinement correcting PE-boundary edges
    void correctRefine( const AMR::EdgeData& extra );

    //! ...
    void updateEdgeData();

    //! Aggregate number of extra edges across all chares
    void matched();

    //! Update old mesh after refinement
    void updateMesh();

    //! Update volume mesh after mesh refinement
    void newVolMesh( const std::unordered_set< std::size_t >& old,
                     const std::unordered_set< std::size_t >& ref );

    //! Update boundary data structures after mesh refinement
    void newBndMesh( const std::unordered_set< std::size_t >& old,
                     const std::unordered_set< std::size_t >& ref );

    //! \brief Generate boundary data structure used to update refined
    //!   boundary faces and nodes assigned to side sets
    BndFaces boundary();

    //! Regenerate boundary faces after mesh refinement step
    void updateBndFaces( const std::unordered_set< std::size_t >& old,
                         const std::unordered_set< std::size_t >& ref,
                         const BndFaces& bnd );

    //! Regenerate boundary nodes after mesh refinement step
    void updateBndNodes( const std::unordered_set< std::size_t >& old,
                         const std::unordered_set< std::size_t >& ref,
                         const BndFaces& bnd );

    //! Evaluate initial conditions (IC) at mesh nodes
    tk::Fields nodeinit( std::size_t npoin,
                         const std::pair< std::vector< std::size_t >,
                                          std::vector< std::size_t > >& esup );

    //! Output mesh to file(s)
    void writeMesh( const std::string& basefilename,
                    uint64_t it,
                    tk::real t,
                    CkCallback c );

    //! Functor to call the resize() member function behind SchemeBase::Proxy
    struct Resize : boost::static_visitor<> {
      const std::vector< std::size_t >& Ginpoel;
      const tk::UnsMesh::Chunk& Chunk;
      const tk::UnsMesh::Coords& Coord;
      const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& AddedNodes;
      const std::unordered_map< int, std::vector< std::size_t > >& Msum;
      const std::map< int, std::vector< std::size_t > > Bface;
      const std::map< int, std::vector< std::size_t > > Bnode;
      const std::vector< std::size_t > Triinpoel;
      Resize(
        const std::vector< std::size_t >& ginpoel,
        const tk::UnsMesh::Chunk& chunk,
        const tk::UnsMesh::Coords& coord,
        const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& addednodes,
        const std::unordered_map< int, std::vector< std::size_t > >& msum,
        const std::map< int, std::vector< std::size_t > >& bface,
        const std::map< int, std::vector< std::size_t > >& bnode,
        const std::vector< std::size_t >& triinpoel )
        : Ginpoel(ginpoel), Chunk(chunk), Coord(coord), AddedNodes(addednodes),
          Msum(msum), Bface(bface), Bnode(bnode), Triinpoel(triinpoel) {}
      template< typename P > void operator()( const P& p ) const {
        p.ckLocal()->resize( Ginpoel, Chunk, Coord, AddedNodes, Msum, Bface,
                             Bnode, Triinpoel );
      }
    };
};

} // inciter::

#endif // Refiner_h
