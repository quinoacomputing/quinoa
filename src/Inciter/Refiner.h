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

    //! Refine mesh
    void refine();

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

    //! Get refinement field data in mesh cells
    std::tuple< std::vector< std::string >,
                std::vector< std::vector< tk::real > > >
    refinementFields() const;

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
      p | m_oldTets;
      p | m_addedNodes;
      p | m_addedTets;
      p | m_prevnTets;
      p | m_grand;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] r Refiner object reference
    friend void operator|( PUP::er& p, Refiner& r ) { r.pup(p); }
    //@}

  private:
    //! Used to associate tet id (value) to a boundary triangle face (key)
    using BndFaceTets =
      std::unordered_map< tk::UnsMesh::Face, std::size_t,
                          tk::UnsMesh::Hash<3>, tk::UnsMesh::Eq<3> >;

    //! Used to associate side set ids (value) to a boundary triangle face (key)
    using BndFaceSides =
      std::unordered_map< tk::UnsMesh::Face, std::unordered_set< int >,
                          tk::UnsMesh::Hash<3>, tk::UnsMesh::Eq<3> >;

    //! Used to associate a unique set of faces (3 node ids) to side sets
    using BndFaceSets = std::unordered_map< int, tk::UnsMesh::FaceSet >;

    //! Used to associate a unique set of nodes to side sets
    using BndNodeSets =
      std::unordered_map< int, std::unordered_set< std::size_t > >;

    //! Boundary face data bundle
    using BndFaceData =
      std::tuple< BndFaceTets, BndFaceSides, BndFaceSets, BndNodeSets >;

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
    std::unordered_map< int, tk::UnsMesh::EdgeSet > m_bndEdges;
    //! \brief Global mesh node IDs bordering the mesh chunk held by fellow
    //!    worker chares associated to their chare IDs
    //! \details msum: mesh chunks surrounding mesh chunks and their neighbor
    //!   points. This is the same data as in Discretization::m_msum, but the
    //!   nodelist is stored as a set.
    std::unordered_map< int, std::unordered_set< std::size_t > > m_msumset;
    //! ...
    std::vector< std::size_t > m_oldTetIdMap;
    //! Local tetrahedron IDs before refinement step
    std::unordered_set< std::size_t > m_oldTets;
    //! Newly added mesh nodes (local id) and their parents (local ids)
    std::unordered_map< std::size_t, tk::UnsMesh::Edge > m_addedNodes;
    //! Newly added mesh cells (local id) and their parent (local ids)
    std::unordered_map< std::size_t, std::size_t > m_addedTets;
    //! Number of tetrahedra in the mesh before refinement
    std::size_t m_prevnTets;
    //! Unique nodes of the mesh before a refinement step with local ids
    std::unordered_set< std::size_t > m_grand;

    //! Generate flat coordinate data from coordinate map
    tk::UnsMesh::Coords flatcoord( const tk::UnsMesh::CoordMap& coordmap );

    //! Output mesh to file before a new step of mesh refinement
    void t0ref();

    //! Generate boundary edges and send them to all chares
    void bndEdges();

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

    //! \brief Generate boundary data structures used to update refined
    //!   boundary faces and nodes of side sets
    BndFaceData boundary();

    //! Regenerate boundary faces after mesh refinement step
    void updateBndFaces( const std::unordered_set< std::size_t >& old,
                         const std::unordered_set< std::size_t >& ref,
                         const BndFaceData& bnd );

    //! Regenerate boundary nodes after mesh refinement step
    void updateBndNodes( const std::unordered_set< std::size_t >& old,
                         const std::unordered_set< std::size_t >& ref,
                         const BndFaceData& bnd );

    //! Evaluate initial conditions (IC) at mesh nodes
    tk::Fields nodeinit( std::size_t npoin,
                         const std::pair< std::vector< std::size_t >,
                                          std::vector< std::size_t > >& esup );

    //! Output mesh to file(s)
    void writeMesh( const std::string& basefilename,
                    uint64_t it,
                    tk::real t,
                    CkCallback c ) const;

    //! Find parents of a list of nodes
    //! \tparam Container Container whose value_type is std::size_t
    //! \param[in] oldmesh Unique nodes of the old mesh using local ids
    //! \param[in] nodes A container of node ids. The ids must have a type of
    //!   std::size_t.
    //! \return A unique set of nodes ids of the parents of the input nodes
    //! \details Parameter 'oldmesh' defines the old mesh, i.e., the mesh that
    //!   is considered as the parent mesh in which if nodes are found during
    //!   the search, they are considered parents of themselves. In other words,
    //!   parent lookup in the h-refinement hierarchy stops at the mesh whose
    //!   unique node ids are passed in in oldmesh.
    template< class Container >
    std::unordered_set< std::size_t >
    parents( const std::unordered_set< std::size_t >& oldmesh,
             const Container& nodes )
    {
      static_assert(
        std::is_same< typename Container::value_type, std::size_t >::value,
        "Container::value_type must be of type std::size_t" );
      std::unordered_set< std::size_t > s;
      for (auto n : nodes) {
        if (oldmesh.find(n) != end(oldmesh)) {
          s.insert( n );
        } else {
          auto p = m_refiner.node_connectivity.get( n );
          s.insert( begin(p), end(p) );
        }
      }
      // Ensure all parent nodes are in the old mesh
      Assert( std::all_of( begin(s), end(s), [ &oldmesh ]( std::size_t j ){
                return oldmesh.find(j) != end(oldmesh); } ),
              "Old face nodes not in old mesh" );
      return s;
    }

    //! Return a set of side set ids for a primitive
    //! \tparam BndSets Type of map of sets we search for the primitive
    //! \tparam Primitive The primitive (face/node) we search for
    //! \note that BndSets::mapped_type == Primitive must be true
    //! \param[in] sets Map of sets we search in
    //! \param[in] p Primitive (face or node) we search for
    //! \return A unique set of side set ids in which the primitive (face/node)
    //!    is found or an empty set if the primitive was not found.
    template< class BndSets, class Primitive >
    std::unordered_set< int >
    sides( const BndSets& sets, const Primitive& p ) {
      static_assert( std::is_same< typename BndSets::mapped_type::value_type,
        Primitive >::value, "Type of primitive (face/node) in map of sets must "
        "be the same as the type of primitive (face/node) that is searched" );
      std::unordered_set< int > ss;
      for (const auto& s : sets)
        if (s.second.find(p) != end(s.second))
          ss.insert( s.first );
      return ss;
    }

    //! \brief Function class to call the resizeAfterREfined() member function
    //!   behind SchemeBase::Proxy
    struct ResizeAfterRefined : boost::static_visitor<> {
      const std::vector< std::size_t >& Ginpoel;
      const tk::UnsMesh::Chunk& Chunk;
      const tk::UnsMesh::Coords& Coord;
      const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& AddedNodes;
      const std::unordered_map< std::size_t, std::size_t >& AddedTets;
      const std::unordered_map< int, std::vector< std::size_t > >& Msum;
      const std::map< int, std::vector< std::size_t > > Bface;
      const std::map< int, std::vector< std::size_t > > Bnode;
      const std::vector< std::size_t > Triinpoel;
      ResizeAfterRefined(
        const std::vector< std::size_t >& ginpoel,
        const tk::UnsMesh::Chunk& chunk,
        const tk::UnsMesh::Coords& coord,
        const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& addednodes,
        const std::unordered_map< std::size_t, std::size_t >& addedtets,
        const std::unordered_map< int, std::vector< std::size_t > >& msum,
        const std::map< int, std::vector< std::size_t > >& bface,
        const std::map< int, std::vector< std::size_t > >& bnode,
        const std::vector< std::size_t >& triinpoel )
        : Ginpoel(ginpoel), Chunk(chunk), Coord(coord), AddedNodes(addednodes),
          AddedTets(addedtets), Msum(msum), Bface(bface), Bnode(bnode),
          Triinpoel(triinpoel) {}
      template< typename P > void operator()( const P& p ) const {
        p.ckLocal()->resizeAfterRefined( Ginpoel, Chunk, Coord, AddedNodes,
          AddedTets, Msum, Bface, Bnode, Triinpoel );
      }
    };

    //! \brief Function class to call the solution() member function
    //!   behind SchemeBase::Proxy
    struct solution : boost::static_visitor< tk::Fields > {
      template< typename P > const tk::Fields& operator()( const P& p ) const {
        return p.ckLocal()->solution();
      }
    };
};

} // inciter::

#endif // Refiner_h
