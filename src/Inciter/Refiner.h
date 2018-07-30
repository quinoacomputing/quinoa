// *****************************************************************************
/*!
  \file      src/Inciter/Refiner.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
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

#include "AMR/mesh_adapter.h"
#include "Inciter/Options/AMRInitial.h"
#include "TaggedTuple.h"
#include "Tags.h"
#include "Callback.h"
#include "UnsMesh.h"
#include "Base/Fields.h"

#include "NoWarning/transporter.decl.h"
#include "NoWarning/refiner.decl.h"

namespace inciter {

//! Mesh refiner for interfacing the mesh refinement library
class Refiner : public CBase_Refiner {

  public:
    //! Constructor
    explicit Refiner( const CProxy_Transporter& transporter,
                      const CProxy_Sorter& sorter,
                      const tk::CProxy_Solver& solver,
                      const Scheme& scheme,
                      const tk::RefinerCallback& cbr,
                      const tk::SorterCallback& cbs,
                      const std::vector< std::size_t >& ginpoel,
                      const tk::UnsMesh::CoordMap& coordmap,
                      const std::map< int, std::vector< std::size_t > >& bface,
                      const std::vector< std::size_t >& triinpoel,
                      const std::map< int, std::vector< std::size_t > >& bnode,
                      int nchare );

    //! Configure Charm++ reduction types
    static void registerReducers();

    //! Receive boundary edges from all PEs (including this one)
    void addBndEdges( CkReductionMsg* msg );

    //! Receive newly added mesh node IDs on our chare boundary
    void addRefBndEdges( int fromch, const tk::UnsMesh::EdgeNodeCoord& ed );

    //! Acknowledge received newly added nodes shared with other chares
    void recvRefBndEdges();

    //! Correct refinement to arrive at conforming mesh across chare boundaries
    void correctref();

    //! Decide wether to continue with another step of mesh refinement
    void nextref();

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_Refiner::pup(p);
      p | m_host;
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
      p | m_triinpoel;
      p | m_bnode;
      p | m_nchare;
      p | m_initref;
      //p | m_refiner;
      p | m_nref;
      p | m_extra;
      p | m_ch;
      p | m_edgenode;
      p | m_edgenodeCh;
      p | m_bndEdges;
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
    //! Linear soilver proxy
    tk::CProxy_Solver m_solver;
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
    std::tuple< std::vector< std::size_t >,
                std::vector< std::size_t >,
                std::unordered_map< std::size_t, std::size_t > > m_el;
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
    //! Boundary face-node connectivity
    std::vector< std::size_t > m_triinpoel;
    //! List of boundary nodes associated to side-set IDs
    std::map< int, std::vector< std::size_t > > m_bnode;
    //! Total number of refiner chares
    int m_nchare;
    //! Initial mesh refinement type list (in reverse order)
    std::vector< ctr::AMRInitialType > m_initref;
    //! Mesh refiner (library) object
    AMR::mesh_adapter_t m_refiner;
    //! Counter during distribution of newly added nodes to chare-boundary edges
    std::size_t m_nref;
    //! Number of chare-boundary newly added nodes that need correction
    std::size_t m_extra;
    //! Chares we share at least a single edge with
    std::unordered_set< int > m_ch;
    //! Map associating global IDs and coordinates of a node added to an edge
    tk::UnsMesh::EdgeNodeCoord m_edgenode;
    //! \brief Map associating global IDs and coordinates of a node added to an
    //!   edge associated to another chare the edge is shared with
    std::unordered_map< int, tk::UnsMesh::EdgeNodeCoord > m_edgenodeCh;
    //! Boundary edges associated to chares we share these edges with
    std::unordered_map< int, EdgeSet > m_bndEdges;

    //! Generate flat coordinate data from coordinate map
    tk::UnsMesh::Coords flatcoord( const tk::UnsMesh::CoordMap& coordmap );

    //! Prepare for next step of mesh refinement
    void start();

    //! Generate boundary edges and send them to all chares
    void bndEdges();

    //! Refine mesh
    void refine();

    //! Communicate refined edges after a refinement step
    void comExtra();

    //! Finish initiel mesh refinement
    void finish();

    //! Do uniform mesh refinement
    void uniformRefine();

    //! Do error-based mesh refinement
    void errorRefine();

    //! Do mesh refinement based on user explicitly tagging edges
    void userRefine();

    //! Do mesh refinement based on tagging edges based on end-point coordinates
    void coordRefine();

    //! Do mesh refinement correcting PE-boundary edges
    void correctRefine( const EdgeSet& extra );

    //! Update mesh after refinement
    void updateMesh();

    //! Update volume mesh after mesh refinement
    void updateVolMesh( const std::unordered_set< std::size_t >& old,
                        const std::unordered_set< std::size_t >& ref );

    //! Update boundary data structures after mesh refinement
    void updateBndMesh( const std::unordered_set< std::size_t >& old,
                        const std::unordered_set< std::size_t >& ref );

    //! Generate boundary data structures after mesh refinement
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
};

} // inciter::

#endif // Refiner_h
