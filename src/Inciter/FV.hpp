// *****************************************************************************
/*!
  \file      src/Inciter/FV.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     FV advances a system of PDEs with the finite volume scheme
  \details   FV advances a system of partial differential equations (PDEs) using
    the finite volume (FV) spatial discretization (on tetrahedron elements).

    There are a potentially large number of FV Charm++ chares created by
    Transporter. Each FV gets a chunk of the full load (part of the mesh)
    and does the same: initializes and advances a number of PDE systems in time.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/fv.ci.
*/
// *****************************************************************************
#ifndef FV_h
#define FV_h

#include <array>
#include <set>
#include <unordered_set>
#include <unordered_map>

#include "DerivedData.hpp"
#include "FaceData.hpp"
#include "ElemDiagnostics.hpp"

#include "NoWarning/fv.decl.h"

namespace inciter {

//! FV Charm++ chare array used to advance PDEs in time with FV
class FV : public CBase_FV {

  public:
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
    FV_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit FV( const CProxy_Discretization& disc,
                 const std::map< int, std::vector< std::size_t > >& bface,
                 const std::map< int, std::vector< std::size_t > >& /* bnode */,
                 const std::vector< std::size_t >& triinpoel );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit FV( CkMigrateMessage* msg ) : CBase_FV( msg ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Return from migration
    void ResumeFromSync() override;

    //! Start sizing communication buffers and setting up ghost data
    void resizeComm();

    //! Receive unique set of faces we potentially share with/from another chare
    void comfac( int fromch, const tk::UnsMesh::FaceSet& infaces );

    //! Receive ghost data on chare boundaries from fellow chare
    void comGhost( int fromch, const GhostData& ghost );

    //! Receive requests for ghost data
    void reqGhost();

    //! Send all of our ghost data to fellow chares
    void sendGhost();

    //! Setup node-neighborhood (esup)
    void nodeNeighSetup();

    //! Receive element-surr-points data on chare boundaries from fellow chare
    void comEsup( int fromch,
      const std::unordered_map< std::size_t, std::vector< std::size_t > >&
        bndEsup,
      const std::unordered_map< std::size_t, std::vector< tk::real > >&
        nodeBndCells );

    //! Configure Charm++ reduction types for concatenating BC nodelists
    static void registerReducers();

    //! Setup: query boundary conditions, output mesh, etc.
    void setup();

    //! Receive total box IC volume and set conditions in box
    void box( tk::real v );

    // Evaluate whether to do load balancing
    void evalLB( int nrestart );

    //! Start time stepping
    void start();

    //! Continue to next time step
    void next();

    //! Receive chare-boundary limiter function data from neighboring chares
    void comlim( int fromch,
                 const std::vector< std::size_t >& tetid,
                 const std::vector< std::vector< tk::real > >& u,
                 const std::vector< std::vector< tk::real > >& prim );

    //! Receive chare-boundary ghost data from neighboring chares
    void comsol( int fromch,
                 const std::vector< std::size_t >& tetid,
                 const std::vector< std::vector< tk::real > >& u,
                 const std::vector< std::vector< tk::real > >& prim );

    //! \brief Receive nodal solution (ofor field output) contributions from
    //!   neighboring chares
    void comnodeout( const std::vector< std::size_t >& gid,
                     const std::vector< std::size_t >& nesup,
                     const std::vector< std::vector< tk::real > >& L );

    //! Optionally refine/derefine mesh
    void refine( const std::vector< tk::real >& l2res );

    //! Receive new mesh from Refiner
    void resizePostAMR(
      const std::vector< std::size_t >& /* ginpoel */,
      const tk::UnsMesh::Chunk& chunk,
      const tk::UnsMesh::Coords& coord,
      const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& /* addedNodes */,
      const std::unordered_map< std::size_t, std::size_t >& addedTets,
      const std::set< std::size_t >& /*removedNodes*/,
      const tk::NodeCommMap& nodeCommMap,
      const std::map< int, std::vector< std::size_t > >& bface,
      const std::map< int, std::vector< std::size_t > >& /* bnode */,
      const std::vector< std::size_t >& triinpoel );

    //! Extract field output going to file
    void extractFieldOutput(
      const std::vector< std::size_t >& /* ginpoel */,
      const tk::UnsMesh::Chunk& chunk,
      const tk::UnsMesh::Coords& coord,
      const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& /* addedNodes */,
      const std::unordered_map< std::size_t, std::size_t >& addedTets,
      const tk::NodeCommMap& nodeCommMap,
      const std::map< int, std::vector< std::size_t > >& bface,
      const std::map< int, std::vector< std::size_t > >& /* bnode */,
      const std::vector< std::size_t >& triinpoel,
      CkCallback c );

    //! Const-ref access to current solution
    //! \return Const-ref to current solution
    const tk::Fields& solution() const { return m_u; }

    //! Compute left hand side
    void lhs();

    //! Unused in FV
    void resized() {}

    //! Compute right hand side and solve system
    void solve( tk::real newdt );

    //! Evaluate whether to continue with next time step
    void step();

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_disc;
      p | m_ncomfac;
      p | m_nadj;
      p | m_ncomEsup;
      p | m_nsol;
      p | m_ninitsol;
      p | m_nlim;
      p | m_nnod;
      p | m_inpoel;
      p | m_coord;
      p | m_fd;
      p | m_u;
      p | m_un;
      p | m_p;
      p | m_geoFace;
      p | m_geoElem;
      p | m_lhs;
      p | m_rhs;
      p | m_nfac;
      p | m_nunk;
      p | m_npoin;
      p | m_ipface;
      p | m_bndFace;
      p | m_ghostData;
      p | m_sendGhost;
      p | m_ghostReq;
      p | m_ghost;
      p | m_exptGhost;
      p | m_recvGhost;
      p | m_diag;
      p | m_stage;
      p | m_bid;
      p | m_uc;
      p | m_pc;
      p | m_initial;
      p | m_expChBndFace;
      p | m_infaces;
      p | m_esup;
      p | m_esupc;
      p | m_elemfields;
      p | m_nodefields;
      p | m_nodefieldsc;
      p | m_outmesh;
      p | m_boxelems;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i FV object reference
    friend void operator|( PUP::er& p, FV& i ) { i.pup(p); }
    //@}

  private:
    //! Local face & tet IDs associated to 3 global node IDs
    //! \details This map stores tetrahedron cell faces (map key) and their
    //!   associated local face ID and inner local tet id adjacent to the face
    //!   (map value). A face is given by 3 global node IDs.
    using FaceMap =
      std::unordered_map< tk::UnsMesh::Face,  // 3 global node IDs
                          std::array< std::size_t, 2 >, // local face & tet ID
                          tk::UnsMesh::Hash<3>,
                          tk::UnsMesh::Eq<3> >;

    //! Storage type for refined mesh used for field output
    struct OutMesh {
      //! Element connectivity, local->global node ids, global->local nodes ids
      tk::UnsMesh::Chunk chunk;
      //! Node coordinates
      tk::UnsMesh::Coords coord;
      //! Triangle element connectivity
      std::vector< std::size_t > triinpoel;
      //! Boundary-face connectivity
      std::map< int, std::vector< std::size_t > > bface;
      //! Node communinaction map
      tk::NodeCommMap nodeCommMap;
      //! \brief Pack/Unpack serialize member function
      //! \param[in,out] p Charm++'s PUP::er serializer object reference
      void pup( PUP::er& p ) {
        p|chunk; p|coord; p|triinpoel; p|bface; p|nodeCommMap;
      }
      //! Destroyer
      void destroy() {
        tk::destroy( std::get<0>(chunk) );
        tk::destroy( std::get<1>(chunk) );
        tk::destroy( std::get<2>(chunk) );
        tk::destroy( coord[0] );
        tk::destroy( coord[1] );
        tk::destroy( coord[2] );
        tk::destroy( triinpoel );
        tk::destroy( bface );
        tk::destroy( nodeCommMap );
      }
    };

    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Counter for face adjacency communication map
    std::size_t m_ncomfac;
    //! Counter signaling that all ghost data have been received
    std::size_t m_nadj;
    //! Counter for element-surr-node adjacency communication map
    std::size_t m_ncomEsup;
    //! Counter signaling that we have received all our solution ghost data
    std::size_t m_nsol;
    //! \brief Counter signaling that we have received all our solution ghost
    //!    data during setup
    std::size_t m_ninitsol;
    //! \brief Counter signaling that we have received all our limiter function
    //!   ghost data
    std::size_t m_nlim;
    //! \brief Counter signaling that we have received all our node solution
    //!   contributions
    std::size_t m_nnod;
    //! Mesh connectivity extended
    std::vector< std::size_t > m_inpoel;
    //! Node coordinates extended
    tk::UnsMesh::Coords m_coord;
    //! Face data
    FaceData m_fd;
    //! Vector of unknown/solution average over each mesh element
    tk::Fields m_u;
    //! Vector of unknown at previous time-step
    tk::Fields m_un;
    //! Vector of primitive quantities over each mesh element
    tk::Fields m_p;
    //! Face geometry
    tk::Fields m_geoFace;
    //! Element geometry
    tk::Fields m_geoElem;
    //! Left-hand side mass-matrix which is a diagonal matrix
    tk::Fields m_lhs;
    //! Vector of right-hand side
    tk::Fields m_rhs;
    //! Counter for number of faces on this chare (including chare boundaries)
    std::size_t m_nfac;
    //! Counter for number of unknowns on this chare (including ghosts)
    std::size_t m_nunk;
    //! Counter for number of nodes on this chare excluding ghosts
    std::size_t m_npoin;
    //! Internal + physical boundary faces (inverse of inpofa)
    tk::UnsMesh::FaceSet m_ipface;
    //! Face & tet IDs associated to global node IDs of the face for each chare
    //! \details This map stores not only the unique faces associated to
    //!   fellow chares, but also a newly assigned local face ID and adjacent
    //!   local tet ID.
    std::unordered_map< int, FaceMap > m_bndFace;
    //! Ghost data associated to chare IDs we communicate with
    std::unordered_map< int, GhostData > m_ghostData;
    //! Elements which are ghosts for other chares associated to those chare IDs
    std::unordered_map< int, std::unordered_set< std::size_t > > m_sendGhost;
    //! Number of chares requesting ghost data
    std::size_t m_ghostReq;
    //! Local element id associated to ghost remote id charewise
    //! \details This map associates the local element id (inner map value) to
    //!    the (remote) element id of the ghost (inner map key) based on the
    //!    chare id (outer map key) this remote element lies in.
    std::unordered_map< int,
      std::unordered_map< std::size_t, std::size_t > > m_ghost;
    //! Expected ghost tet ids (used only in DEBUG)
    std::set< std::size_t > m_exptGhost;
    //! Received ghost tet ids (used only in DEBUG)
    std::set< std::size_t > m_recvGhost;
    //! Diagnostics object
    ElemDiagnostics m_diag;
    //! Runge-Kutta stage counter
    std::size_t m_stage;
    //! Map local ghost tet ids (value) and zero-based boundary ids (key)
    std::unordered_map< std::size_t, std::size_t > m_bid;
    //! Solution receive buffers for ghosts only
    std::array< std::vector< std::vector< tk::real > >, 3 > m_uc;
    //! Primitive-variable receive buffers for ghosts only
    std::array< std::vector< std::vector< tk::real > >, 3 > m_pc;
    //! 1 if starting time stepping, 0 if during time stepping
    std::size_t m_initial;
    //! Unique set of chare-boundary faces this chare is expected to receive
    tk::UnsMesh::FaceSet m_expChBndFace;
    //! Incoming communication buffer during chare-boundary face communication
    std::unordered_map< int, tk::UnsMesh::FaceSet > m_infaces;
    //! Elements (value) surrounding point (key) data-structure
    std::map< std::size_t, std::vector< std::size_t > > m_esup;
    //! Communication buffer for esup data-structure
    std::map< std::size_t, std::vector< std::size_t > > m_esupc;
    //! Elem output fields
    std::vector< std::vector< tk::real > > m_elemfields;
    //! Node output fields
    std::vector< std::vector< tk::real > > m_nodefields;
    //! Receive buffer for communication of node output fields
    //! \details Key: global node id, value: output fields and number of
    //!   elements surrounding the node
    std::unordered_map< std::size_t, std::pair< std::vector< tk::real >,
                                                std::size_t > > m_nodefieldsc;
    //! Storage for refined mesh used for field output
    OutMesh m_outmesh;
    //! Element ids at which box ICs are defined by user (multiple boxes)
    std::vector< std::unordered_set< std::size_t > > m_boxelems;

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }

    //! Compute partial boundary surface integral and sum across all chares
    void bndIntegral();

    //! Compute chare-boundary faces
    void bndFaces();

    //! Perform leak test on chare-boundary faces
    bool leakyAdjacency();

    //! Check if esuf of chare-boundary faces matches
    bool faceMatch();

    //! Verify that all chare-boundary faces have been received
    bool receivedChBndFaces();

    //! Check if entries in inpoel, inpofa and node-triplet are consistent
    std::size_t
    nodetripletMatch( const std::array< std::size_t, 2 >& id,
                      const tk::UnsMesh::Face& t );

    //! Find any chare for face (given by 3 global node IDs)
    int findchare( const tk::UnsMesh::Face& t );

    //! Setup own ghost data on this chare
    void setupGhost();

    //! Continue after face adjacency communication map completed on this chare
    void faceAdj();

    //! Continue after node adjacency communication map completed on this chare
    void adj();

    //! Fill elements surrounding a face along chare boundary
    void addEsuf( const std::array< std::size_t, 2 >& id, std::size_t ghostid );

    //! Fill elements surrounding a element along chare boundary
    void addEsuel( const std::array< std::size_t, 2 >& id,
                   std::size_t ghostid,
                   const tk::UnsMesh::Face& t );

    void addEsup();

    //! Fill face geometry data along chare boundary
    void addGeoFace( const tk::UnsMesh::Face& t,
                     const std::array< std::size_t, 2 >& id );

    //! Output mesh field data
    void writeFields( CkCallback c );

    //! Compute solution reconstructions
    void reco();

    //! Compute limiter function
    void lim();

    //! Compute time step size
    void dt();

    //! Evaluate whether to continue with next time step stage
    void stage();

    //! Evaluate whether to save checkpoint/restart
    void evalRestart();

    //! Evaluate solution on incomping (a potentially refined) mesh
    std::tuple< tk::Fields, tk::Fields, tk::Fields, tk::Fields >
    evalSolution(
      const std::vector< std::size_t >& inpoel,
      const tk::UnsMesh::Coords& coord,
      const std::unordered_map< std::size_t, std::size_t >& addedTets );

    //! Decide wether to output field data
    bool fieldOutput() const;

    //! Decide if we write field output using a refined mesh
    bool refinedOutput() const;

    //! Start preparing fields for output to file
    void startFieldOutput( CkCallback c );
};

} // inciter::

#endif // FV_h
