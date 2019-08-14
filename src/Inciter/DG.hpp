// *****************************************************************************
/*!
  \file      src/Inciter/DG.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     DG advances a system of PDEs with the discontinuous Galerkin scheme
  \details   DG advances a system of partial differential equations (PDEs) using
    discontinuous Galerkin (DG) finite element (FE) spatial discretization (on
    tetrahedron elements) combined with Runge-Kutta (RK) time stepping.

    There are a potentially large number of DG Charm++ chares created by
    Transporter. Each DG gets a chunk of the full load (part of the mesh)
    and does the same: initializes and advances a number of PDE systems in time.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/dg.ci.
*/
// *****************************************************************************
#ifndef DG_h
#define DG_h

#include <array>
#include <set>
#include <unordered_set>
#include <unordered_map>

#include "DerivedData.hpp"
#include "FaceData.hpp"
#include "ElemDiagnostics.hpp"

#include "NoWarning/dg.decl.h"

namespace inciter {

//! DG Charm++ chare array used to advance PDEs in time with DG+RK
class DG : public CBase_DG {

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
    DG_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit DG( const CProxy_Discretization& disc,
                 const std::map< int, std::vector< std::size_t > >& bface,
                 const std::map< int, std::vector< std::size_t > >& /* bnode */,
                 const std::vector< std::size_t >& triinpoel );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit DG( CkMigrateMessage* msg ) : CBase_DG( msg ) {}
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

    //! Configure Charm++ reduction types for concatenating BC nodelists
    static void registerReducers();

    //! Setup: query boundary conditions, output mesh, etc.
    void setup();

    // Evaluate whether to do load balancing
    void evalLB();

    //! Continue to next time step
    void next();

    //! Receive chare-boundary limiter function data from neighboring chares
    void comlim( int fromch,
                 const std::vector< std::size_t >& tetid,
                 const std::vector< std::vector< tk::real > >& u,
                 const std::vector< std::vector< tk::real > >& prim,
                 const std::vector< std::size_t >& ndof );

    //! Receive chare-boundary ghost data from neighboring chares
    void comsol( int fromch,
                 std::size_t fromstage,
                 const std::vector< std::size_t >& tetid,
                 const std::vector< std::vector< tk::real > >& u,
                 const std::vector< std::vector< tk::real > >& prim,
                 const std::vector< std::size_t >& ndof );

    //! Optionally refine/derefine mesh
    void refine();

    //! Receive new mesh from refiner
    void resizePostAMR(
      const std::vector< std::size_t >& /* ginpoel */,
      const tk::UnsMesh::Chunk& chunk,
      const tk::UnsMesh::Coords& coord,
      const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& /* addedNodes */,
      const std::unordered_map< std::size_t, std::size_t >& addedTets,
      const std::unordered_map< int, std::vector< std::size_t > >& msum,
      const std::map< int, std::vector< std::size_t > >& bface,
      const std::map< int, std::vector< std::size_t > >& /* bnode */,
      const std::vector< std::size_t >& triinpoel );

    //! Const-ref access to current solution
    //! \return Const-ref to current solution
    const tk::Fields& solution() const { return m_u; }

    //! Compute left hand side
    void lhs();

    //! Compute limiter function
    void lim();

    //! Const-ref access to current solution
    //! \param[in,out] u Reference to update with current solution
    void solution( tk::Fields& u ) const { u = m_u; }

    //! Unused in DG
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
      p | m_nsol;
      p | m_ninitsol;
      p | m_nlim;
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
      p | m_ncoord;
      p | m_msumset;
      p | m_ipface;
      p | m_bndFace;
      p | m_ghostData;
      p | m_ghostReq;
      p | m_ghost;
      p | m_exptGhost;
      p | m_recvGhost;
      p | m_diag;
      p | m_stage;
      p | m_ndof;
      p | m_bid;
      p | m_uc;
      p | m_pc;
      p | m_ndofc;
      p | m_initial;
      p | m_expChBndFace;
      p | m_infaces;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i DG object reference
    friend void operator|( PUP::er& p, DG& i ) { i.pup(p); }
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

    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Counter for face adjacency communication map
    std::size_t m_ncomfac;
    //! Counter signaling that all ghost data have been received
    std::size_t m_nadj;
    //! Counter signaling that we have received all our solution ghost data
    std::size_t m_nsol;
    //! \brief Counter signaling that we have received all our solution ghost
    //!    data during setup
    std::size_t m_ninitsol;
    //! Counter signaling that we have received all our limiter function ghost data
    std::size_t m_nlim;
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
    std::size_t m_ncoord;
    //! \brief Global mesh node IDs bordering the mesh chunk held by fellow
    //!    worker chares associated to their chare IDs
    //! \details msum: mesh chunks surrounding mesh chunks and their neighbor
    //!   points. This is the same data as in Discretization::m_msum, but the
    //!   nodelist is stored as a set.
    std::unordered_map< int, std::unordered_set< std::size_t > > m_msumset;
    //! Internal + physical boundary faces (inverse of inpofa)
    tk::UnsMesh::FaceSet m_ipface;
    //! Face & tet IDs associated to global node IDs of the face for each chare
    //! \details This map stores not only the unique faces associated to
    //!   fellow chares, but also a newly assigned local face ID and adjacent
    //!   local tet ID.
    std::unordered_map< int, FaceMap > m_bndFace;
    //! Ghost data associated to chare IDs we communicate with
    std::unordered_map< int, GhostData > m_ghostData;
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
    //! Vector of local number of degrees of freedom for each element
    std::vector< std::size_t > m_ndof;
    //! Map local ghost tet ids (value) and zero-based boundary ids (key)
    std::unordered_map< std::size_t, std::size_t > m_bid;
    //! Solution receive buffers for ghosts only
    std::array< std::vector< std::vector< tk::real > >, 2 > m_uc;
    //! Primitive-variable receive buffers for ghosts only
    std::array< std::vector< std::vector< tk::real > >, 2 > m_pc;
    //! \brief Number of degrees of freedom (for p-adaptive) receive buffers
    //!   for ghosts only
    std::array< std::vector< std::size_t >, 2 > m_ndofc;
    //! 1 if starting time stepping, 0 if during time stepping
    int m_initial;
    //! Unique set of chare-boundary faces this chare is expected to receive
    tk::UnsMesh::FaceSet m_expChBndFace;
    //! Incoming communication buffer during chare-boundary face communication
    std::unordered_map< int, tk::UnsMesh::FaceSet > m_infaces;

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
    void adj();

    //! Fill elements surrounding a face along chare boundary
    void addEsuf( const std::array< std::size_t, 2 >& id, std::size_t ghostid );

    //! Fill elements surrounding a element along chare boundary
    void addEsuel( const std::array< std::size_t, 2 >& id,
                   std::size_t ghostid,
                   const tk::UnsMesh::Face& t );

    //! Fill face geometry data along chare boundary
    void addGeoFace( const tk::UnsMesh::Face& t,
                     const std::array< std::size_t, 2 >& id );

    //! Output mesh and particle fields to files
    void out();

    //! Output mesh-based fields to file
    void writeFields( CkCallback c ) const;

    //! Compute time step size
    void dt();

    //! Evaluate whether to continue with next time step stage
    void stage();

    //! Evaluate whether to save checkpoint/restart
    void evalRestart();

    //! p-refine all elements that are adjacent to p-refined elements
    void propagate_ndof();
};

} // inciter::

#endif // DG_h
