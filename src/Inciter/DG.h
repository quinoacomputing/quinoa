// *****************************************************************************
/*!
  \file      src/Inciter/DG.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
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

    #### Call graph ####
    The following is a directed acyclic graph (DAG) that outlines the
    asynchronous algorithm implemented in this class The detailed discussion of
    the algorithm is given in the Charm++ interface file transporter.ci. On the
    DAG orange fills denote global synchronization points that contain or
    eventually lead to global reductions. Dashed lines are potential shortcuts
    that allow jumping over some of the task-graph under some circumstances or
    optional code paths (taken, e.g., only in DEBUG mode). See the detailed
    discussion in dg.ci.
    \dot
    digraph "DG SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      OwnGhost [ label="OwnGhost"
               tooltip="own ghost data computed"
               URL="\ref inciter::DG::setupGhost"];
      ReqGhost [ label="ReqGhost"
               tooltip="all of ghost data have been requested from us"
               URL="\ref inciter::DG::reqGhost"];
      OwnGhost -> sendGhost [ style="solid" ];
      ReqGhost -> sendGhost [ style="solid" ];
      OwnSol [ label="OwnSol"
               tooltip="own solution/unknown data computed"
               URL="\ref inciter::DG::advance"];
      ComSol [ label="ComSol"
               tooltip="communicated (ghost) solution/unknown data received"
               URL="\ref inciter::DG::comsol"];
      OwnSol -> Solve [ style="solid" ];
      ComSol -> Solve [ style="solid" ];
    }
    \enddot
    \include Inciter/dg.ci

*/
// *****************************************************************************
#ifndef DG_h
#define DG_h

#include <array>
#include <set>
#include <unordered_set>
#include <unordered_map>

#include "DerivedData.h"
#include "FaceData.h"
#include "ElemDiagnostics.h"

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
                 const tk::CProxy_Solver&,
                 const FaceData& fd );

    //! Migrate constructor
    explicit DG( CkMigrateMessage* ) {}

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
    void setup( tk::real v );

    //! Compute time step size
    void dt();

    //! Receive chare-boundary ghost data from neighboring chares
    void comsol( int fromch,
                 const std::vector< std::size_t >& tetid,
                 const std::vector< std::vector< tk::real > >& u );

    //! Evaluate whether to continue with next step
    void eval();

    //! Advance equations to next time step
    void advance( tk::real newdt );

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_DG::pup(p);
      p | m_ncomfac;
      p | m_nadj;
      p | m_nsol;
      p | m_itf;
      p | m_disc;
      p | m_fd;
      p | m_u;
      p | m_un;
      p | m_vol;
      p | m_geoFace;
      p | m_geoElem;
      p | m_lhs;
      p | m_rhs;
      p | m_nfac;
      p | m_nunk;
      p | m_msumset;
      p | m_esuelTet;
      p | m_ipface;
      p | m_bndFace;
      p | m_ghostData;
      p | m_ghostReq;
      p | m_exptNbface;
      p | m_ghost;
      p | m_exptGhost;
      p | m_recvGhost;
      p | m_diag;
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
                          tk::UnsMesh::FaceHasher,
                          tk::UnsMesh::FaceEq >;

    //! Counter for face adjacency communication map
    std::size_t m_ncomfac;
    //! Counter signaling that all ghost data have been received
    std::size_t m_nadj;
    //! Counter signaling that we have received all our solution ghost data
    std::size_t m_nsol;
    //! Field output iteration count
    uint64_t m_itf;
    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Face data
    FaceData m_fd;
    //! Vector of unknown/solution average over each mesh element
    tk::Fields m_u;
    //! Vector of unknown at previous time-step
    tk::Fields m_un;
    //! Total mesh volume
    tk::real m_vol;
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
    //! \brief Global mesh node IDs bordering the mesh chunk held by fellow
    //!    worker chares associated to their chare IDs
    //! \details msum: mesh chunks surrounding mesh chunks and their neighbor
    //!   points. This is the same data as in Discretization::m_msum, but the
    //!   nodelist is stored as a set.
    std::unordered_map< int, std::unordered_set< std::size_t > > m_msumset;
    //! Elements surrounding elements with -1 at boundaries, see genEsuelTet()
    std::vector< int > m_esuelTet;
    //! Internal + physical boundary faces (inverse of inpofa)
    tk::UnsMesh::FaceSet m_ipface;
    //! Face * tet IDs associated to global node IDs of the face for each chare
    //! \details This maps stores not only the unique faces associated to
    //!   fellow chares, but also a newly assigned local face ID and adjacent
    //!   local tet ID.
    std::unordered_map< int, FaceMap > m_bndFace;
    //! Ghost data associated to chare IDs we communicate with
    std::unordered_map< int, GhostData > m_ghostData;
    //! Number of chares requesting ghost data
    std::size_t m_ghostReq;
    //! Expected number of boundary faces (used only in DEBUG)
    std::size_t m_exptNbface;
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

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }

    //! Perform leak test on mesh partition
    bool leakyPartition();

    //! Perform leak test on chare-boundary faces
    bool leakyAdjacency();

    //! Find any chare for face (given by 3 global node IDs)
    int findchare( const tk::UnsMesh::Face& t );

    //! Setup own ghost data on this chare
    void setupGhost();

    //! Convert chare-node adjacency map to hold sets instead of vectors
    std::unordered_map< int, std::unordered_set< std::size_t > >
    msumset() const;

    //! Continue after face adjacency communication map completed on this chare
    void adj();

    //! Fill elements surrounding a face along chare boundary
    void addEsuf( const std::array< std::size_t, 2 >& id, std::size_t ghostid );

    //! Fill face geometry data along chare boundary
    void addGeoFace( const tk::UnsMesh::Face& t,
                     const std::array< std::size_t, 2 >& id );

    //! Compute left hand side
    void lhs();

    //! Compute right hand side and solve system
    void solve();

    //! Output mesh and particle fields to files
    void out();

    //! Output mesh-based fields to file
    void writeFields( tk::real time );
};

} // inciter::

#endif // DG_h
