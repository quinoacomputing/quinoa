// *****************************************************************************
/*!
  \file      src/Inciter/Ghosts.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Declarations file for generating ghost data structures
  \details   Declarations file for asynchronous distributed
             ghost data structures using Charm++.

    There are a potentially large number of Ghosts Charm++ chares.
    Each Ghosts chare gets a chunk of the full load, due to partiting the mesh.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality.
*/
// *****************************************************************************
#ifndef Ghosts_h
#define Ghosts_h

#include "Fields.hpp"
#include "FaceData.hpp"
#include "Discretization.hpp"

#include "NoWarning/ghosts.decl.h"

namespace inciter {

//! Ghosts Charm++ chare array used to determine ghost data structures
class Ghosts : public CBase_Ghosts {

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
    Ghosts_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit
    Ghosts( const CProxy_Discretization& disc,
      const std::map< int, std::vector< std::size_t > >& bface,
      const std::vector< std::size_t >& triinpoel,
      std::size_t nunk );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit Ghosts( CkMigrateMessage* ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

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
    //! Counter for number of unknowns on this chare (including ghosts)
    std::size_t m_nunk;
    //! Mesh connectivity extended
    std::vector< std::size_t > m_inpoel;
    //! Node coordinates extended
    tk::UnsMesh::Coords m_coord;
    //! Face data
    FaceData m_fd;
    //! Face geometry
    tk::Fields m_geoFace;
    //! Element geometry
    tk::Fields m_geoElem;
    //! Counter for number of faces on this chare (including chare boundaries)
    std::size_t m_nfac;
    //! Face & tet IDs associated to global node IDs of the face for each chare
    //! \details This map stores not only the unique faces associated to
    //!   fellow chares, but also a newly assigned local face ID and adjacent
    //!   local tet ID.
    std::unordered_map< int, FaceMap > m_bndFace;
    //! Elements which are ghosts for other chares associated to those chare IDs
    std::unordered_map< int, std::unordered_set< std::size_t > > m_sendGhost;
    //! Local element id associated to ghost remote id charewise
    //! \details This map associates the local element id (inner map value) to
    //!    the (remote) element id of the ghost (inner map key) based on the
    //!    chare id (outer map key) this remote element lies in.
    std::unordered_map< int,
      std::unordered_map< std::size_t, std::size_t > > m_ghost;
    //! Expected ghost tet ids (used only in DEBUG)
    std::set< std::size_t > m_exptGhost;
    //! Map local ghost tet ids (value) and zero-based boundary ids (key)
    std::unordered_map< std::size_t, std::size_t > m_bid;
    //! Elements (value) surrounding point (key) data-structure
    std::map< std::size_t, std::vector< std::size_t > > m_esup;

    //1 Start setup of communication maps for cell-centered schemes
    void startCommSetup();

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

    /** @name Pack/unpack (Charm++ serialization) routines */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_disc;
      p | m_nunk;
      p | m_inpoel;
      p | m_coord;
      p | m_fd;
      p | m_geoFace;
      p | m_geoElem;
      p | m_nfac;
      p | m_bndFace;
      p | m_sendGhost;
      p | m_ghost;
      p | m_exptGhost;
      p | m_bid;
      p | m_esup;
      p | m_ncomfac;
      p | m_nadj;
      p | m_ncomEsup;
      p | m_ipface;
      p | m_ghostData;
      p | m_ghostReq;
      p | m_initial;
      p | m_expChBndFace;
      p | m_infaces;
      p | m_esupc;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] a Ghosts object reference
    friend void operator|( PUP::er& p, Ghosts& a ) { a.pup(p); }
    ///@}

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

    //! Counter for face adjacency communication map
    std::size_t m_ncomfac;
    //! Counter signaling that all ghost data have been received
    std::size_t m_nadj;
    //! Counter for element-surr-node adjacency communication map
    std::size_t m_ncomEsup;
    //! Internal + physical boundary faces (inverse of inpofa)
    tk::UnsMesh::FaceSet m_ipface;
    //! Ghost data associated to chare IDs we communicate with
    std::unordered_map< int, GhostData > m_ghostData;
    //! Number of chares requesting ghost data
    std::size_t m_ghostReq;
    //! 1 if starting time stepping, 0 if during time stepping
    std::size_t m_initial;
    //! Unique set of chare-boundary faces this chare is expected to receive
    tk::UnsMesh::FaceSet m_expChBndFace;
    //! Incoming communication buffer during chare-boundary face communication
    std::unordered_map< int, tk::UnsMesh::FaceSet > m_infaces;
    //! Communication buffer for esup data-structure
    std::map< std::size_t, std::vector< std::size_t > > m_esupc;

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }

    //! Compute chare-boundary faces
    void bndFaces();

    //! Setup own ghost data on this chare
    void setupGhost();

    //! Continue after face adjacency communication map completed on this chare
    void faceAdj();

    //! Compute partial boundary surface integral and sum across all chares
    void bndIntegral();

    //! Continue after node adjacency communication map completed on this chare
    void adj();

    //! Perform leak-test on chare boundary faces
    bool leakyAdjacency();

    //! Check if esuf of chare-boundary faces matches
    bool faceMatch();

    //! Verify that all chare-boundary faces have been received
    bool receivedChBndFaces();

    //! Find any chare for face (given by 3 global node IDs)
    int findchare( const tk::UnsMesh::Face& t );

    //! Check if entries in inpoel, inpofa and node-triplet are consistent
    std::size_t nodetripletMatch(
      const std::array< std::size_t, 2 >& id,
      const tk::UnsMesh::Face& t );

    //! Fill elements surrounding a face along chare boundary
    void addEsuf(
      const std::array< std::size_t, 2 >& id,
      std::size_t ghostid );

    //! Fill elements surrounding a element along chare boundary
    void addEsuel(
      const std::array< std::size_t, 2 >& id,
      std::size_t ghostid,
      const tk::UnsMesh::Face& t );

    //! Fill face-geometry data along chare boundary
    void addGeoFace(
      const tk::UnsMesh::Face& t,
      const std::array< std::size_t, 2 >& id );
};

} // inciter::

#endif // Ghosts_h
