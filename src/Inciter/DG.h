// *****************************************************************************
/*!
  \file      src/Inciter/DG.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     DG advances a system of PDEs with the discontinuous Galerkin scheme
  \details   DG advances a system of partial differential equations (PDEs) using
    discontinuous Galerkin (DG) finite element (FE) spatial discretization (on
    tetrahedron elements) combined with Runge-Kutta (RK) time stepping.

    There are a potentially large number of DG Charm++ chares created by
    Transporter. Each DG gets a chunk of the full load (part of the mesh)
    and does the same: initializes and advances a number of PDE systems in time.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication.
*/
// *****************************************************************************
#ifndef DG_h
#define DG_h

#include <array>
#include <unordered_set>
#include <unordered_map>

#include "DerivedData.h"
#include "FaceData.h"

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
                 const tk::CProxy_Solver& solver,
                 const FaceData& fd );

    //! Migrate constructor
    explicit DG( CkMigrateMessage* ) {}

    //! Receive ghost data on chare boundaries from fellow chare
    void comadj( int fromch, const GhostData& ghost );

    //! Configure Charm++ reduction types for concatenating BC nodelists
    static void registerReducers();

    //! Setup: query boundary conditions, output mesh, etc.
    void setup( tk::real v );

    //! Compute time step size
    void dt();

    //! Advance equations to next time step
    void advance( tk::real newdt );

    //! Evaluate whether to continue with next step
    void eval();

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_DG::pup(p);
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
      p | m_msumset;
      p | m_ghost;
      p | m_chBndFace;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i DG object reference
    friend void operator|( PUP::er& p, DG& i ) { i.pup(p); }
    //@}

  private:
    //! Node ID triplet denoting a tetrahedron face
    using Triplet = std::array< std::size_t, 3 >;
    // Hash functor for node Triplet
    struct TripletHasher {
      std::size_t operator()( const Triplet& key ) const {
        return std::hash< std::size_t >()( key[0] ) ^
               std::hash< std::size_t >()( key[1] ) ^
               std::hash< std::size_t >()( key[2] );
      }
    };
   //! \brief Key-equal function for node triplet in which the order matters but
   //!   not the positions
    struct TripletEq {
      bool operator()( const Triplet& l, const Triplet& r ) const {
        return (l[0] == r[0] && l[1] == r[1] && l[2] == r[2]) ||
               (l[0] == r[1] && l[1] == r[2] && l[2] == r[0]) ||
               (l[0] == r[2] && l[1] == r[0] && l[2] == r[1]);
      }
    };
    //! Node ID triplets representing a tetrahedron face
    using Faces = std::unordered_set< Triplet, TripletHasher, TripletEq >;
    //! Tetrahedron face ID associated to node ID triplet
    using FaceIDs =
      std::unordered_map< Triplet, std::size_t, TripletHasher, TripletEq >;

    tk::CProxy_Solver m_solver;
    //! Counter for face adjacency communication map
    std::size_t m_nadj;
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
    //! \brief Global mesh node IDs bordering the mesh chunk held by fellow
    //!    worker chares associated to their chare IDs
    //! \details msum: mesh chunks surrounding mesh chunks and their neighbor
    //!   points. This is the same data as in Discretization::m_msum, but the
    //!   nodelist is stored as a set.
    std::unordered_map< int, std::unordered_set< std::size_t > > m_msumset;
    //! Local element id associated to ghost remote id
    //! \details This map associates the local element id (map value) to the
    //!    (remote) element id of the ghost (map key).
    std::unordered_map< std::size_t, std::size_t > m_ghost;
    //! ...
    FaceIDs m_chBndFace;

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }

    //! ...
    std::unordered_map< int, std::unordered_set< std::size_t > >
    msumset( const std::vector< std::size_t >& inpofa ) const;

    //! Continue after face adjacency communication map is complete on this chare
    void adj();

    //! Compute left hand side
    void lhs();

    //! Compute right hand side
    void rhs();

    //! Time stepping
    void solve( tk::real deltat );

    //! Prepare for next step
    void next();
 
    //! Output mesh and particle fields to files
    void out();

    //! Compute diagnostics, e.g., residuals
    bool diagnostics();

    //! Output mesh-based fields to file
    void writeFields( tk::real time );
};

} // inciter::

#endif // DG_h
