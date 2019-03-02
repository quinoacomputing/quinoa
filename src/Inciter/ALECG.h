// *****************************************************************************
/*!
  \file      src/Inciter/ALECG.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     ALECG for a PDE system with continuous Galerkin + ALE + RK
  \details   ALECG advances a system of partial differential equations (PDEs)
    using a continuous Galerkin (CG) finite element (FE) spatial discretization
    (using linear shapefunctions on tetrahedron elements) combined with a
    Runge-Kutta (RK) time stepping scheme in the arbitrary Eulerian-Lagrangian
    reference frame.

    There are a potentially large number of ALECG Charm++ chares created by
    Transporter. Each ALECG gets a chunk of the full load (part of the mesh)
    and does the same: initializes and advances a number of PDE systems in time.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality.
*/
// *****************************************************************************
#ifndef ALECG_h
#define ALECG_h

#include <vector>
#include <map>

#include "QuinoaConfig.h"
#include "Types.h"
#include "Fields.h"
#include "DerivedData.h"
#include "FluxCorrector.h"
#include "NodeDiagnostics.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "FaceData.h"

#include "NoWarning/alecg.decl.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! ALECG Charm++ chare array used to advance PDEs in time with ALECG+RK
class ALECG : public CBase_ALECG {

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
    ALECG_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit ALECG( const CProxy_Discretization& disc, const FaceData& fd );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit ALECG( CkMigrateMessage* ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ custom reduction types initiated from this chare array
    static void registerReducers();

    //! Setup: query boundary conditions, output mesh, etc.
    void setup( tk::real v );

    // Initially compute left hand side diagonal matrix
    void init();

    //! Send own chare-boundary data to neighboring chares
    void sendinit(){}

    //! Advance equations to next time step
    void advance( tk::real newdt );

    //! Compute left-hand side of transport equations
    void lhs();

    //! Receive contributions to left-hand side matrix on chare-boundaries
    void comlhs( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& L );

    //! Receive contributions to right-hand side vector on chare-boundaries
    void comrhs( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& R );

    //! Update solution at the end of time step
    void update( const tk::Fields& a );

    //! Signal the runtime system that diagnostics have been computed
    void diag();

    //! Optionally refine/derefine mesh
    void refine();

    //! Receive new mesh from refiner
    void resize( const tk::UnsMesh::Chunk& chunk,
                 const tk::UnsMesh::Coords& coord,
                 const tk::Fields& u,
                 const std::unordered_map< int,
                         std::vector< std::size_t > >& msum,
                 const std::map< int, std::vector< std::size_t > >& bnode );

    //! Const-ref access to current solution
    //! \param[in,out] u Reference to update with current solution
    void solution( tk::Fields& u ) const { u = m_u; }

    //! Resizing data sutrctures after mesh refinement has been completed
    void resized();

    //! Evaluate whether to continue with next time step
    void step();

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_disc;
      p | m_initial;
      p | m_nsol;
      p | m_nlhs;
      p | m_nrhs;
      p | m_fd;
      p | m_u;
      p | m_du;
      p | m_lhs;
      p | m_rhs;
      p | m_lhsc;
      p | m_rhsc;
      p | m_vol;
      p | m_diag;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i ALECG object reference
    friend void operator|( PUP::er& p, ALECG& i ) { i.pup(p); }
    //@}

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! True if starting time stepping, false if during time stepping
    bool m_initial;
    //! Counter for high order solution vector nodes updated
    std::size_t m_nsol;
    //! Counter for left-hand side matrix (vector) nodes updated
    std::size_t m_nlhs;
    //! Counter for right-hand side vector nodes updated
    std::size_t m_nrhs;
    //! Face data
    FaceData m_fd;
    //! Unknown/solution vector at mesh nodes
    tk::Fields m_u;
    //! Unknown/solution vector increment (high order)
    tk::Fields m_du;
    //! Lumped lhs mass matrix
    tk::Fields m_lhs;
    //! Right-hand side vector (for the high order system)
    tk::Fields m_rhs;
    //! Receive buffers for communication
    std::vector< std::vector< tk::real > > m_lhsc, m_rhsc;
    //! Total mesh volume
    tk::real m_vol;
    //! Diagnostics object
    NodeDiagnostics m_diag;

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }

    //! Size communication buffers
    void resizeComm();

    //! Output mesh and particle fields to files
    void out();

    //! Output mesh-based fields to file
    void writeFields( CkCallback c );

    //! The own and communication portion of the left-hand side is complete
    void lhsdone();

    //! Combine own and communicated contributions to left hand side
    void lhsmerge();

    //! Compute righ-hand side vector of transport equations
    void rhs();

    //! Start time stepping
    void start();

    //! Solve low and high order diagonal systems
    void solve();

    //! Compute time step size
    void dt();
};

} // inciter::

#endif // ALECG_h
