// *****************************************************************************
/*!
  \file      src/Inciter/DiagCG.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     DiagCG for a PDE system with continuous Galerkin without a matrix
  \details   DiagCG advances a system of partial differential equations (PDEs)
    using continuous Galerkin (CG) finite element (FE) spatial discretization
    (using linear shapefunctions on tetrahedron elements) combined with a time
    stepping scheme that is equivalent to the Lax-Wendroff (LW) scheme within
    the unstructured-mesh FE context and treats discontinuities with
    flux-corrected transport (FCT). The left-hand side (lumped-mass) matrix
    is diagonal thus this scheme does not use a matrix-based linear solver.

    There are a potentially large number of DiagCG Charm++ chares created by
    Transporter. Each DiagCG gets a chunk of the full load (part of the mesh)
    and does the same: initializes and advances a number of PDE systems in time.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/diagcg.ci.

    #### Call graph ####
    The following is a directed acyclic graph (DAG) that outlines the
    asynchronous algorithm implemented in this class The detailed discussion of
    the algorithm is given in the Charm++ interface file transporter.ci, which
    also repeats the graph below using ASCII graphics. On the DAG orange
    fills denote global synchronization points that contain or eventually lead
    to global reductions. Dashed lines are potential shortcuts that allow
    jumping over some of the task-graph under some circumstances or optional
    code paths (taken, e.g., only in DEBUG mode). See the detailed discussion in
    diagcg.ci.
    \dot
    digraph "DiagCG SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      OwnLhs [ label="OwnLhs"
               tooltip="own contributions to the left hand side lumped-mass
                        matrix computed"
               URL="\ref inciter::DiagCG::lhs"];
      ComLhs [ label="ComLhs"
               tooltip="contributions to the left hand side lumped-mass matrix
                        communicated"
               URL="\ref inciter::DiagCG::comlhs"];
      OwnRhs [ label="OwnRhs"
               tooltip="own contributions to the right hand side computed"
               URL="\ref inciter::DiagCG::rhs"];
      ComRhs [ label="ComRhs"
               tooltip="contributions to the right hand side communicated"
               URL="\ref inciter::DiagCG::comrhs"];
      OwnDif [ label="OwnDif"
               tooltip="own contributions to the mass diffusion rhs computed"
               URL="\ref inciter::DiagCG::rhs"];
      ComDif [ label="ComDif"
               tooltip="contributions to the mass diffusion rhs communicated"
               URL="\ref inciter::DiagCG::comdif"];
      Start [ label="Ver" tooltip="start time stepping"
              URL="\ref inciter::DiagCG::start"];
      Solve [ label="Ver" tooltip="solve diagonal systems"
              URL="\ref inciter::DiagCG::solve"];
      OwnLhs -> Start [ style="solid" ];
      ComLhs -> Start [ style="solid" ];
      OwnRhs -> Solve [ style="solid" ];
      ComRhs -> Solve [ style="solid" ];
      OwnDif -> Solve [ style="solid" ];
      ComDif -> Solve [ style="solid" ];
    }
    \enddot
    \include Inciter/diagcg.ci
*/
// *****************************************************************************
#ifndef DiagCG_h
#define DiagCG_h

#include <vector>
#include <map>
#include <unordered_set>

#include "QuinoaConfig.h"
#include "Types.h"
#include "Fields.h"
#include "DerivedData.h"
#include "VectorReducer.h"
#include "FluxCorrector.h"
#include "NodeDiagnostics.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "FaceData.h"

#include "NoWarning/diagcg.decl.h"

namespace tk {
  class ExodusIIMeshWriter;
  class RootMeshWriter;
}

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! DiagCG Charm++ chare array used to advance PDEs in time with DiagCG+LW+FCT
class DiagCG : public CBase_DiagCG {

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
    DiagCG_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit DiagCG( const CProxy_Discretization& disc,
                     const tk::CProxy_Solver&,
                     const FaceData& );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    explicit DiagCG( CkMigrateMessage* ) : m_diag( *Disc() ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ custom reduction types initiated from this chare array
    static void registerReducers();

    //! Setup: query boundary conditions, output mesh, etc.
    void setup( tk::real v );

    //! Compute time step size
    void dt();

    //! Advance equations to next time step
    void advance( tk::real newdt );

    //! Receive contributions to left-hand side matrix on chare-boundaries
    void comlhs( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& L );

    //! Receive contributions to right-hand side vector on chare-boundaries
    void comrhs( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& R );

    //!  Receive contributions to RHS mass diffusion on chare-boundaries
    void comdif( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& D );

    //! Verify that solution does not change at Dirichlet boundary conditions
    bool correctBC( const tk::Fields& a );

    //! Prepare for next step    
    void next( const tk::Fields& a );

    //! Evaluate whether to continue with next step
    void eval();

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_DiagCG::pup(p);
      p | m_itf;
      p | m_nsol;
      p | m_nlhs;
      p | m_nrhs;
      p | m_ndif;
      p | m_disc;
      p | m_solver;
      p | m_side;
      p | m_u;
      p | m_ul;
      p | m_du;
      p | m_dul;
      p | m_ue;
      p | m_lhs;
      p | m_rhs;
      p | m_dif;
      p | m_lhsc;
      p | m_rhsc;
      p | m_difc;
      p | m_vol;
      p | m_diag;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i DiagCG object reference
    friend void operator|( PUP::er& p, DiagCG& i ) { i.pup(p); }
    //@}

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

    //! Field output iteration count
    uint64_t m_itf;
    //! Counter for high order solution vector nodes updated
    std::size_t m_nsol;
    //! Counter for left-hand side matrix (vector) nodes updated
    std::size_t m_nlhs;
    //! Counter for right-hand side vector nodes updated
    std::size_t m_nrhs;
    //! Counter for right-hand side masss-diffusion vector nodes updated
    std::size_t m_ndif;
    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Linear system merger and solver proxy
    tk::CProxy_Solver m_solver;
    //! Map associating local node IDs to side set IDs
    std::map< int, std::vector< std::size_t > > m_side;
    //! Unknown/solution vector at mesh nodes
    tk::Fields m_u;
    //! Unknown/solution vector at mesh nodes (low orderd)
    tk::Fields m_ul;
    //! Unknown/solution vector increment (high order)
    tk::Fields m_du;
    //! Unknown/solution vector increment (low order)
    tk::Fields m_dul;
    //! Unknown/solution vector at mesh cells
    tk::Fields m_ue;
    //! Lumped lhs mass matrix
    tk::Fields m_lhs;
    //! Right-hand side vector (for the high order system)
    tk::Fields m_rhs;
    //! Mass diffusion right-hand side vector (for the low order system)
    tk::Fields m_dif;
    //! Boundary conditions evaluated and assigned to mesh node IDs
    //! \details Vector of pairs of bool and boundary condition value associated
    //!   to meshnode IDs at which the user has set Dirichlet boundary
    //!   conditions for all PDEs integrated. The bool indicates whether the BC
    //!   is set at the node for that component the if true, the real value is
    //!   the increment (from t to dt) in the BC specified for a component.
    std::unordered_map< std::size_t,
      std::vector< std::pair< bool, tk::real > > > m_bc;
    //! Receive buffers for communication
    std::vector< std::vector< tk::real > > m_lhsc, m_rhsc, m_difc;
    //! Total mesh volume
    tk::real m_vol;
    //! Diagnostics object
    NodeDiagnostics m_diag;

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }

    //! Output mesh and particle fields to files
    void out();

    //! Output mesh-based fields to file
    void writeFields( tk::real time );

    //! \brief Extract node IDs from side set node lists and match to
    //    user-specified boundary conditions
    void bc();

    //! Compute left-hand side of transport equations
    void lhs();

    //! Compute righ-hand side vector of transport equations
    void rhs();

    //! Start time stepping
    void start();

    //! Solve low and high order diagonal systems
    void solve();
};

} // inciter::

#endif // DiagCG_h
