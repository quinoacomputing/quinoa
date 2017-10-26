// *****************************************************************************
/*!
  \file      src/Inciter/CG.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     CG advances a system of PDEs with the continuous Galerkin scheme
  \details   CG advances a system of partial differential equations (PDEs) using
    continuous Galerkin (CG) finite element (FE) spatial discretization (using
    linear shapefunctions on tetrahedron elements) combined with a time stepping
    scheme that is equivalent to the Lax-Wendroff (LW) scheme within the
    unstructured-mesh FE context and treats discontinuities with flux-corrected
    transport (FCT).

    There are a potentially large number of CG Charm++ chares created by
    Transporter. Each CG gets a chunk of the full load (part of the mesh)
    and does the same: initializes and advances a number of PDE systems in time.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/cg.ci.

    #### Call graph ####
    The following is a directed acyclic graph (DAG) that outlines the
    asynchronous algorithm implemented in this class The detailed discussion of
    the algorithm is given in the Charm++ interface file transporter.ci, which
    also repeats the graph below using ASCII graphics. On the DAG orange
    fills denote global synchronization points that contain or eventually lead
    to global reductions. Dashed lines are potential shortcuts that allow
    jumping over some of the task-graph under some circumstances or optional
    code paths (taken, e.g., only in DEBUG mode). See the detailed discussion in
    cg.ci.
    \dot
    digraph "CG SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      OwnAEC [ label="OwnAEC"
               tooltip="own contributions to the antidiffusive element
                        contributions computed"
               URL="\ref inciter::CG::aec"];
      ComAEC [ label="ComAEC"
               tooltip="contributions to the antidiffusive element contributions
                        communicated"
               URL="\ref inciter::CG::comaec"];
      OwnALW [ label="OwnALW"
               tooltip="own contributions to the maximum and minimum unknowns of
                        elements surrounding nodes computed"
               URL="\ref inciter::CG::alw"];
      ComALW [ label="ComALW"
               tooltip="contributions to the the maximum and minimum unknowns of
                        elements surrounding nodes communicated"
               URL="\ref inciter::CG::comalw"];
      Ver [ label="Ver" tooltip="verify antidiffusive element contributions"
            URL="\ref inciter::CG::verify"];
      OwnLim [ label="OwnLim"
               tooltip="compute limited antidiffusive element contributions"
               URL="\ref inciter::CG::lim"];
      ComLim [ label="ComLim"
               tooltip="contributions to the limited antidiffusive element
                        contributions communicated"
               URL="\ref inciter::CG::comlim"];
      Apply [ label="Apply"
              tooltip="apply limited antidiffusive element contributions"
              URL="\ref inciter::CG::limit"];
      OwnAEC -> Ver [ style="dashed" ];
      OwnALW -> Ver [ style="dashed" ];
      OwnAEC -> OwnLim [ style="solid" ];
      ComAEC -> OwnLim [ style="solid" ];
      OwnALW -> OwnLim [ style="solid" ];
      ComALW -> OwnLim [ style="solid" ];
      OwnAEC -> ComLim [ style="solid" ];
      ComAEC -> ComLim [ style="solid" ];
      OwnALW -> ComLim [ style="solid" ];
      ComALW -> ComLim [ style="solid" ];
      OwnLim -> Apply [ style="solid" ];
      ComLim -> Apply [ style="solid" ];
    }
    \enddot
    \include Inciter/cg.ci
*/
// *****************************************************************************
#ifndef CG_h
#define CG_h

#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>
#include <cstring>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <set>

#include "QuinoaConfig.h"
#include "Types.h"
#include "Fields.h"
#include "DerivedData.h"
#include "VectorReducer.h"
#include "FluxCorrector.h"
#include "Inciter/InputDeck/InputDeck.h"

#include "NoWarning/cg.decl.h"

namespace tk {
  class ExodusIIMeshWriter;
  class RootMeshWriter;
}

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! CG Charm++ chare array used to advance PDEs in time with CG+LW+FCT
class CG : public CBase_CG {

  public:
    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wunused-parameter"
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
    CG_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit CG( const CProxy_Discretization& disc,
                 const tk::CProxy_Solver& solver );

    //! Migrate constructor
    explicit CG( CkMigrateMessage* ) {}

    //! \brief Read mesh node coordinates and optionally add new edge-nodes in
    //!   case of initial uniform refinement
    void coord();

    //! \brief Setup rows, query boundary conditions, output
    //!    mesh, etc.
    void setup( tk::real v );

    //! Request owned node IDs on which a Dirichlet BC is set by the user
    void requestBCs();

    //! Look up and return old node ID for new one
    void oldID( int frompe, const std::vector< std::size_t >& newids );

    //! Compute time step size
    void dt();

    //! \brief Set ICs, compute initial time step size, output initial field
    //!   data, compute left-hand-side matrix
    void init();

    //! Update high order solution vector
    void updateSol( //solMsg* m );
                    const std::vector< std::size_t >& gid,
                    const std::vector< tk::real >& sol );

    //! Update low order solution vector
    void updateLowSol( //solMsg* m );
                       const std::vector< std::size_t >& gid,
                       const std::vector< tk::real >& sol );

    //! Advance equations to next time step
    void advance( uint64_t it, tk::real t, tk::real newdt );

    //! Receive sums of antidiffusive element contributions on chare-boundaries
    void comaec( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& P );

    //! \brief Receive contributions to the maxima and minima of unknowns of all
    //!   elements surrounding mesh nodes on chare-boundaries
    void comalw( const std::vector< std::size_t >& gid,
                    const std::vector< std::vector< tk::real > >& Q );

    //! \brief Receive contributions of limited antidiffusive element
    //!   contributions on chare-boundaries
    void comlim( const std::vector< std::size_t >& gid,
                    const std::vector< std::vector< tk::real > >& A );

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_CG::pup(p);
      p | m_itf;
      p | m_nhsol;
      p | m_nlsol;
      p | m_naec;
      p | m_nalw;
      p | m_nlim;
      p | m_disc;
      p | m_solver;
      p | m_fluxcorrector;
      p | m_u;
      p | m_ul;
      p | m_du;
      p | m_dul;
      p | m_ue;
      p | m_p;
      p | m_q;
      p | m_a;
      p | m_lhsd;
      p | m_lhso;
      p | m_pc;
      p | m_qc;
      p | m_ac;
      p | m_vol;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i CG object reference
    friend void operator|( PUP::er& p, CG& i ) { i.pup(p); }
    //@}

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

    //! Field output iteration count
    uint64_t m_itf;
    //! Counter for high order solution nodes updated
    std::size_t m_nhsol;
    //! Counter for low order solution nodes updated
    std::size_t m_nlsol;
    //! \brief Number of chares from which we received antidiffusive element
    //!   contributions on chare boundaries
    std::size_t m_naec;
    //! \brief Number of chares from which we received maximum and minimum
    //!   unknowns of elements surrounding nodes on chare boundaries
    std::size_t m_nalw;
    //! \brief Number of chares from which we received limited antidiffusion
    //!   element contributiones on chare boundaries
    std::size_t m_nlim;
    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Linear system merger and solver proxy
    tk::CProxy_Solver m_solver;
    //! \brief Map associating local node IDs to side set IDs for all side sets
    //!   read from mesh file (not only those the user sets BCs on)
    std::map< int, std::vector< std::size_t > > m_side;
    //! Flux corrector performing FCT
    FluxCorrector m_fluxcorrector;
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
    //! Flux-corrected transport data structures
    tk::Fields m_p, m_q, m_a;
    //! Sparse matrix sotring the diagonals and off-diagonals of nonzeros
    tk::Fields m_lhsd, m_lhso;
    //! Receive buffers for FCT
    std::vector< std::vector< tk::real > > m_pc, m_qc, m_ac;
    //! Total mesh volume
    tk::real m_vol;

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

    //! Search particle ina single mesh cell
    bool parinel( std::size_t p, std::size_t e, std::array< tk::real, 4 >& N );

    //! Compute and sum antidiffusive element contributions (AEC) to mesh nodes
    void aec();

    //! \brief Compute the maximum and minimum unknowns of all elements
    //!   surrounding nodes
    void alw();

    //! \brief Verify antidiffusive element contributions up to linear solver
    //!   convergence
    void verify();

    //! Compute the limited antidiffusive element contributions
    void lim();

    //! Apply limited antidiffusive element contributions
    void apply();

    //! Compute diagnostics, e.g., residuals
    void diagnostics();

    //! Verify that solution does not change at Dirichlet boundary conditions
    bool correctBC();
};

} // inciter::

#endif // CG_h
