// *****************************************************************************
/*!
  \file      src/Inciter/DG.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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
#include "Ghosts.hpp"

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
                 const CProxy_Ghosts& ghostsproxy,
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

    //! Configure Charm++ reduction types for concatenating BC nodelists
    static void registerReducers();

    //! Resize solution vectors after extension due to Ghosts and continue setup
    void resizeSolVectors();

    //! Setup: query boundary conditions, output mesh, etc.
    void setup();

    //! Receive total box IC volume and set conditions in box
    void box( tk::real v, const std::vector< tk::real >& blkvols );

    // Evaluate whether to do load balancing
    void evalLB( int nrestart );

    //! Start time stepping
    void start();

    //! Continue to next time step
    void next();

    //!  Receive chare-boundary updated solution from neighboring chares
    void comrefine( int fromch,
                    const std::vector< std::size_t >& tetid,
                    const std::vector< std::size_t >& ndof );

    //  Receive chare-boundary refined ndof from neighboring chares
    void comsmooth( int fromch,
                    const std::vector< std::size_t >& tetid,
                    const std::vector< std::size_t >& ndof );

    //! Receive chare-boundary limiter function data from neighboring chares
    void comlim( int fromch,
                 const std::vector< std::size_t >& tetid,
                 const std::vector< std::vector< tk::real > >& u,
                 const std::vector< std::vector< tk::real > >& prim );

    //! Receive chare-boundary reconstructed data from neighboring chares
    void comreco( int fromch,
                  const std::vector< std::size_t >& tetid,
                  const std::vector< std::vector< tk::real > >& u,
                  const std::vector< std::vector< tk::real > >& prim );

    //! Receive chare-boundary ghost data from neighboring chares
    void comsol( int fromch,
                 std::size_t fromstage,
                 const std::vector< std::size_t >& tetid,
                 const std::vector< std::vector< tk::real > >& u,
                 const std::vector< std::vector< tk::real > >& prim,
                 const std::vector< std::size_t >& ndof );

    //! Receive contributions to nodal gradients on chare-boundaries
    void
    comnodalExtrema( const std::vector< std::size_t >& gid,
                     const std::vector< std::vector< tk::real > >& G1,
                     const std::vector< std::vector< tk::real > >& G2 );

    //! Initialize the vector of nodal extrema
    void resizeNodalExtremac();

    //! Compute the nodal extrema of ref el derivatives for chare-boundary nodes
    void evalNodalExtrmRefEl(
      const std::size_t ncomp,
      const std::size_t nprim,
      const std::size_t ndof_NodalExtrm,
      const std::vector< std::size_t >& bndel,
      const std::vector< std::size_t >& inpoel,
      const std::vector< std::size_t >& gid,
      const std::unordered_map< std::size_t, std::size_t >& bid,
      const tk::Fields& U,
      const tk::Fields& P,
      std::vector< std::vector<tk::real> >& uNodalExtrm,
      std::vector< std::vector<tk::real> >& pNodalExtrm );

    //! \brief Receive nodal solution (ofor field output) contributions from
    //!   neighboring chares
    void comnodeout( const std::vector< std::size_t >& gid,
                     const std::vector< std::size_t >& nesup,
                     const std::vector< std::vector< tk::real > >& Lu,
                     const std::vector< std::vector< tk::real > >& Lp );

    //! Optionally refine/derefine mesh
    void refine( const std::vector< tk::real >& l2res );

    //! Receive new mesh from Refiner
    void resizePostAMR(
      const std::vector< std::size_t >& /* ginpoel */,
      const tk::UnsMesh::Chunk& chunk,
      const tk::UnsMesh::Coords& coord,
      const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& /* addedNodes */,
      const std::unordered_map< std::size_t, std::size_t >& addedTets,
      const std::set< std::size_t >& removedNodes,
      const std::unordered_map< std::size_t, std::size_t >& amrNodeMap,
      const tk::NodeCommMap& nodeCommMap,
      const std::map< int, std::vector< std::size_t > >& bface,
      const std::map< int, std::vector< std::size_t > >& /* bnode */,
      const std::vector< std::size_t >& triinpoel,
      const std::unordered_map< std::size_t, std::set< std::size_t > >&
        elemblockid );

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

    //! Unused in DG
    void resized() {}

    //! (no-op)
    void transferSol() {}

    //! (no-op)
    void advance( tk::real, tk::real ) {}

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
      p | m_ghosts;
      p | m_ndof_NodalExtrm;
      p | m_nsol;
      p | m_ninitsol;
      p | m_nlim;
      p | m_nnod;
      p | m_nrefine;
      p | m_nsmooth;
      p | m_nreco;
      p | m_nnodalExtrema;
      p | m_u;
      p | m_un;
      p | m_p;
      p | m_geoElem;
      p | m_lhs;
      p | m_mtInv;
      p | m_uNodalExtrm;
      p | m_pNodalExtrm;
      p | m_uNodalExtrmc;
      p | m_pNodalExtrmc;
      p | m_rhs;
      p | m_rhsprev;
      p | m_stiffrhs;
      p | m_stiffrhsprev;
      p | m_stiffEqIdx;
      p | m_nonStiffEqIdx;
      p | m_nstiffeq;
      p | m_nnonstiffeq;
      p | m_npoin;
      p | m_diag;
      p | m_stage;
      p | m_ndof;
      p | m_numEqDof;
      p | m_uc;
      p | m_pc;
      p | m_ndofc;
      p | m_initial;
      p | m_uElemfields;
      p | m_pElemfields;
      p | m_uNodefields;
      p | m_pNodefields;
      p | m_uNodefieldsc;
      p | m_pNodefieldsc;
      p | m_outmesh;
      p | m_boxelems;
      p | m_shockmarker;
      p | m_rho0mat;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i DG object reference
    friend void operator|( PUP::er& p, DG& i ) { i.pup(p); }
    //@}

  private:
    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Distributed Ghosts proxy
    CProxy_Ghosts m_ghosts;
    //! \brief Degree of freedom for nodal extrema vector. When DGP1 is applied,
    //!   there is one degree of freedom for cell average variable. When DGP2 is
    //!   applied, the degree of freedom is 4 which refers to cell average and
    //!   gradients in three directions
    std::size_t m_ndof_NodalExtrm;
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
    //! Counter signaling that we have received all refined ndof
    std::size_t m_nrefine;
    //! Counter signaling that we have received all smoothed ndof
    std::size_t m_nsmooth;
    //! Counter signaling that we have received all our reconstructed ghost data
    std::size_t m_nreco;
    //! \brief Counter signaling that we have received all our nodal extrema from
    //!   ghost chare partitions
    std::size_t m_nnodalExtrema;
    //! Counters signaling how many stiff and non-stiff equations in the system
    std::size_t m_nstiffeq, m_nnonstiffeq;
    //! Vector of unknown/solution average over each mesh element
    tk::Fields m_u;
    //! Vector of unknown at previous time-step
    tk::Fields m_un;
    //! Vector of primitive quantities over each mesh element
    tk::Fields m_p;
    //! Element geometry
    tk::Fields m_geoElem;
    //! Left-hand side mass-matrix which is a diagonal matrix
    tk::Fields m_lhs;
    //! Vector of right-hand side
    tk::Fields m_rhs;
    //! Vector of previous right-hand side values used in the IMEX-RK scheme
    tk::Fields m_rhsprev;
    //! Vector of right-hand side for stiff equations
    tk::Fields m_stiffrhs;
    //! Vector of previous right-hand side values for stiff equations
    tk::Fields m_stiffrhsprev;
    //! Vectors that indicates which equations are stiff and non-stiff
    std::vector< std::size_t > m_stiffEqIdx, m_nonStiffEqIdx;
    //! Inverse of Taylor mass-matrix
    std::vector< std::vector< tk::real > > m_mtInv;
    //! Vector of nodal extrema for conservative variables
    std::vector< std::vector<tk::real> > m_uNodalExtrm;
    //! Vector of nodal extrema for primitive variables
    std::vector< std::vector<tk::real> > m_pNodalExtrm;
    //! Buffer for vector of nodal extrema for conservative variables
    std::unordered_map< std::size_t, std::vector< tk::real > > m_uNodalExtrmc;
    //! Buffer for vector of nodal extrema for primitive variables
    std::unordered_map< std::size_t, std::vector< tk::real > > m_pNodalExtrmc;
    //! Counter for number of nodes on this chare excluding ghosts
    std::size_t m_npoin;
    //! Diagnostics object
    ElemDiagnostics m_diag;
    //! Runge-Kutta stage counter
    std::size_t m_stage;
    //! Vector of local number of degrees of freedom for each element
    std::vector< std::size_t > m_ndof;
    //! Vector of number of degrees of freedom for each PDE equation/component
    std::vector< std::size_t > m_numEqDof;
    //! Solution receive buffers for ghosts only
    std::array< std::vector< std::vector< tk::real > >, 3 > m_uc;
    //! Primitive-variable receive buffers for ghosts only
    std::array< std::vector< std::vector< tk::real > >, 3 > m_pc;
    //! \brief Number of degrees of freedom (for p-adaptive) receive buffers
    //!   for ghosts only
    std::array< std::vector< std::size_t >, 3 > m_ndofc;
    //! 1 if starting time stepping, 0 if during time stepping
    std::size_t m_initial;
    //! Solution elem output fields
    tk::Fields m_uElemfields;
    //! Primitive elem output fields
    tk::Fields m_pElemfields;
    //! Solution nodal output fields
    tk::Fields m_uNodefields;
    //! Primitive nodal output fields
    tk::Fields m_pNodefields;
    //! Receive buffer for communication of solution node fields
    //! \details Key: global node id, value: output fields and number of
    //!   elements surrounding the node
    std::unordered_map< std::size_t, std::pair< std::vector< tk::real >,
                                                std::size_t > > m_uNodefieldsc;
    //! Receive buffer for communication of primitive quantity node fields
    //! \details Key: global node id, value: output fields and number of
    //!   elements surrounding the node
    std::unordered_map< std::size_t, std::pair< std::vector< tk::real >,
                                                std::size_t > > m_pNodefieldsc;
    //! Storage for refined mesh used for field output
    Ghosts::OutMesh m_outmesh;
    //! Element ids at which box ICs are defined by user (multiple boxes)
    std::vector< std::unordered_set< std::size_t > > m_boxelems;
    //! Shock detection marker for field output
    std::vector< std::size_t > m_shockmarker;
    //! Initial densities of each material
    std::vector< tk::real > m_rho0mat;

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }

    //! Access bound Discretization class pointer
    Ghosts* myGhosts() const {
      Assert( m_ghosts[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_ghosts[ thisIndex ].ckLocal();
    }

    //! Output mesh field data
    void writeFields(
      CkCallback c,
      const std::unordered_map< std::size_t, std::size_t >& addedTets );

    //! Add the protective layer for ndof refinement
    void refine();

    //! Smooth the refined ndof distribution
    void smooth();

    //! Compute solution reconstructions
    void reco();

    //! Compute nodal extrema at chare-boundary nodes
    void nodalExtrema();

    //! Compute limiter function
    void lim();

    //! Compute time step size
    void dt();

    //! Evaluate whether to continue with next time step stage
    void stage();

    //! Evaluate whether to save checkpoint/restart
    void evalRestart();

    //! p-refine all elements that are adjacent to p-refined elements
    void refine_ndof();

    //! Smooth the refined ndof distribution to avoid zigzag refinement
    void smooth_ndof();

    //! Decide wether to output field data
    bool fieldOutput() const;

    //! Decide if we write field output using a refined mesh
    bool refinedOutput() const;

    //! Start preparing fields for output to file
    void startFieldOutput( CkCallback c );

    //! Compute the integration step for IMEX-RK
    void imex_integrate();
};

} // inciter::

#endif // DG_h
