// *****************************************************************************
/*!
  \file      src/Inciter/OversetFE.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     OversetFE for a PDE system with continuous Galerkin FE + RK
  \details   OversetFE advances a system of partial differential equations
    using a continuous Galerkin (CG) finite element (FE) spatial discretization
    (using linear shapefunctions on tetrahedron elements) combined with a
    Runge-Kutta (RK) time stepping scheme and overset grids.

    There are a potentially large number of OversetFE Charm++ chares created by
    Transporter. Each OversetFE gets a chunk of the full load (part of the mesh)
    and does the same: initializes and advances a number of PDE systems in time.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality.
*/
// *****************************************************************************
#ifndef OversetFE_h
#define OversetFE_h

#include <vector>
#include <map>

#include "Types.hpp"
#include "Fields.hpp"
#include "Table.hpp"
#include "DerivedData.hpp"
#include "FluxCorrector.hpp"
#include "NodeDiagnostics.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Ghosts.hpp"

#include "NoWarning/oversetfe.decl.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! OversetFE Charm++ chare array used to advance PDEs in time with OversetFE+RK
class OversetFE : public CBase_OversetFE {

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
    OversetFE_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit OversetFE( const CProxy_Discretization& disc,
      const CProxy_Ghosts& ghostsproxy,
      const std::map< int, std::vector< std::size_t > >& bface,
      const std::map< int, std::vector< std::size_t > >& bnode,
      const std::vector< std::size_t >& triinpoel );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit OversetFE( CkMigrateMessage* msg ) : CBase_OversetFE( msg ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ custom reduction types initiated from this chare array
    static void registerReducers();

    //! Return from migration
    void ResumeFromSync() override;

    //! Start setup for solution
    void setup();

    //! Receive total box IC volume and set conditions in box
    void box( tk::real v, const std::vector< tk::real >& blkvols );

    // Start time stepping
    void start();

    //! Advance equations to next time step
    void advance( tk::real newdt );

    //! Compute left-hand side of transport equations
    void lhs();

    //! Transfer solution from O to B
    void transferOtoB();

    //! Receive contributions to duual-face normals on chare boundaries
    void comdfnorm(
      const std::unordered_map< tk::UnsMesh::Edge,
              std::array< tk::real, 3 >,
              tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> >& dfnorm );

    //! Receive boundary point normals on chare-boundaries
    void comnorm( const std::unordered_map< int,
      std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& innorm );

    //! Receive mesh block information for nodes on chare-boundaries
    void comblk( const std::vector< std::size_t >& gid,
      const std::vector< std::set< std::size_t > >& mb );

    //! Receive contributions to gradients on chare-boundaries
    void comChBndGrad( const std::vector< std::size_t >& gid,
                       const std::vector< std::vector< tk::real > >& G );

    //! Receive contributions to right-hand side vector on chare-boundaries
    void comrhs( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& R );

    //! Update solution at the end of time step
    void update( const tk::Fields& a );

    //! Optionally refine/derefine mesh
    void refine( const std::vector< tk::real >& l2res );

    //! Extract field output to file
    void extractFieldOutput(
      const std::vector< std::size_t >& /* ginpoel */,
      const tk::UnsMesh::Chunk& /*chunk*/,
      const tk::UnsMesh::Coords& /*coord*/,
      const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& /* addedNodes */,
      const std::unordered_map< std::size_t, std::size_t >& /*addedTets*/,
      const tk::NodeCommMap& /*nodeCommMap*/,
      const std::map< int, std::vector< std::size_t > >& /*bface*/,
      const std::map< int, std::vector< std::size_t > >& /* bnode */,
      const std::vector< std::size_t >& /*triinpoel*/,
      CkCallback /*c*/ ) {}

    //! Const-ref access to current solution
    //! \return Const-ref to current solution
    const tk::Fields& solution() const { return m_u; }

    //! Evaluate whether to continue with next time step
    void step();

    // Evaluate whether to do load balancing
    void evalLB( int nrestart );

    //! Evaluate whether to continue with next time step stage
    void stage();

    //! Continue to next time step
    void next();

    //! Size communication buffers (no-op)
    void resizeComm() {}

    //! Setup node-neighborhood (no-op)
    void nodeNeighSetup() {}

    //! Receive new mesh from Refiner (no-op)
    void resizePostAMR(
      const std::vector< std::size_t >&,
      const tk::UnsMesh::Chunk&,
      const tk::UnsMesh::Coords&,
      const std::unordered_map< std::size_t, tk::UnsMesh::Edge >&,
      const std::unordered_map< std::size_t, std::size_t >&,
      const std::set< std::size_t >&,
      const std::unordered_map< std::size_t, std::size_t >&,
      const tk::NodeCommMap&,
      const std::map< int, std::vector< std::size_t > >&,
      const std::map< int, std::vector< std::size_t > >&,
      const std::vector< std::size_t >&,
      const std::unordered_map< std::size_t, std::set< std::size_t > >& ) {}

    //! Resizing data structures after mesh refinement has completed (no-op)
    void resized() {}

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_disc;
      p | m_nsol;
      p | m_ngrad;
      p | m_nrhs;
      p | m_nbnorm;
      p | m_ndfnorm;
      p | m_nmblk;
      p | m_bnode;
      p | m_bface;
      p | m_triinpoel;
      p | m_bndel;
      p | m_dfnorm;
      p | m_dfnormc;
      p | m_dfn;
      p | m_esup;
      p | m_psup;
      p | m_u;
      p | m_uc;
      p | m_un;
      p | m_rhs;
      p | m_rhsc;
      p | m_chBndGrad;
      p | m_dirbc;
      p | m_chBndGradc;
      p | m_blank;
      p | m_diag;
      p | m_bnorm;
      p | m_bnormc;
      p | m_symbcnodes;
      p | m_farfieldbcnodes;
      p | m_symbctri;
      p | m_timedepbcnodes;
      p | m_timedepbcFn;
      p | m_stage;
      p | m_boxnodes;
      p | m_edgenode;
      p | m_edgeid;
      p | m_dtp;
      p | m_tp;
      p | m_finished;
      p | m_movedmesh;
      p | m_nusermeshblk;
      p | m_nodeblockid;
      p | m_nodeblockidc;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i OversetFE object reference
    friend void operator|( PUP::er& p, OversetFE& i ) { i.pup(p); }
    //@}

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Counter for high order solution vector nodes updated
    std::size_t m_nsol;
    //! Counter for nodal gradients updated
    std::size_t m_ngrad;
    //! Counter for right-hand side vector nodes updated
    std::size_t m_nrhs;
    //! Counter for receiving boundary point normals
    std::size_t m_nbnorm;
    //! Counter for receiving dual-face normals on chare-boundary edges
    std::size_t m_ndfnorm;
    //! Counter for mesh block nodes updated
    std::size_t m_nmblk;
    //! Boundary node lists mapped to side set ids used in the input file
    std::map< int, std::vector< std::size_t > > m_bnode;
    //! Boundary face lists mapped to side set ids used in the input file
    std::map< int, std::vector< std::size_t > > m_bface;
    //! Boundary triangle face connecitivity where BCs are set by user
    std::vector< std::size_t > m_triinpoel;
    //! Elements along mesh boundary
    std::vector< std::size_t > m_bndel;
    //! Dual-face normals along edges
    std::unordered_map< tk::UnsMesh::Edge, std::array< tk::real, 3 >,
                        tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > m_dfnorm;
    //! Receive buffer for dual-face normals along chare-boundary edges
    std::unordered_map< tk::UnsMesh::Edge, std::array< tk::real, 3 >,
                     tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > m_dfnormc;
    //! Streamable dual-face normals
    std::vector< tk::real > m_dfn;
    //! El;ements surrounding points
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > m_esup;
    //! Points surrounding points
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > m_psup;
    //! Unknown/solution vector at mesh nodes
    tk::Fields m_u;
    //! \brief Copy of unknown/solution vector at mesh nodes for m2m transfer,
    //!   appended with a solution-transfer-flag. This flag indicates
    //!   appropriate solution transfers for the overset procedure.
    //!   Value 0: Do not transfer solution, 1: transfer sol, 2: blank nodes
    //! TODO: avoid creating this copy
    tk::Fields m_uc;
    //! Unknown/solution vector at mesh nodes at previous time
    tk::Fields m_un;
    //! Right-hand side vector (for the high order system)
    tk::Fields m_rhs;
    //! Receive buffer for communication of the right hand side
    //! \details Key: global node id, value: rhs for all scalar components per
    //!   node.
    std::unordered_map< std::size_t, std::vector< tk::real > > m_rhsc;
    //! Nodal gradients at chare-boundary nodes
    tk::Fields m_chBndGrad;
    //! Boundary conditions evaluated and assigned to local mesh node IDs
    //! \details Vector of pairs of bool and boundary condition value associated
    //!   to local mesh node IDs at which the user has set Dirichlet boundary
    //!   conditions for all PDEs integrated. The bool indicates whether the BC
    //!   is set at the node for that component the if true, the real value is
    //!   the increment (from t to dt) in the BC specified for a component.
    std::unordered_map< std::size_t,
      std::vector< std::pair< bool, tk::real > > > m_dirbc;
    //! Receive buffer for communication of the nodal gradients
    //! \details Key: chare id, value: gradients for all scalar components per
    //!   node
    std::unordered_map< std::size_t, std::vector< tk::real > > m_chBndGradc;
    //! Blanking coefficient for overset, indicating hole in the background mesh
    std::vector< tk::real > m_blank;
    //! Diagnostics object
    NodeDiagnostics m_diag;
    //! Face normals in boundary points associated to side sets
    //! \details Key: local node id, value: unit normal and inverse distance
    //!   square between face centroids and points, outer key: side set id
    std::unordered_map< int,
      std::unordered_map< std::size_t, std::array< tk::real, 4 > > > m_bnorm;
    //! \brief Receive buffer for communication of the boundary point normals
    //!   associated to side sets
    //! \details Key: global node id, value: normals (first 3 components),
    //!   inverse distance squared (4th component), outer key, side set id
    std::unordered_map< int,
      std::unordered_map< std::size_t, std::array< tk::real, 4 > > > m_bnormc;
    //! Unique set of nodes at which symmetry BCs are set
    std::unordered_set< std::size_t > m_symbcnodes;
    //! Unique set of nodes at which farfield BCs are set
    std::unordered_set< std::size_t > m_farfieldbcnodes;
    //! Vector with 1 at symmetry BC boundary triangles
    std::vector< int > m_symbctri;
    //! \brief Unique set of nodes at which time dependent BCs are set
    //    for each time dependent BC
    std::vector< std::unordered_set< std::size_t > > m_timedepbcnodes;
    //! \brief User defined discrete function of time used in the time dependent
    //    BCs associated with (index in vector) the number of distinct time
    //    dependent BCs specified. This index is the same as the index in
    //    m_timedepbcnodes.
    std::vector< tk::Table<5> > m_timedepbcFn;
    //! Runge-Kutta stage counter
    std::size_t m_stage;
    //! Mesh node ids at which user-defined box ICs are defined (multiple boxes)
    std::vector< std::unordered_set< std::size_t > > m_boxnodes;
    //! Local node IDs of edges
    std::vector< std::size_t > m_edgenode;
    //! Edge ids in the order of access
    std::vector< std::size_t > m_edgeid;
    //! Time step size for each mesh node
    std::vector< tk::real > m_dtp;
    //! Physical time for each mesh node
    std::vector< tk::real > m_tp;
    //! True in the last time step
    int m_finished;
    //! True if overset mesh moved
    int m_movedmesh;
    //! Number of mesh-blocks with user-defined ICs
    std::size_t m_nusermeshblk;
    //! Local node ids associated with mesh block ids
    std::unordered_map< std::size_t, std::set< std::size_t > > m_nodeblockid;
    //! Receive buffer for communication of the mesh block ids
    //! \details Key: mesh block id, value: set of global node ids for nodes
    //!   in this mesh block.
    std::unordered_map< std::size_t, std::set< std::size_t > > m_nodeblockidc;

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }

    //! Compute normal of dual-mesh associated to edge
    std::array< tk::real, 3 >
    edfnorm( const tk::UnsMesh::Edge& edge,
             const std::unordered_map< tk::UnsMesh::Edge,
                     std::vector< std::size_t >,
                     tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> >& esued ) const;

    //! Compute chare-boundary edges
    void bndEdges();

    //! Start (re-)computing boundare point-, and dual-face normals
    void norm();

    //! Compute dual-face normals associated to edges
    void dfnorm();

    //! Compute boundary point normals
    void
    bnorm( const std::unordered_map< int,
             std::unordered_set< std::size_t > >& bcnodes );

    //! Set flags informing solution transfer decisions
    void setTransferFlags(
      std::size_t dirn );

    //! \brief Finish computing dual-face and boundary point normals and apply
    //!   boundary conditions on the initial conditions
    void normfinal();

    //! Continue setup for solution, after communication for mesh blocks
    void continueSetup();

    //! Output mesh field data and continue to next time step
    void out();

    //! Output mesh-based fields to file
    void writeFields( CkCallback c );

    //! Combine own and communicated contributions to normals
    void mergelhs();

    //! Compute gradients
    void chBndGrad();

    //! Compute righ-hand side vector of transport equations
    void rhs();

    //! Advance systems of equations
    void solve();

    //! Compute time step size
    void dt();

    //! Transfer solution to other solver and mesh if coupled
    void transfer();

    //! Evaluate whether to save checkpoint/restart
    void evalRestart();

    //! Query/update boundary-conditions-related data structures from user input
    void getBCNodes();

    // \brief Apply the transferred solution to the solution vector based on
    //   transfer flags previously set up
    void applySolTransfer( std::size_t dirn );

    //! Apply boundary conditions
    void BC();
};

} // inciter::

#endif // OversetFE_h
