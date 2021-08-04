// *****************************************************************************
/*!
  \file      src/Inciter/ALECG.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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

    ALE time-stepping is performed in an unsplit fashion, as opposed to
    Lagrange + remap. See also J. Waltz, N.R. Morgan, T.R. Canfield, M.R.J.
    Charest, L.D. Risinger, J.G. Wohlbier, A three-dimensional finite element
    arbitrary Lagrangian–Eulerian method for shock hydrodynamics on unstructured
    grids, Computers & Fluids, 92: 172-187, 2014.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality.
*/
// *****************************************************************************
#ifndef ALECG_h
#define ALECG_h

#include <vector>
#include <map>

#include "Types.hpp"
#include "Fields.hpp"
#include "DerivedData.hpp"
#include "FluxCorrector.hpp"
#include "NodeDiagnostics.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

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
    explicit ALECG( const CProxy_Discretization& disc,
                    const std::map< int, std::vector< std::size_t > >& bface,
                    const std::map< int, std::vector< std::size_t > >& bnode,
                    const std::vector< std::size_t >& triinpoel );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit ALECG( CkMigrateMessage* msg ) : CBase_ALECG( msg ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ custom reduction types initiated from this chare array
    static void registerReducers();

    //! Return from migration
    void ResumeFromSync() override;

    //! Size communication buffers (no-op)
    void resizeComm() {}

    //! Setup node-neighborhood (no-op)
    void nodeNeighSetup() {}

    //! Start setup for solution
    void setup();

    //! Receive total box IC volume and set conditions in box
    void box( tk::real v );

    // Start time stepping
    void start();

    //! Advance equations to next time step
    void advance( tk::real newdt );

    //! Apply mesh velocity smoother boundary conditions for ALE mesh motion
    void meshvelbc( tk::real maxv = 0.0 );

    //! Compute left-hand side of transport equations
    void lhs();

    //! Receive contributions to duual-face normals on chare boundaries
    void comdfnorm(
      const std::unordered_map< tk::UnsMesh::Edge,
              std::array< tk::real, 3 >,
              tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> >& dfnorm );

    //! Receive boundary point normals on chare-boundaries
    void comnorm( const std::unordered_map< int,
      std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& innorm );

    //! Receive contributions to gradients on chare-boundaries
    void comChBndGrad( const std::vector< std::size_t >& gid,
                       const std::vector< std::vector< tk::real > >& G );

    //! Receive contributions to right-hand side vector on chare-boundaries
    void comrhs( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& R );

    //! Receive contributions to vorticity on chare-boundaries
    void comvort( const std::vector< std::size_t >& gid,
                  const std::vector< std::array< tk::real, 3 > >& v );

    //! Receive contributions to velocity divergence on chare-boundaries
    void comdiv( const std::vector< std::size_t >& gid,
                 const std::vector< tk::real >& v );

    //! Receive contributions to scalar potential gradient on chare-boundaries
    void compot( const std::vector< std::size_t >& gid,
                 const std::vector< std::array< tk::real, 3 > >& v );

    //! Receive contributions to ALE mesh force on chare-boundaries
    void commeshforce( const std::vector< std::size_t >& gid,
                       const std::vector< std::array< tk::real, 3 > >& w );

    //! Update solution at the end of time step
    void update( const tk::Fields& a );

    //! Optionally refine/derefine mesh
    void refine( const std::vector< tk::real >& l2res );

    //! Receive new mesh from Refiner
    void resizePostAMR(
      const std::vector< std::size_t >& ginpoel,
      const tk::UnsMesh::Chunk& chunk,
      const tk::UnsMesh::Coords& coord,
      const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& addedNodes,
      const std::unordered_map< std::size_t, std::size_t >& addedTets,
      const std::set< std::size_t >& removedNodes,
      const tk::NodeCommMap& nodeCommMap,
      const std::map< int, std::vector< std::size_t > >& bface,
      const std::map< int, std::vector< std::size_t > >& bnode,
      const std::vector< std::size_t >& triinpoel );

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

    //! Resizing data sutrctures after mesh refinement has been completed
    void resized();

    //! Solve mesh velocity linear solve for ALE mesh motion
    void applied( CkDataMsg* msg = nullptr );

    //! Mesh velocity smoother linear solver converged
    void smoothed( CkDataMsg* msg = nullptr );

    //! Compute the gradient of the scalar potential for ALE
    void helmholtz( CkDataMsg* msg = nullptr );

    //! Evaluate whether to continue with next time step
    void step();

    // Evaluate whether to do load balancing
    void evalLB( int nrestart );

    //! Evaluate whether to continue with next time step stage
    void stage();

    //! Continue to next time step
    void next();

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_disc;
      p | m_initial;
      p | m_nsol;
      p | m_ngrad;
      p | m_nrhs;
      p | m_nvort;
      p | m_ndiv;
      p | m_npot;
      p | m_nwf;
      p | m_nbnorm;
      p | m_ndfnorm;
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
      p | m_un;
      p | m_w;
      p | m_vel;
      p | mp;
      p | m_wf;
      p | m_veldiv;
      p | m_veldivc;
      p | m_gradpot;
      p | m_gradpotc;
      p | m_wfc;
      p | m_rhs;
      p | m_rhsc;
      p | m_chBndGrad;
      p | m_dirbc;
      p | m_chBndGradc;
      p | m_diag;
      p | m_bnorm;
      p | m_bnormn;
      p | m_bnormc;
      p | m_symbcnodes;
      p | m_farfieldbcnodes;
      p | m_meshveldirbcnodes;
      p | m_meshvelsymbcnodes;
      p | m_symbctri;
      p | m_spongenodes;
      p | m_stage;
      p | m_boxnodes;
      p | m_edgenode;
      p | m_edgeid;
      p | m_dtp;
      p | m_tp;
      p | m_finished;
      p | m_newmesh;
      p | m_coordn;
      p | m_vorticity;
      p | m_vorticityc;
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
    //! 1 if starting time stepping, 0 if during time stepping
    std::size_t m_initial;
    //! Counter for high order solution vector nodes updated
    std::size_t m_nsol;
    //! Counter for nodal gradients updated
    std::size_t m_ngrad;
    //! Counter for right-hand side vector nodes updated
    std::size_t m_nrhs;
    //! Counter for communicating the vorticity for ALE
    std::size_t m_nvort;
    //! Counter for communicating the velocity divergence for ALE
    std::size_t m_ndiv;
    //! Counter for communicating the gradient of the scalar potential for ALE
    std::size_t m_npot;
    //! Counter for communicating the mesh force for ALE
    std::size_t m_nwf;
    //! Counter for receiving boundary point normals
    std::size_t m_nbnorm;
    //! Counter for receiving dual-face normals on chare-boundary edges
    std::size_t m_ndfnorm;
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
    //! Unknown/solution vector at mesh nodes at previous time
    tk::Fields m_un;
    //! Mesh velocity for ALE mesh motion
    tk::Fields m_w;
    //! Mesh force for ALE mesh motion
    tk::Fields m_wf;
    std::vector< tk::real > mp;
    //! Fluid velocity for ALE mesh motion
    tk::UnsMesh::Coords m_vel;
    //! Fluid velocity divergence for ALE mesh motion
    std::vector< tk::real > m_veldiv;
    //! Receive buffer for communication of the velocity divergence for ALE
    //! \details Key: global node id, value: divergence in nodes
    std::unordered_map< std::size_t, tk::real > m_veldivc;
    //! Gradient of the scalar potentinal for ALE mesh motion
    tk::UnsMesh::Coords m_gradpot;
    //! Receive buffer for the gradient of the scalar potential for ALE
    //! \details Key: global node id, value: scalar potential gradient in nodes
    std::unordered_map< std::size_t, std::array< tk::real, 3 > > m_gradpotc;
    //! Receive buffer for the mesh force for ALE
    //! \details Key: global node id, value: mesh force in nodes
    std::unordered_map< std::size_t, std::array< tk::real, 3 > > m_wfc;
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
    //! Diagnostics object
    NodeDiagnostics m_diag;
    //! Face normals in boundary points associated to side sets
    //! \details Key: local node id, value: unit normal and inverse distance
    //!   square between face centroids and points, outer key: side set id
    std::unordered_map< int,
      std::unordered_map< std::size_t, std::array< tk::real, 4 > > > m_bnorm;
    //! Face normals in boundary points at time n for ALE
    decltype(m_bnorm) m_bnormn;
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
    //! Unique set of nodes at which ALE mesh velocity Dirichlet BCs are set
    std::unordered_set< std::size_t > m_meshveldirbcnodes;
    //! Unique set of nodes at which ALE mesh velocity symmetry BCs are set
    std::unordered_set< std::size_t > m_meshvelsymbcnodes;
    //! Vector with 1 at symmetry BC boundary triangles
    std::vector< int > m_symbctri;
    //! Unique set of nodes at which sponge parameters are set
    std::unordered_set< std::size_t > m_spongenodes;
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
    //! State indicating the reason we are recomputing the normals
    int m_newmesh;
    //! Mesh coordinates at the time n for ALE
    tk::UnsMesh::Coords m_coordn;
    //! Vorticity for ALE
    tk::UnsMesh::Coords m_vorticity;
    //! Receive buffer for communication of the vorticity for ALE
    //! \details Key: global node id, value: vorticity in nodes
    std::unordered_map< std::size_t, std::array< tk::real, 3 > > m_vorticityc;

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

    //! \brief Finish computing dual-face and boundary point normals and apply
    //!   boundary conditions on the initial conditions
    void normfinal();

    //! Output mesh and particle fields to files
    void out();

    //! Output mesh-based fields to file
    void writeFields( CkCallback c );

    //! Combine own and communicated contributions to normals
    void merge();

    //! Compute gradients
    void chBndGrad();

    //! Start computing new mesh veloctity for ALE mesh motion
    void meshvel();

    //! Assign new mesh veloctity for ALE mesh motion
    void assignMeshvel();

    //! Compute righ-hand side vector of transport equations
    void rhs();

    //! Advance systems of equations
    void solve();

    //! Continue after ALE mesh movement
    void ale();

    //! Compute time step size
    void dt();

    //! Transfer solution to other solver and mesh if coupled
    void transfer();

    //! Evaluate whether to save checkpoint/restart
    void evalRestart();

    //! Query boundary conditions from user input
    void queryBC();

    //! Apply boundary conditions
    void BC();

    //! Multiply solution with mesh volume
    void volumetric( tk::Fields& u );

    //! Divide solution with mesh volume
    void conserved( tk::Fields& u );

    //! Finalize computing fluid vorticity for ALE
    void vorticity();

    //! Finalize computing the velocity divergence for ALE
    void veldiv();

    //! Finalize computing the scalar potential gradient for ALE
    void gradpot();

    //! Apply mesh force
    void meshforce();
};

} // inciter::

#endif // ALECG_h
