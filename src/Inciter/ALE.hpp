// *****************************************************************************
/*!
  \file      src/Inciter/ALE.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Declarations file for distributed ALE mesh motion
  \details   Declarations file for asynchronous distributed
             arbitrary Lagrangian-Eulerian (ALE) mesh motion using Charm++.

    There are a potentially large number of ALE Charm++ chares.
    Each ALE chare gets a chunk of the full load, due to partiting the mesh.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality.
*/
// *****************************************************************************
#ifndef ALE_h
#define ALE_h

#include "ConjugateGradients.hpp"
#include "Inciter/Options/MeshVelocitySmoother.hpp"
#include "Options/UserTable.hpp"
#include "Fields.hpp"
#include "Table.hpp"

#include "NoWarning/ale.decl.h"

namespace inciter {

//! ALE Charm++ chare array used to perform arbitrary ALE mesh movement
class ALE : public CBase_ALE {

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
    ALE_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit
    ALE( const tk::CProxy_ConjugateGradients& conjugategradientsproxy,
         const std::array< std::vector< tk::real >, 3 >& coord,
         const std::vector< std::size_t >& inpoel,
         const std::vector< std::size_t >& gid,
         const std::unordered_map< std::size_t, std::size_t >& lid,
         const tk::NodeCommMap& nodecommmap );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit ALE( CkMigrateMessage* ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Solve linear system to smooth ALE mesh velocity
    void solve( CkCallback c );

    //! Query the solution of the Conjugrate Gradients linear solver
    bool converged() const;

    //! Start computing new mesh velocity for ALE mesh motion
    void start(
      const tk::UnsMesh::Coords vel,
      const std::vector< tk::real >& soundspeed,
      CkCallback done,
      const std::array< std::vector< tk::real >, 3 >& coord,
      const tk::UnsMesh::Coords coordn,
      const std::vector< tk::real >& vol0,
      const std::vector< tk::real >& vol,
      const std::unordered_map< int,
        std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& bnorm,
      std::size_t initial,
      std::size_t it,
      tk::real t,
      tk::real adt );

    //! Query the mesh velocity
    const tk::Fields& meshvel() const { return m_w; }

    //! \brief Query mesh velocity boundary conditions node lists and node list
    //!    at which ALE moves boundaries
    void meshvelBnd( const std::map< int, std::vector< std::size_t > >& bface,
                     const std::map< int, std::vector< std::size_t > >& bnode,
                     const std::vector< std::size_t >& triinpoel );

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
    void comfor( const std::vector< std::size_t >& gid,
                 const std::vector< std::array< tk::real, 3 > >& w );

    //! Apply mesh velocity smoother boundary conditions for ALE mesh motion
    void meshvelbc( tk::real maxv = 0.0 );

    //! Solve mesh velocity linear solve for ALE mesh motion
    void applied( CkDataMsg* msg = nullptr );

    //! Mesh velocity smoother linear solver converged
    void solved( CkDataMsg* msg = nullptr );

    //! Compute the gradient of the scalar potential for ALE
    void helmholtz( CkDataMsg* msg = nullptr );

    /** @name Pack/unpack (Charm++ serialization) routines */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_conjugategradients;
      p | m_done;
      p | m_soundspeed;
      p | m_nvort;
      p | m_ndiv;
      p | m_npot;
      p | m_nwf;
      p | m_nodeCommMap;
      p | m_lid;
      p | m_coord0;
      p | m_coord;
      p | m_inpoel;
      p | m_vol0;
      p | m_vol;
      p | m_it;
      p | m_t;
      p | m_adt;
      p | m_w;
      p | m_wf;
      p | m_wfc;
      p | m_veldiv;
      p | m_veldivc;
      p | m_gradpot;
      p | m_gradpotc;
      p | m_vorticity;
      p | m_vorticityc;
      p | m_bnorm;
      p | m_meshveldirbcnodes;
      p | m_meshvelsymbcnodes;
      p | m_move;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] a ALE object reference
    friend void operator|( PUP::er& p, ALE& a ) { a.pup(p); }
    ///@}

  private:
    //! Distributed conjugrate gradients solver proxy
    tk::CProxy_ConjugateGradients m_conjugategradients;
    //! Function call to continue with when mesh velocity has been computed
    CkCallback m_done;
    //! Speed of sound in mesh nodes
    std::vector< tk::real > m_soundspeed;
    //! Counter for communicating the vorticity for ALE
    std::size_t m_nvort;
    //! Counter for communicating the velocity divergence for ALE
    std::size_t m_ndiv;
    //! Counter for communicating the gradient of the scalar potential for ALE
    std::size_t m_npot;
    //! Counter for communicating the mesh force for ALE
    std::size_t m_nwf;
    //! \brief Global mesh node IDs bordering the mesh chunk held by fellow
    //!   Discretization chares associated to their chare IDs
    tk::NodeCommMap m_nodeCommMap;
    //! Local node ids associated to the global ones of owned elements
    std::unordered_map< std::size_t, std::size_t > m_lid;
    //! Mesh coordinates at the time 0 for ALE
    tk::UnsMesh::Coords m_coord0;
    //! Mesh point coordinates
    tk::UnsMesh::Coords m_coord;
    //! Element connectivity
    std::vector< std::size_t > m_inpoel;
    //! Mesh element volumes at t=t0
    std::vector< tk::real > m_vol0;
    //! Volume of nodes
    //! \details This is the volume of the mesh associated to nodes of owned
    //!   elements (sum of surrounding cell volumes / 4) with contributions from
    //!   other chares on chare-boundaries
    std::vector< tk::real > m_vol;
    //! Iteration count
    std::size_t m_it;
    //! Physics time
    tk::real m_t;
    //! alpha*dt of the Runge-Kutta time step
    tk::real m_adt;
    //! Mesh velocity for ALE mesh motion
    tk::Fields m_w;
    //! Mesh force for ALE mesh motion
    tk::Fields m_wf;
    //! Receive buffer for the mesh force for ALE
    //! \details Key: global node id, value: mesh force in nodes
    std::unordered_map< std::size_t, std::array< tk::real, 3 > > m_wfc;
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
    //! Vorticity for ALE
    tk::UnsMesh::Coords m_vorticity;
    //! Receive buffer for communication of the vorticity for ALE
    //! \details Key: global node id, value: vorticity in nodes
    std::unordered_map< std::size_t, std::array< tk::real, 3 > > m_vorticityc;
    //! Face normals in boundary points associated to side sets
    //! \details Key: local node id, value: unit normal and inverse distance
    //!   square between face centroids and points, outer key: side set id
    std::unordered_map< int,
      std::unordered_map< std::size_t, std::array< tk::real, 4 > > > m_bnorm;
    //! Unique set of nodes at which ALE mesh velocity Dirichlet BCs are set
    std::unordered_set< std::size_t > m_meshveldirbcnodes;
    //! Unique set of nodes at which ALE mesh velocity symmetry BCs are set
    std::unordered_set< std::size_t > m_meshvelsymbcnodes;
    //! Data structure storing configuration for moving boundaries with ALE
    //! \details Tuple: 0: user-defined function type (i.e., how it should be
    //!    interpreted), (1) user-defined function for, and (2) unique set of
    //!    boundary nodes for each move ... end input file block (vector).
    std::vector<
      std::tuple< tk::ctr::UserTableType,
                  tk::Table<3>,
                  std::unordered_set< std::size_t > >
    > m_move;

    //! Generate {A,x,b} for Laplacian mesh velocity smoother
    std::tuple< tk::CSR, std::vector< tk::real >, std::vector< tk::real > >
    Laplacian( std::size_t ncomp,
               const std::array< std::vector< tk::real >, 3 >& coord ) const;

    //! Initialize user-defined functions for ALE moving sides
    decltype(m_move) moveCfg();

    //! Find Dirichlet BCs on mesh velocity with prescribed movement
    bool move( std::size_t i ) const;

    //! Finalize computing fluid vorticity and velocity divergence for ALE
    void mergevel();

    //! Finalize computing the scalar potential gradient for ALE
    void gradpot();

    //! Compute mesh force for the ALE mesh velocity
    void startforce();

    //! Apply mesh force
    void meshforce();
};

} // inciter::

#endif // ALE_h
