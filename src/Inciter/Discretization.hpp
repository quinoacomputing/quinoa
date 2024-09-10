// *****************************************************************************
/*!
  \file      src/Inciter/Discretization.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \details   Data and functionality common to all discretization schemes
     The Discretization class contains data and functionality common to all
     discretization schemes.
*/
// *****************************************************************************
#ifndef Discretization_h
#define Discretization_h

#include <brigand/algorithms/for_each.hpp>

#include "Types.hpp"
#include "Timer.hpp"
#include "Fields.hpp"
#include "PUPUtil.hpp"
#include "PDFReducer.hpp"
#include "UnsMesh.hpp"
#include "CommMap.hpp"
#include "History.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "M2MTransfer.hpp"

#include "NoWarning/discretization.decl.h"
#include "NoWarning/refiner.decl.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! \brief Discretization Charm++ chare array holding common functinoality to
//!   all discretization schemes
class Discretization : public CBase_Discretization {

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
    Discretization_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit
      Discretization(
        std::size_t meshid,
        const std::vector< CProxy_Discretization >& disc,
        const CProxy_ALE& aleproxy,
        const tk::CProxy_ConjugateGradients& conjugategradientsproxy,
        const CProxy_Transporter& transporter,
        const tk::CProxy_MeshWriter& meshwriter,
        const tk::UnsMesh::CoordMap& coordmap,
        const tk::UnsMesh::Chunk& el,
        const tk::CommMaps& msum,
        const std::map< int, std::vector< std::size_t > >& bface,
        const std::vector< std::size_t >& triinpoel,
        const std::unordered_map< std::size_t, std::set< std::size_t > >&
          elemblockid,
        int nc );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit Discretization( CkMigrateMessage* ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ reduction types
    static void registerReducers();

    //! Register mesh with mesh-transfer lib
    void addMesh();

    //! Start computing new mesh veloctity for ALE mesh motion
    void meshvelStart(
      const tk::UnsMesh::Coords vel,
      const std::vector< tk::real >& soundspeed,
      const std::unordered_map< int,
        std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& bnorm,
      tk::real adt,
      CkCallback done ) const;

    //! Query the mesh velocity
    const tk::Fields& meshvel() const;

   //! \brief Query ALE mesh velocity boundary condition node lists and node
   //!   lists at which ALE moves boundaries
   void meshvelBnd(
     const std::map< int, std::vector< std::size_t > >& bface,
     const std::map< int, std::vector< std::size_t > >& bnode,
     const std::vector< std::size_t >& triinpoel ) const;

    //! Assess and record mesh velocity linear solver convergence
    void meshvelConv();

    //! \brief Our mesh has been registered with the mesh-to-mesh transfer
    //!   library (if coupled to other solver)
    void transferInit();

    //! Finish setting up communication maps and solution transfer callbacks
    void comfinal();

    //! Start solution transfer (if coupled)
    void transfer(
      tk::Fields& u,
      std::size_t dirn,
      CkCallback cb );

    //! Solution transfer from background to overset mesh completed (from ExaM2M)
    void to_complete();

    //! Solution transfer from overset to background mesh completed (from ExaM2M)
    void from_complete();

    //! Solution transfer completed (from dest Discretization)
    void transfer_complete_from_dest();

    //! Solution transfer completed (one-way)
    void transfer_complete();

    //! Resize mesh data structures after mesh refinement
    void resizePostAMR(
      const tk::UnsMesh::Chunk& chunk,
      const tk::UnsMesh::Coords& coord,
      const std::unordered_map< std::size_t, std::size_t >& amrNodeMap,
      const tk::NodeCommMap& nodeCommMap,
      const std::set< std::size_t >& removedNodes,
      const std::unordered_map< std::size_t, std::set< std::size_t > >&
        elemblockid );

    //! Get ready for (re-)computing/communicating nodal volumes
    void startvol();

    //! Sum mesh volumes to nodes, start communicating them on chare-boundaries
    void vol();

    //! Set Refiner Charm++ proxy
    void setRefiner( const CProxy_Refiner& ref );

    //! Collect nodal volumes across chare boundaries
    void comvol( const std::vector< std::size_t >& gid,
                 const std::vector< tk::real >& nodevol );

    //! Sum mesh volumes and contribute own mesh volume to total volume
    void totalvol();

    //! Compute mesh cell statistics
    void stat( tk::real mesh_volume );

    //! Compute total box IC volume
    void
    boxvol( const std::vector< std::unordered_set< std::size_t > >& nodes,
      const std::unordered_map< std::size_t, std::set< std::size_t > >& nodeblk,
      std::size_t nuserblk );

    /** @name Accessors */
    ///@{
    //! Coordinates accessor as const-ref
    const tk::UnsMesh::Coords& Coord() const { return m_coord; }
    //! Coordinates accessor as reference
    tk::UnsMesh::Coords& Coord() { return m_coord; }
    //! Coordinates at time n accessor as const-ref
    const tk::UnsMesh::Coords& Coordn() const { return m_coordn; }

    //! Global ids accessors as const-ref
    const std::vector< std::size_t >& Gid() const { return m_gid; }

    //! Local ids accessors as const-ref
     const std::unordered_map< std::size_t, std::size_t >& Lid() const
    { return m_lid; }

    //! Tetrahedron element connectivity (with local ids) accessors as const-ref
    const std::vector< std::size_t >& Inpoel() const { return m_inpoel; }

    //! Mesh chunk accessor as const-ref
    const tk::UnsMesh::Chunk& Chunk() const { return m_el; }

    //! Mesh block id accessor as const-ref
    const std::unordered_map< std::size_t, std::set< std::size_t > >&
      ElemBlockId() const { return m_elemblockid; }

    //! Total mesh volume accessor
    tk::real meshvol() const { return m_meshvol; }

    //! Nodal mesh volume accessors const-ref
    const std::vector< tk::real >& V() const { return m_v; }

    //! Nodal mesh volumes at current time step accessors as const-ref
    const std::vector< tk::real >& Vol() const { return m_vol; }
    //! Nodal mesh volumes at previous time step accessors as const-ref
    const std::vector< tk::real >& Voln() const { return m_voln; }
    //! Nodal mesh volumes at previous time step accessors as ref
    std::vector< tk::real >& Voln() { return m_voln; }
    //! Element mesh volumes at t=t0 accessors as const-ref
    const std::vector< tk::real >& Vol0() const { return m_vol0; }
    //! Element mesh velocity accessor as const-ref
    const tk::Fields& MeshVel() const { return m_meshvel; }
    //! Element mesh velocity accessor as ref
    tk::Fields& MeshVel() { return m_meshvel; }

    //! Set 'initial' flag
    //! \param[in] i Value to put in 'initial'
    void Initial( std::size_t i ) { m_initial = i; }
    //! Query 'initial' flag
    //! \return True during setup, false durign time stepping
    bool Initial() const { return m_initial; }

    //! Update coordinates at time n
    void UpdateCoordn() { m_coordn = m_coord; }

    //! History points data accessor as const-ref
    const std::vector< HistData >& Hist() const { return m_histdata; }

    //! Box volume accessor
    tk::real& Boxvol() { return m_boxvol; }

    //! Block volume accessor
    std::vector< tk::real >& MeshBlkVol() { return m_meshblkvol; }

    //! Mesh ID accessor
    std::size_t MeshId() const { return m_meshid; }

    //! Time step size accessor
    tk::real Dt() const { return m_dt; }
    //! Time step size at previous time step accessor
    tk::real Dtn() const { return m_dtn; }
    //! Physical time accessor
    tk::real T() const { return m_t; }
    //! Iteration count accessor
    uint64_t It() const { return m_it; }

    //! Non-const-ref refinement iteration count accessor
    uint64_t& Itr() { return m_itr; }
    //! Non-const-ref field-output iteration count accessor
    uint64_t& Itf() { return m_itf; }

    //! Non-const-ref number of restarts accessor
    int& Nrestart() { return m_nrestart; }

    //! Timer accessor as const-ref
    const tk::Timer& Timer() const { return m_timer; }
    //! Timer accessor as non-const-ref
    tk::Timer& Timer() { return m_timer; }

    //! Accessor to flag indicating if the mesh was refined as a value
    int refined() const { return m_refined; }
    //! Accessor to flag indicating if the mesh was refined as non-const-ref
    int& refined() { return m_refined; }

    //! Transporter proxy accessor as const-ref
    const CProxy_Transporter& Tr() const { return m_transporter; }
    //! Transporter proxy accessor as non-const-ref
    CProxy_Transporter& Tr() { return m_transporter; }

    //! Access bound Refiner class pointer
    Refiner* Ref() const {
      Assert( m_refiner[ thisIndex ].ckLocal() != nullptr,
              "Refiner ckLocal() null" );
      return m_refiner[ thisIndex ].ckLocal();
    }

    //! Access Discretization proxy for a mesh
    CProxy_Discretization coupled( std::size_t meshid ) const {
      Assert( meshid < m_disc.size(),
              "No proxy for mesh ID " + std::to_string(meshid) );
      return m_disc[ meshid ];
    }

    //! Const-ref accessor to solver/mesh transfer configuration
    const std::vector< Transfer >& Transfers() const { return m_transfer; }

    //! Boundary node ids accessor as const-ref
    const std::unordered_map< std::size_t, std::size_t >& Bid() const
    { return m_bid; }

    //! Node communication map accessor as const-ref
    const tk::NodeCommMap& NodeCommMap() const { return m_nodeCommMap; }

    //! Edge communication map accessor as const-ref
    const tk::EdgeCommMap& EdgeCommMap() const { return m_edgeCommMap; }
    //@}

    //! Set time step size
    void setdt( tk::real newdt );

    //! Prepare for next step
    void next();

    //! Otput one-liner status report
    void status();

    //! Construct history output filename
    std::string histfilename( const std::string& id,
                              std::streamsize precision );

    //! Output headers for time history files (one for each point)
    void histheader( std::vector< std::string >&& names );

    //! Output time history for a time step
    void history( std::vector< std::vector< tk::real > >&& data );

    //! Output mesh and fields data (solution dump) to file(s)
    void write( const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const std::map< int, std::vector< std::size_t > >& bface,
                const std::map< int, std::vector< std::size_t > >& bnode,
                const std::vector< std::size_t >& triinpoel,
                const std::vector< std::string>& elemfieldnames,
                const std::vector< std::string>& nodefieldnames,
                const std::vector< std::string>& elemsurfnames,
                const std::vector< std::string>& nodesurfnames,
                const std::vector< std::vector< tk::real > >& elemfields,
                const std::vector< std::vector< tk::real > >& nodefields,
                const std::vector< std::vector< tk::real > >& elemsurfs,
                const std::vector< std::vector< tk::real > >& nodesurfs,
                CkCallback c );

    //! Zero grind-timer
    void grindZero();

    //! Detect if just returned from a checkpoint and if so, zero timers
    bool restarted( int nrestart );

    //! Remap mesh data due to new local ids
    void remap( const std::unordered_map< std::size_t, std::size_t >& map );

    //! \brief Function object for querying the node ids that belong to side
    //!   sets of the same type, called for each PDE type
    template< typename... tags >
    struct SidesetNodes {

      const std::map< int, std::vector< std::size_t > >& m_bface;
      const std::vector< std::size_t >& m_triinpoel;
      std::unordered_map< int, std::unordered_set< std::size_t > >& m_nodes;
      std::size_t m_mid;

      explicit
        SidesetNodes( const std::map< int, std::vector< std::size_t > >& bface,
                      const std::vector< std::size_t >& triinpoel,
                      std::unordered_map< int,
                        std::unordered_set< std::size_t > >& nodes,
                      std::size_t mid )
        : m_bface(bface), m_triinpoel(triinpoel), m_nodes(nodes), m_mid(mid)
      {
        const auto& bc =
          g_inputdeck.template get< tag::bc >();
        std::vector< std::size_t > ss;
        for (const auto& bci : bc) {
          const auto& bcm = bci.get< tag::mesh >();
          for (const auto& im : bcm) {
            // only if this bc is meant for current mesh
            // collect sidesets for this mesh with this bc type
            if (im-1 == m_mid) {
              ss.insert( ss.end(), bci.template get< tags... >().begin(),
                bci.template get< tags... >().end() );
            }
          }
        }
        for (const auto& s : ss) {
          auto k = m_bface.find(static_cast<int>(s));
          if (k != end(m_bface)) {
            auto& n = m_nodes[ k->first ];  // associate set id
            for (auto f : k->second) {      // face ids on side set
              n.insert( m_triinpoel[f*3+0] );
              n.insert( m_triinpoel[f*3+1] );
              n.insert( m_triinpoel[f*3+2] );
            }
          }
        }
      }
    };

    //! \brief Query nodes that belong to side sets of the same type for all
    //!   PDE types
    //! \tparam tags Tags addressing the location of a vector of vectors of
    //!   side set ids in the input deck
    //! \param[in] bface Boundary-faces mapped to side set ids
    //! \param[in] triinpoel Boundary-face connectivity
    //! \return Node ids that belong side sets of the same type (value),
    //!    associated to sides set id (key)
    template< typename... tags >
    std::unordered_map< int, std::unordered_set< std::size_t > >
    bcnodes( const std::map< int, std::vector< std::size_t > >& bface,
             const std::vector< std::size_t >& triinpoel ) const
    {
      std::unordered_map< int, std::unordered_set< std::size_t > > nodes;
      SidesetNodes< tags... >( bface, triinpoel, nodes, this->MeshId() );
      return nodes;
    }

    //! Find elements along our mesh chunk boundary
    std::vector< std::size_t > bndel() const;

    //! Decide if field output iteration count interval is hit
    bool fielditer() const;

    //! Decide if field output physics time interval is hit
    bool fieldtime() const;

    //! Decide if physics time falls into a field output time range
    bool fieldrange() const;

    //! Decide if history output iteration count interval is hit
    bool histiter() const;

    //! Decide if history output physics time interval is hit
    bool histtime() const;

    //! Decide if physics time falls into a history output time range
    bool histrange() const;

    //! Decide if this is the last time step
    bool finished() const;

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_meshid;
      p | m_transfer_complete;
      p | m_transfer;
      p | m_mytransfer;
      p | m_disc;
      p | m_nchare;
      p | m_it;
      p | m_itr;
      p | m_itf;
      p | m_initial;
      p | m_t;
      p | m_lastDumpTime;
      p | m_physFieldFloor;
      p | m_physHistFloor;
      p | m_rangeFieldFloor;
      p | m_rangeHistFloor;
      p | m_dt;
      p | m_dtn;
      p | m_nvol;
      p | m_nxfer;
      p | m_ale;
      p | m_transporter;
      p | m_meshwriter;
      p | m_refiner;
      p | m_el;
      if (p.isUnpacking()) {
        m_inpoel = std::get< 0 >( m_el );
        m_gid = std::get< 1 >( m_el );
        m_lid = std::get< 2 >( m_el );
      }
      p | m_coord;
      p | m_coordn;
      p | m_nodeCommMap;
      p | m_edgeCommMap;
      p | m_meshvol;
      p | m_v;
      p | m_vol;
      p | m_volc;
      p | m_voln;
      p | m_vol0;
      p | m_boxvol;
      p | m_meshblkvol;
      p | m_bid;
      p | m_timer;
      p | m_refined;
      p( reinterpret_cast<char*>(&m_prevstatus), sizeof(Clock::time_point) );
      p | m_nrestart;
      p | m_histdata;
      p | m_nsrc;
      p | m_ndst;
      p | m_meshvel;
      p | m_meshvel_converged;
      p | m_bface;
      p | m_triinpoel;
      p | m_elemblockid;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i Discretization object reference
    friend void operator|( PUP::er& p, Discretization& i ) { i.pup(p); }
    //@}

  private:
    // Shorthand for clock, setting an internal clock type
    using Clock = std::chrono::high_resolution_clock;

    //! Mesh ID
    std::size_t m_meshid;
    //! \brief Charm++ callback of the function to call after a mesh-to-mesh
    //!   solution transfer (to-and-fro) is complete
    CkCallback m_transfer_complete;
    //! Solution/mesh transfer (coupling) information coordination propagation
    //! \details This has the same size with the same src/dst information on
    //!   all solvers.
    std::vector< Transfer > m_transfer;
    //! My solution transfer/mesh (coupling) information
    //! \details This is a subset of m_transfer, holding only those entries
    //!   that this solver is involved in (either a source or a destination).
    std::vector< Transfer > m_mytransfer;
    //! Discretization proxies (one per mesh)
    std::vector< CProxy_Discretization > m_disc;
    //! Total number of Discretization chares
    int m_nchare;
    //! Iteration count
    uint64_t m_it;
    //! Iteration count with mesh refinement
    //! \details Used as the restart sequence number {RS} in saving output in
    //!    an ExodusII sequence
    //! \see https://www.paraview.org/Wiki/Restarted_Simulation_Readers
    uint64_t m_itr;
    //! Field output iteration count without mesh refinement
    //! \details Counts the number of field outputs to file during two
    //!   time steps with mesh efinement
    uint64_t m_itf;
    //! Flag that is nonzero during setup and zero during time stepping
    std::size_t m_initial;
    //! Physical time
    tk::real m_t;
    //! Physics time at last field output
    tk::real m_lastDumpTime;
    //! Recent floor of physics time divided by field output interval time
    tk::real m_physFieldFloor;
    //! Recent floor of physics time divided by history output interval time
    tk::real m_physHistFloor;
    //! Recent floors of physics time divided by field output time for ranges
    tk::real m_rangeFieldFloor;
    //! Recent floors of physics time divided by history output time for ranges
    tk::real m_rangeHistFloor;
    //! Physical time step size
    tk::real m_dt;
    //! Physical time step size at the previous time step
    tk::real m_dtn;
    //! \brief Number of chares from which we received nodal volume
    //!   contributions on chare boundaries
    std::size_t m_nvol;
    //! \brief Number of chares from which we received solution transfers
    //!   contributions on chare boundaries
    std::size_t m_nxfer;
    //! Distributed ALE proxy
    CProxy_ALE m_ale;
    //! Transporter proxy
    CProxy_Transporter m_transporter;
    //! Mesh writer proxy
    tk::CProxy_MeshWriter m_meshwriter;
    //! Mesh refiner proxy
    CProxy_Refiner m_refiner;
    //! \brief Elements of the mesh chunk we operate on
    //! \details Initialized by the constructor. The first vector is the element
    //!   connectivity (local IDs), the second vector is the global node IDs of
    //!   owned elements, while the third one is a map of global->local node
    //!   IDs.
    tk::UnsMesh::Chunk m_el;
    //! Alias to element connectivity
    std::vector< std::size_t >& m_inpoel = std::get<0>( m_el );
    //! Alias to global node IDs of owned elements
    std::vector< std::size_t >& m_gid = std::get<1>( m_el );
    //! \brief Alias to local node ids associated to the global ones of owned
    //!    elements
    std::unordered_map< std::size_t, std::size_t >& m_lid = std::get<2>( m_el );
    //! Mesh point coordinates
    tk::UnsMesh::Coords m_coord;
    //! Mesh coordinates at the time n for ALE
    tk::UnsMesh::Coords m_coordn;
    //! \brief Global mesh node IDs bordering the mesh chunk held by fellow
    //!   Discretization chares associated to their chare IDs
    tk::NodeCommMap m_nodeCommMap;
    //! \brief Edges with global node IDs bordering the mesh chunk held by
    //!   fellow Discretization chares associated to their chare IDs
    tk::EdgeCommMap m_edgeCommMap;
    //! Total mesh volume
    tk::real m_meshvol;
    //! Nodal mesh volumes
    //! \details This is the volume of the mesh associated to nodes of owned
    //!   elements (sum of surrounding cell volumes / 4) without contributions
    //!   from other chares on chare-boundaries
    std::vector< tk::real > m_v;
    //! Volume of nodes
    //! \details This is the volume of the mesh associated to nodes of owned
    //!   elements (sum of surrounding cell volumes / 4) with contributions from
    //!   other chares on chare-boundaries
    std::vector< tk::real > m_vol;
    //! Receive buffer for volume of nodes (with global node id as key)
    //! \details This is a communication buffer used to compute the volume of
    //!   the mesh associated to nodes of owned elements (sum of surrounding
    //!   cell volumes / 4) with contributions from other chares on
    //!   chare-boundaries.
    std::unordered_map< std::size_t, tk::real > m_volc;
    //! Volume of nodes at previous time step
    //! \details This is the volume of the mesh associated to nodes of owned
    //!   elements (sum of surrounding cell volumes / 4) with contributions from
    //!   other chares on chare-boundaries at the previous time step stage
    std::vector< tk::real > m_voln;
    //! Mesh element volumes at t=t0
    std::vector< tk::real > m_vol0;
    //! Volume of user-defined box IC
    tk::real m_boxvol;
    //! Volumes of mesh-blocks with user-defined ICs
    std::vector< tk::real > m_meshblkvol;
    //! \brief Local chare-boundary mesh node IDs at which we receive
    //!   contributions associated to global mesh node IDs of elements we
    //!   contribute to
    std::unordered_map< std::size_t, std::size_t > m_bid;
    //! Timer measuring a time step
    tk::Timer m_timer;
    //! 1 if mesh was refined in a time step, 0 if it was not
    int m_refined;
    //! Time point storing clock state at status()
    Clock::time_point m_prevstatus;
    //! Number of times restarted
    int m_nrestart;
    //! Data at history point locations
    std::vector< HistData > m_histdata;
    //! Number of transfers requested as a source
    std::size_t m_nsrc;
    //! Number of transfers requested as a destination
    std::size_t m_ndst;
    //! Mesh velocity if ALE is not enabled
    tk::Fields m_meshvel;
    //! \brief True if all stages of the time step converged the mesh velocity
    //!   linear solve in ALE
    bool m_meshvel_converged;
    //! Boundary faces side-set information
    std::map< int, std::vector< std::size_t > > m_bface;
    //! Triangle face connecitivity
    std::vector< std::size_t > m_triinpoel;
    //! Local tet ids associated with mesh block ids
    std::unordered_map< std::size_t, std::set< std::size_t > > m_elemblockid;

    //! Generate the Bid data-structure based on the node communication-map
    std::unordered_map< std::size_t, std::size_t > genBid();

    //! Generate {A,x,b} for Laplacian mesh velocity smoother
    std::tuple< tk::CSR, std::vector< tk::real >, std::vector< tk::real > >
    Laplacian( std::size_t ncomp ) const;

    //! Set mesh coordinates based on coordinates map
    tk::UnsMesh::Coords setCoord( const tk::UnsMesh::CoordMap& coordmap );

    //! Determine if communication of mesh transfer callbacks is complete
    bool transferCallbacksComplete() const;
};

} // inciter::

#endif // Discretization_h
