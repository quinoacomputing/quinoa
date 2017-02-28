// *****************************************************************************
/*!
  \file      src/Inciter/Carrier.h
  \author    J. Bakosi
  \date      Tue 28 Feb 2017 03:58:54 PM MST
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Carrier advances a system of transport equations
  \details   Carrier advances a system of transport equations. There are a
    potentially large number of Carrier Charm++ chares created by Transporter.
    Each carrier gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances a system of systems of PDEs in time.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/carrier.ci.

    #### Call graph ####
    The following is a directed acyclic graph (DAG) that outlines the
    asynchronous algorithm implemented in this class The detailed discussion of
    the algorithm is given in the Charm++ interface file transporter.ci, which
    also repeats the graph below using ASCII graphics. On the DAG orange
    fills denote global synchronization points that contain or eventually lead
    to global reductions. Dashed lines are potential shortcuts that allow
    jumping over some of the task-graph under some circumstances or optional
    code paths (taken, e.g., only in DEBUG mode). See the detailed discussion in
    carrier.ci.
    \dot
    digraph "Carrier SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      OwnAEC [ label="OwnAEC"
               tooltip="own contributions to the antidiffusive element
                        contributions computed"
               URL="\ref inciter::Carrier::aec"];
      ComAEC [ label="ComAEC"
               tooltip="contributions to the antidiffusive element contributions
                        communicated"
               URL="\ref inciter::Carrier::comaec"];
      OwnALW [ label="OwnALW"
               tooltip="own contributions to the maximum and minimum unknowns of
                        elements surrounding nodes computed"
               URL="\ref inciter::Carrier::alw"];
      ComALW [ label="ComALW"
               tooltip="contributions to the the maximum and minimum unknowns of
                        elements surrounding nodes communicated"
               URL="\ref inciter::Carrier::comalw"];
      Ver [ label="Ver" tooltip="verify antidiffusive element contributions"
            URL="\ref inciter::Carrier::verify"];
      OwnLim [ label="OwnLim"
               tooltip="compute limited antidiffusive element contributions"
               URL="\ref inciter::Carrier::lim"];
      ComLim [ label="ComLim"
               tooltip="contributions to the limited antidiffusive element
                        contributions communicated"
               URL="\ref inciter::Carrier::comlim"];
      Apply [ label="Apply"
              tooltip="apply limited antidiffusive element contributions"
              URL="\ref inciter::Carrier::limit"];
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
    \include Inciter/carrier.ci
*/
// *****************************************************************************
#ifndef Carrier_h
#define Carrier_h

#include <array>
#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>
#include <cstring>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <set>

#include "Types.h"
#include "Fields.h"
#include "Tracker.h"
#include "DerivedData.h"
#include "VectorReducer.h"
#include "FluxCorrector.h"
#include "Inciter/InputDeck/InputDeck.h"

#include "NoWarning/transporter.decl.h"
#include "NoWarning/carrier.decl.h"
#include "NoWarning/particlewriter.decl.h"

namespace tk { class ExodusIIMeshWriter; }

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern CkReduction::reducerType VerifyBCMerger;

//! Carrier Charm++ chare array used to advance transport equations in time
class Carrier : public CBase_Carrier {

  private:
    using TransporterProxy = CProxy_Transporter;
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Transporter,
                                                       CProxy_Carrier,
                                                       AuxSolverLumpMassDiff >;
    using ParticleWriterProxy = tk::CProxy_ParticleWriter< TransporterProxy >;

  public:
    using Array = CBase_Carrier;

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wunused-parameter"
    #elif defined(__GNUC__)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-parameter"
      #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(__INTEL_COMPILER)
      #pragma warning( push )
      #pragma warning( disable: 1478 )
    #endif
    // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
    // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
    Carrier_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(__GNUC__)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit
      Carrier( const TransporterProxy& transporter,
               const LinSysMergerProxy& lsm,
               const ParticleWriterProxy& pw,
               const std::vector< std::size_t >& conn,
               const std::unordered_map< int,
                       std::vector< std::size_t > >& msum,
               const std::unordered_map< std::size_t, std::size_t >& nodemap,
               const tk::UnsMesh::EdgeNodes& edgenodemap,
               int ncarr );

    //! Migrate constructor
    explicit Carrier( CkMigrateMessage* ) {}

    //! \brief Configure Charm++ reduction types
    //! \details Since this is a [nodeinit] routine, see carrier.ci, the
    //!   Charm++ runtime system executes the routine exactly once on every
    //!   logical node early on in the Charm++ init sequence. Must be static as
    //!   it is called without an object. See also: Section "Initializations at
    //!   Program Startup" at in the Charm++ manual
    //!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
    static void registerReducers() {
      VerifyBCMerger = CkReduction::addReducer( tk::mergeVector );
    }

    //! \brief Read mesh node coordinates, sum mesh volumes to nodes, and start
    //!   communicating them on chare-boundaries
    void vol();

    //! Collect nodal volumes across chare boundaries
    void comvol( const std::vector< std::size_t >& gid,
                 const std::vector< tk::real >& V );

    //! \brief Setup rows, query boundary conditions, generate particles, output
    //!    mesh, etc.
    void setup();

    //! Request owned node IDs on which a Dirichlet BC is set by the user
    void requestBCs();

    //! Look up and return old node ID for new one
    void oldID( int frompe, const std::vector< std::size_t >& newids );

    //! \brief Set ICs, compute initial time step size, output initial field
    //!   data, compute left-hand-side matrix
    void init();

    //! Update high order solution vector
    void updateSol( const std::vector< std::size_t >& gid,
                    const std::vector< tk::real >& sol );

    //! Update low order solution vector
    void updateLowSol( const std::vector< std::size_t >& gid,
                       const std::vector< tk::real >& sol );

    //! Advance equations to next stage in multi-stage time stepping
    void advance( uint8_t stage, tk::real newdt, uint64_t it, tk::real t );

    //! Compute time step size
    void dt();

    //! Find particles missing by the requestor and make those found ours
    void findpar( int fromch,
                  const std::vector< std::size_t >& miss,
                  const std::vector< std::vector< tk::real > >& ps )
    { m_tracker.findpar( thisProxy, m_coord, m_inpoel, fromch, miss, ps ); }

    //! Receive particle indices found elsewhere (by fellow neighbors)
    void foundpar( const std::vector< std::size_t >& found ) {
      m_tracker.foundpar( m_transporter, thisProxy, m_msum, this, thisIndex,
                          found );
    }

    //! Find particles missing by the requestor and make those found ours
    void collectpar( int fromch,
                     const std::vector< std::size_t >& miss,
                     const std::vector< std::vector< tk::real > >& ps )
    { m_tracker.collectpar( thisProxy, m_coord, m_inpoel, fromch, miss, ps ); }

    //! Collect particle indices found elsewhere (by far fellows)
    void collectedpar( const std::vector< std::size_t >& found )
    { m_tracker.collectedpar( m_transporter, this, found, m_ncarr ); }

    //! Extract the velocity at cell nodes
    std::array< std::array< tk::real, 4 >, 3 > velocity( std::size_t e );

    //! Output mesh and particle fields to files
    void out();

    //! Output particles fields to file
    void doWriteParticles();

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
      CBase_Carrier::pup(p);
      p | m_it;
      p | m_itf;
      p | m_t;
      p | m_dt;
      p | m_stage;
      p | m_nvol;
      p | m_nhsol;
      p | m_nlsol;
      p | m_naec;
      p | m_nalw;
      p | m_nlim;
      p | m_ncarr;
      p | m_outFilename;
      p | m_transporter;
      p | m_linsysmerger;
      p | m_particlewriter;
      p | m_fluxcorrector;
      p | m_nodemap;
      p | m_edgenodes;
      p | m_el;
      if (p.isUnpacking()) {
        m_inpoel = std::get< 0 >( m_el );
        m_gid = std::get< 1 >( m_el );
        m_lid = std::get< 2 >( m_el );
      }
      p | m_coord;
      p | m_psup;
      p | m_u;
      p | m_ul;
      p | m_uf;
      p | m_ulf;
      p | m_du;
      p | m_dul;
      p | m_p;
      p | m_q;
      p | m_a;
      p | m_lhsd;
      p | m_lhso;
      p | m_msum;
      p | m_vol;
      p | m_bid;
      p | m_pc;
      p | m_qc;
      p | m_ac;
      p | m_tracker;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i Carrier object reference
    //! \author J. Bakosi
    friend void operator|( PUP::er& p, Carrier& i ) { i.pup(p); }
    //@}

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

    //! Iteration count
    uint64_t m_it;
    //! Field output iteration count
    uint64_t m_itf;
    //! Physical time
    tk::real m_t;
    //! Physical time step size
    tk::real m_dt;
    //! Stage in multi-stage time stepping
    uint8_t m_stage;
    //! \brief Number of chares from which we received nodal volume
    //!   contributions on chare boundaries
    std::size_t m_nvol;
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
    //! Total number of carrier chares
    std::size_t m_ncarr;
    //! Output filename
    std::string m_outFilename;
    //! Transporter proxy
    TransporterProxy m_transporter;
    //! Linear system merger proxy
    LinSysMergerProxy m_linsysmerger;
    //! Particle writer proxy
    ParticleWriterProxy m_particlewriter;
    //! \brief Map associating old node IDs (as in file) to new node IDs (as in
    //!   producing contiguous-row-id linear system contributions)
    std::unordered_map< std::size_t, std::size_t > m_nodemap;
    //! \brief Maps associating node node IDs to edges (a pair of old node IDs)
    //!   for only the nodes newly added as a result of initial uniform
    //!   refinement.
    tk::UnsMesh::EdgeNodes m_edgenodes;
    //! \brief Elements of the mesh chunk we operate on
    //! \details Initialized by the constructor. The first vector is the element
    //!   connectivity (local IDs), the second vector is the global node IDs of
    //!   owned elements, while the third one is a map of global->local node
    //!   IDs.
    std::tuple< std::vector< std::size_t >,
                std::vector< std::size_t >,
                std::unordered_map< std::size_t, std::size_t > > m_el;
    //! Alias to element connectivity in m_el
    std::vector< std::size_t > m_inpoel = std::get< 0 >( m_el );
    //! Alias to global node IDs of owned elements in m_el
    std::vector< std::size_t > m_gid = std::get< 1 >( m_el );
    //! \brief Alias to local node ids associated to the global ones of owned
    //!    elements in m_el
    std::unordered_map< std::size_t, std::size_t > m_lid = std::get< 2 >( m_el );
    //! Mesh point coordinates
    std::array< std::vector< tk::real >, 3 > m_coord;
    //! Flux corrector performing FCT
    FluxCorrector m_fluxcorrector;
    //! Points surrounding points of our chunk of the mesh
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > m_psup;
    //! Unknown/solution vectors: global mesh point row ids and values
    tk::Fields m_u, m_ul, m_uf, m_ulf, m_du, m_dul, m_p, m_q, m_a;
    //! Sparse matrix sotring the diagonals and off-diagonals of nonzeros
    tk::Fields m_lhsd, m_lhso;
    //! \brief Global mesh node IDs bordering the mesh chunk held by fellow
    //!   Carrier chares associated to their chare IDs
    //! \details msum: mesh chunks surrounding mesh chunks and their neighbor
    //!   points
    std::unordered_map< int, std::vector< std::size_t > > m_msum;
    //! Volume of nodes of owned elements (sum of surrounding cell volumes / 4 )
    std::vector< tk::real > m_vol;
    //! \brief Local chare-boundary mesh node IDs at which we receive
    //!   contributions associated to global mesh node IDs of mesh elements we
    //!   contribute to
    std::unordered_map< std::size_t, std::size_t > m_bid;
    //! Receive buffers for FCT
    std::vector< std::vector< tk::real > > m_pc, m_qc, m_ac;
    //! Particle tracker
    tk::Tracker m_tracker;

    //! Query old node IDs for a list of new node IDs
    std::vector< std::size_t > old( const std::vector< std::size_t >& newids );

    //! \brief Extract node IDs from side set node lists and match to
    //    user-specified boundary conditions
    void bc();

    //! Read coordinates of mesh nodes given
    void readCoords();

    //! \brief Add coordinates of mesh nodes newly generated to edge-mid points
    //!    during initial refinement
    void addEdgeNodeCoords();

    //! Compute left-hand side of transport equations
    void lhs();

    //! Compute righ-hand side vector of transport equations
    void rhs( tk::real mult, const tk::Fields& sol );

    //! Output chare element blocks to output file
    void writeMesh();

    //! Output solution to file
    void writeSolution( const tk::ExodusIIMeshWriter& ew,
                        uint64_t it,
                        const std::vector< std::vector< tk::real > >& u ) const;

    //! Output mesh-based fields metadata to file
    void writeMeta() const;

    //! Output mesh-based fields to file
    void writeFields( tk::real time );

    //! Search particle ina single mesh cell
    bool parinel( std::size_t p, std::size_t e, std::array< tk::real, 4 >& N );

    //! Search particles in our chunk of the mesh
    void track();

    //! Advance particles based on velocity from mesh cell
    void advanceParticle( std::size_t i,
                          std::size_t e,
                          const std::array< tk::real, 4>& N );

    //! Apply boundary conditions to particles
    void applyParBC( std::size_t i );

    //! Try to find particles and add those found to the list of ours
    std::vector< std::size_t > addpar( const std::vector< std::size_t >& miss,
                                       const tk::Particles& ps );

    //! Output number of particles we will write to file in this step
    void writeParticles();

    //! Compute and sum antidiffusive element contributions (AEC) to mesh nodes
    void aec( const tk::Fields& Un );

    //! \brief Compute the maximum and minimum unknowns of all elements
    //!   surrounding nodes
    void alw( const tk::Fields& Un, const tk::Fields& Ul );

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

#endif // Carrier_h
