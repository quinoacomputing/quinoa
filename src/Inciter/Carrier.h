// *****************************************************************************
/*!
  \file      src/Inciter/Carrier.h
  \author    J. Bakosi
  \date      Fri 09 Sep 2016 02:35:37 PM MDT
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
      AEC [ label="AEC"
            tooltip="communication antidiffusive element contributions"
            URL="\ref inciter::Carrier::aec"];
      Allowed [ label="Allowed"
            tooltip="communication allowed solution increments and decrements"
            URL="\ref inciter::Carrier::allowed"];
      Limit [ label="Limit" color="#e6851c" style="filled"
            tooltip="perform limiting"
            URL="\ref inciter::Carrier::limit"];
      AEC -> Limit [ style="solid" ];
      Allowed -> Limit [ style="solid" ];
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
#include <set>

#include "Types.h"
#include "MeshNodes.h"
#include "Particles.h"
#include "DerivedData.h"
#include "VectorReducer.h"
#include "MeshNodeMerger.h"
#include "FluxCorrector.h"
#include "Inciter/InputDeck/InputDeck.h"

#include "NoWarning/transporter.decl.h"
#include "NoWarning/carrier.decl.h"
#include "NoWarning/particlewriter.decl.h"

namespace tk { class ExodusIIMeshWriter; }

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern CkReduction::reducerType VerifyBCMerger;
extern CkReduction::reducerType MeshNodeMerger;

//! Carrier Charm++ chare array used to advance transport equations in time
class Carrier : public CBase_Carrier {

  private:
    using TransporterProxy = CProxy_Transporter;
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Transporter,
                                                       CProxy_Carrier >;
    using ParticleWriterProxy = tk::CProxy_ParticleWriter< TransporterProxy >;

  public:
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
               const std::unordered_map< std::size_t, std::size_t >& cid,
               int ncarr );

    #if defined(__GNUC__)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Weffc++"
    #endif

    //! Migrate constructor
    explicit Carrier( CkMigrateMessage* ) {}

    #if defined(__GNUC__)
      #pragma GCC diagnostic pop
    #endif

    //! \brief Configure Charm++ reduction types
    //! \details Since this is a [nodeinit] routine, see carrier.ci, the
    //!   Charm++ runtime system executes the routine exactly once on every
    //!   logical node early on in the Charm++ init sequence. Must be static as
    //!   it is called without an object. See also: Section "Initializations at
    //!   Program Startup" at in the Charm++ manual
    //!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
    static void registerReducers() {
      VerifyBCMerger = CkReduction::addReducer( tk::mergeVector );
      MeshNodeMerger = CkReduction::addReducer( mergeMeshNodes< std::size_t > );
    }

    //! \brief Starts collecting global mesh node IDs bordering the mesh chunk
    //!   held by fellow Carrier chares associated to their chare IDs
    void comm();

    //! \brief Reduction target finishing collecting global mesh node IDs
    //!   bordering the mesh chunk held by fellow Carrier chares associated to
    //!   their chare IDs
    void msum( CkReductionMsg* msg );

    //! Initialize mesh IDs, element connectivity, coordinates
    void setup();

    //! Request owned node IDs on which a Dirichlet BC is set by the user
    void requestBCs();

    //! Look up and return old node ID for new one
    void oldID( int frompe, const std::vector< std::size_t >& newids );

    //! Look up boundary condition values at node IDs for all PDEs
    void bcval( int frompe, const std::vector< std::size_t >& nodes );

    //! Initialize communication and mesh data
    void init( tk::real dt );

    //! Update high order solution vector
    void updateHighSol( const std::vector< std::size_t >& gid,
                        const std::vector< tk::real >& sol );

    //! Update low order solution vector
    void updateLowSol( const std::vector< std::size_t >& gid,
                       const std::vector< tk::real >& sol );

    //! Advance equations to next stage in multi-stage time stepping
    void advance( uint8_t stage, tk::real dt, uint64_t it, tk::real t );

    //! Generates particles into mesh cells
    void genpar();

    //! Find particles missing by the requestor and make those found ours
    void findpar( int fromch,
                  const std::vector< std::size_t >& miss,
                  const tk::Particles& ps );

    //! Receive particle indices found elsewhere (by fellow neighbors)
    void foundpar( const std::vector< std::size_t >& found );

    //! Find particles missing by the requestor and make those found ours    
    void collectpar( int fromch,
                     const std::vector< std::size_t >& miss,
                     const tk::Particles& ps );

    //! Collect particle indices found elsewhere (by far fellows)
    void collectedpar( const std::vector< std::size_t >& found );

    //! Output mesh and particle fields to files
    void out();

    //! Output particles fields to file
    void doWriteParticles();

    //! Finish summing antidiffusive element contributions on chare-boundaries
    void finishaec( int fromch,
                    const std::vector< std::size_t >& gid,
                    const std::vector< std::vector< tk::real > >& P );

    //! \brief Acknowledge received antidiffusive element contributions on
    //!   chare-boundaries
    void recaec();

    //! Finish computing the allowed solution values on chare-boundaries
    void finishallowed( int fromch,
                        const std::vector< std::size_t >& gid,
                        const std::vector< std::vector< tk::real > >& A );

    //! \brief Acknowledge received allowed solution increments and decrements
    //!   on chare-boundaries
    void recallowed();

    //! Perform limiting as the final step of flux-corrected transport (FCT)
    void limit();

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_Carrier::pup(p);
      p | m_it;
      p | m_itf;
      p | m_t;
      p | m_stage;
      p | m_nhsol;
      p | m_nlsol;
      p | m_naec;
      p | m_nalw;
      p | m_nchpar;
      p | m_ncarr;
      p | m_outFilename;
      p | m_transporter;
      p | m_linsysmerger;
      p | m_particlewriter;
      p | m_fluxcorrector;
      p | m_cid;
      p | m_el;
      if (p.isUnpacking()) {
        m_inpoel = m_el.first;
        m_gid = m_el.second;
      }
      p | m_lid;
      p | m_coord;
      p | m_psup;
      p | m_uh;
      p | m_ul;
      p | m_uhf;
      p | m_ulf;
      p | m_uhn;
      p | m_uln;
      p | m_up;
      p | m_p;
      p | m_q;
      p | m_lhsd;
      p | m_lhso;
      p | m_particles;
      p | m_msum;
      p | m_parmiss;
      p | m_parelse;
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
    //! Stage in multi-stage time stepping
    uint8_t m_stage;
    //! Counter for high order solution nodes updated
    std::size_t m_nhsol;
    //! Counter for low order solution nodes updated
    std::size_t m_nlsol;
    //! \brief Number of chares received our anitdiffusive element contributions
    //!   on chare boundaries
    std::size_t m_naec;
    //! \brief Number of chares received our allowed solution increments and
    //!   decrements element contributions on chare boundaries
    std::size_t m_nalw;
    //! Number of chares we received particles from
    std::size_t m_nchpar;
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
    std::unordered_map< std::size_t, std::size_t > m_cid;
    //! \brief Elements of the mesh chunk we operate on
    //! \details Initialized by the constructor. The first vector is the element
    //!   connectivity (local IDs), while the second vector is the global node
    //!   IDs of owned elements.
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > m_el;
    //! Alias to element connectivity in m_el
    decltype(m_el.first)& m_inpoel = m_el.first;
    //! Alias to global node IDs of owned elements in m_el
    decltype(m_el.second)& m_gid = m_el.second;
    //!< Local node ids associated to the global ones of owned elements
    std::unordered_map< std::size_t, std::size_t > m_lid;
    //! Mesh point coordinates
    std::array< std::vector< tk::real >, 3 > m_coord;
    //! Flux corrector performing FCT
    FluxCorrector m_fluxcorrector;
    //! Points surrounding points of our chunk of the mesh
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > m_psup;
    //! Elements surrounding points of elements of mesh chunk we operate on
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
      m_esupel;
    //! Unknown/solution vectors: global mesh point row ids and values
    tk::MeshNodes m_uh, m_ul, m_uhf, m_ulf, m_uhn, m_uln, m_up, m_p, m_q;
    //! Sparse matrix sotring the diagonals and off-diagonals of nonzeros
    tk::MeshNodes m_lhsd, m_lhso;
    //! Particle properties
    tk::Particles m_particles;
    //! Element ID in which a particle has last been found for all particles
    std::vector< std::size_t > m_elp;
    //! \brief Global mesh node IDs bordering the mesh chunk held by fellow
    //!   Carrier chares associated to their chare IDs
    //! \details msum: mesh chunks surrounding mesh chunks and their neighbor
    //!   points
    std::unordered_map< int, std::vector< std::size_t > > m_msum;
    //! Indicies of particles not found here (missing)
    std::set< std::size_t > m_parmiss;
    //! Indicies of particles not found here but found by fellows
    decltype(m_parmiss) m_parelse;

    //! Send off global row IDs to linear system merger, setup global->local IDs
    void setupIds();

    //! Extract node IDs from element side sets and match to BCs
    std::vector< std::size_t > queryBCs();

    //! Query old node IDs for a list of new node IDs
    std::vector< std::size_t > old( const std::vector< std::size_t >& newids );

    //! Send node list to our LinSysMerger branch which is then used to set BCs
    void sendBCs( const std::vector< std::size_t >& bc );

    //! Read coordinates of mesh nodes given
    void readCoords();

    //! Compute left-hand side of transport equations
    void lhs();

    //! Compute righ-hand side vector of transport equations
    void rhs( tk::real mult,
              tk::real dt,
              const tk::MeshNodes& sol,
              tk::MeshNodes& rhs );

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

    //! Compute and sum antidiffusive element contributions to mesh nodes
    void aec( const tk::MeshNodes& Un, const tk::MeshNodes& Uh );

    //! Compute the allowed solution values at mesh nodes
    void allowed( const tk::MeshNodes& Ul, const tk::MeshNodes& Un );
};

} // inciter::

#endif // Carrier_h
