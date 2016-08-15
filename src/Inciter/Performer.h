// *****************************************************************************
/*!
  \file      src/Inciter/Performer.h
  \author    J. Bakosi
  \date      Mon 15 Aug 2016 10:23:40 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Performer advances a system of systems of PDEs
  \details   Performer advances a system of systems of PDEs. There are a
    potentially large number of Performer Charm++ chares created by Transporter.
    Each performer gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances a system of systems of PDEs in time.
*/
// *****************************************************************************
#ifndef Performer_h
#define Performer_h

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
#include "Inciter/InputDeck/InputDeck.h"

#include "NoWarning/transporter.decl.h"
#include "NoWarning/performer.decl.h"
#include "NoWarning/particlewriter.decl.h"

namespace tk { class ExodusIIMeshWriter; }

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern CkReduction::reducerType VerifyBCMerger;
extern CkReduction::reducerType MeshNodeMerger;

//! Performer Charm++ chare used to advance a PDE in time
class Performer : public CBase_Performer {

  private:
    using TransporterProxy = CProxy_Transporter;
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Transporter,
                                                       CProxy_Performer >;
    using ParticleWriterProxy = tk::CProxy_ParticleWriter< TransporterProxy >;

  public:
    //! Constructor
    explicit
      Performer( const CProxy_Transporter& transporter,
                 const LinSysMergerProxy& lsm,
                 const ParticleWriterProxy& pw,
                 const std::vector< std::size_t >& conn,
                 const std::unordered_map< std::size_t, std::size_t >& cid,
                 int nperf );

    #if defined(__GNUC__)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Weffc++"
    #endif

    //! Migrate constructor
    explicit Performer( CkMigrateMessage* ) {}

    #if defined(__GNUC__)
      #pragma GCC diagnostic pop
    #endif

    //! \brief Configure Charm++ reduction types
    //! \details Since this is a [nodeinit] routine, see performer.ci, the
    //!   Charm++ runtime system executes the routine exactly once on every
    //!   logical node early on in the Charm++ init sequence. Must be static as
    //!   it is called without an object. See also: Section "Initializations at
    //!   Program Startup" at in the Charm++ manual
    //!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
    static void registerReducers() {
      VerifyBCMerger = CkReduction::addReducer( tk::mergeVector );
      MeshNodeMerger = CkReduction::addReducer( mergeMeshNodes< std::size_t > );
    }

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

    //! Update solution vector
    void updateSolution( const std::vector< std::size_t >& gid,
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

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_Performer::pup(p);
      p | m_it;
      p | m_itf;
      p | m_t;
      p | m_stage;
      p | m_nsol;
      p | m_nchpar;
      p | m_nperf;
      p | m_outFilename;
      p | m_transporter;
      p | m_linsysmerger;
      p | m_particlewriter;
      p | m_cid;
      p | m_el;
      if (p.isUnpacking()) {
        m_inpoel = m_el.first;
        m_gid = m_el.second;
      }
      p | m_lid;
      p | m_coord;
      p | m_psup;
      p | m_u;
      p | m_uf;
      p | m_un;
      p | m_up;
      p | m_lhsd;
      p | m_lhso;
      p | m_particles;
      p | m_msum;
      p | m_parmiss;
      p | m_parelse;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i Performer object reference
    //! \author J. Bakosi
    friend void operator|( PUP::er& p, Performer& i ) { i.pup(p); }
    //@}

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

    uint64_t m_it;                       //!< Iteration count
    uint64_t m_itf;                      //!< Field output iteration count
    tk::real m_t;                        //!< Physical time
    uint8_t m_stage;                     //!< Stage in multi-stage time stepping
    std::size_t m_nsol;                  //!< Counter for solution nodes updated
    std::size_t m_nchpar;                //!< Numbr of chares recvd partcls from
    std::size_t m_nperf;                 //!< Total number of performer chares
    std::string m_outFilename;           //!< Output filename
    TransporterProxy m_transporter;      //!< Transporter proxy
    LinSysMergerProxy m_linsysmerger;    //!< Linear system merger proxy
    ParticleWriterProxy m_particlewriter;//!< Particle writer proxy
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
    //! Points surrounding points of our chunk of the mesh
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > m_psup;
    //! Elements surrounding points of elements of mesh chunk we operate on
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
      m_esupel;
    //! Unknown/solution vector: global mesh point row ids and values
    tk::MeshNodes m_u, m_uf, m_un, m_up;
    //! Sparse matrix sotring the diagonals and off-diagonals of nonzeros
    tk::MeshNodes m_lhsd, m_lhso;
    //! Particle properties
    tk::Particles m_particles;
    //! Element ID in which a particle has last been found for all particles
    std::vector< std::size_t > m_elp;
    //! Fellow Performer chare indices holding neighboring mesh chunks
    std::vector< int > m_msum;
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

    //! Compute left-hand side matrix of PDE
    void lhs();

    //! Compute righ-hand side vector of PDE
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
};

} // inciter::

#endif // Performer_h
