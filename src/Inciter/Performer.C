// *****************************************************************************
/*!
  \file      src/Inciter/Performer.C
  \author    J. Bakosi
  \date      Wed 03 Aug 2016 12:49:08 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Performer advances a PDE
  \details   Performer advances a PDE. There are a potentially
    large number of Performer Charm++ chares created by Conductor. Each
    performer gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances a PDE in time.
*/
// *****************************************************************************

#include <string>
#include <cmath>
#include <array>

#include <gm19.h>

#include "Performer.h"
#include "Tracker.h"
#include "LinSysMerger.h"
#include "Vector.h"
#include "Reader.h"
#include "ContainerUtil.h"
#include "UnsMesh.h"
#include "Reorder.h"
#include "ExodusIIMeshReader.h"
#include "ExodusIIMeshWriter.h"
#include "ParticleWriter.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "DerivedData.h"
#include "PDE.h"
#include "RNGSSE.h"

// Force the compiler to not instantiate the template below as it is
// instantiated in LinSys/LinSysMerger.C (only required on mac)
extern template class tk::LinSysMerger< inciter::CProxy_Conductor,
                                        inciter::CProxy_Performer >;

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< PDE > g_pdes;

//! \brief Charm++ reducers used by Performer
//! \details These variables are defined here in the .C file and declared as
//!   extern in Performer.h. If instead one defines it in the header (as static),
//!   a new version of the variable is created any time the header file is
//!   included, yielding no compilation nor linking errors. However, that leads
//!   to runtime errors, since Performer::registerReducers(), a Charm++
//!   "initnode" entry method, *may* fill one while contribute() may use the
//!   other (unregistered) one. Result: undefined behavior, segfault, and
//!   formatting the internet ...
CkReduction::reducerType VerifyBCMerger;
CkReduction::reducerType MeshNodeMerger;

} // inciter::

using inciter::Performer;

Performer::Performer(
  const ConductorProxy& conductor,
  const LinSysMergerProxy& lsm,
  const TrackerProxy& tracker,
  const ParticleWriterProxy& pw,
  const std::vector< std::size_t >& conn,
  const std::unordered_map< std::size_t, std::size_t >& cid,
  int nperf )
:
  m_it( 0 ),
  m_itf( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_stage( 0 ),
  m_nsol( 0 ),
  m_nperf( nperf ),
  m_outFilename( g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "." +
                 std::to_string( thisIndex ) ),
  m_conductor( conductor ),
  m_linsysmerger( lsm ),
  m_tracker( tracker ),
  m_particlewriter( pw ),
  m_cid( cid ),
  m_el( tk::global2local( conn ) ),     // fills m_inpoel and m_gid
  m_lid(),
  m_coord(),
  m_psup( tk::genPsup( m_inpoel, 4, tk::genEsup(m_inpoel,4) ) ),
  m_u( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_uf( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_un( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_up( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_lhsd( m_psup.second.size()-1, g_inputdeck.get< tag::component >().nprop() ),
  m_lhso( m_psup.first.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_particles( g_inputdeck.get< tag::param, tag::compns, tag::npar >() *
               m_inpoel.size()/4, 3 ),
  m_sum()
// *****************************************************************************
//  Constructor
//! \param[in] conductor Host (Conductor) proxy
//! \param[in] lsm Linear system merger (LinSysMerger) proxy
//! \param[in] tracker Passive tracker proxy
//! \param[in] conn Vector of mesh element connectivity owned (global IDs)
//! \param[in] cid Map associating old node IDs (as in file) to new node IDs (as
//!   in producing contiguous-row-id linear system contributions)
//! \param[in] nper Total number of Performer chares
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( m_psup.second.size()-1 == m_gid.size(),
          "Number of mesh points and number of global IDs unequal" );

  // Register ourselves with the linear system merger
  m_linsysmerger.ckLocalBranch()->checkin();
  // Send number of partciles we will contribute to particle writer
  if (g_inputdeck.get< tag::param, tag::compns, tag::npar >() > 0) {
    m_particlewriter.ckLocalBranch()->npar( m_particles.nunk() );
    // Signal back to host that we are done with sending our number of particles
    contribute(
        CkCallback( CkReductionTarget(Conductor,nparcomplete), m_conductor ) );
    // Collect fellow chare IDs our mesh neighbors
    std::vector< std::pair< int, std::unordered_set< std::size_t > > >
      meshnodes{ { thisIndex, { begin(m_gid), end(m_gid) } } };
    auto stream = serialize( meshnodes );
    CkCallback cb( CkIndex_Performer::msum(nullptr), thisProxy );
    contribute( stream.first, stream.second.get(), MeshNodeMerger, cb );
  }
}

void
Performer::msum( CkReductionMsg* msg )
// *****************************************************************************
// Reduction target finishing collecting fellow chare mesh node IDs
//! \param[in] msg Serialized aggregated mesh node IDs categorized by chares
//! \author J. Bakosi
// *****************************************************************************
{
  std::vector< std::pair< int, std::unordered_set< std::size_t > > > meshnodes;
  PUP::fromMem creator( msg->getData() );
  creator | meshnodes;
  delete msg;

  for (const auto& c : meshnodes)
    for (auto i : m_gid)
      if (c.first != thisIndex && c.second.find(i) != c.second.end()) {
        m_sum.push_back( c.first );
        break;
      }

//   std::cout << thisIndex << ": ";
//   for (auto i : m_sum) std::cout << i << ' ';
//   std::cout << '\n';

  contribute(
     CkCallback( CkReductionTarget(Conductor,msumcomplete), m_conductor ) );
}

void
Performer::setup()
// *****************************************************************************
// Initialize mesh IDs, element connectivity, coordinates
//! \author J. Bakosi
// *****************************************************************************
{
  // Send off global row IDs to linear system merger, setup global->local IDs
  setupIds();
  // Send node IDs from element side sets matched to user-specified BCs
  sendBCs( queryBCs() );
  // Read coordinates of owned and received mesh nodes
  readCoords();
  // Generate particles
  genPar();
  // Output chare mesh to file
  writeMesh();
  // Output fields metadata to output file
  writeMeta();
}

void
Performer::setupIds()
// *****************************************************************************
// Send off global row IDs to linear system merger, setup global->local IDs
//! \author J. Bakosi
// *****************************************************************************
{
  // Send off global row IDs to linear system merger
  m_linsysmerger.ckLocalBranch()->charerow( thisIndex, m_gid );
  // Associate local node IDs to global ones
  for (std::size_t i=0; i<m_gid.size(); ++i) m_lid[ m_gid[i] ] = i;
}

std::vector< std::size_t >
Performer::queryBCs()
// *****************************************************************************
//  Extract nodes IDs from side sets node lists and match to boundary conditions
//! \return List of owned node IDs on which a Dirichlet BC is set by the user
//! \details Boundary conditions, mathematically speaking, are applied on
//!   finite surfaces. These finite surfaces are given by element sets. This
//!   function takes the node lists mapped to side set IDs and extracts only
//!   those nodes that we own that the user has specified boundary condition on
//!   and returns the 'new' (as in producing contiguous-row-id linear system
//!   contributions, see also Partitioner.h) node IDs.
//! \author J. Bakosi
// *****************************************************************************
{
  // Access all side sets from LinSysMerger
  auto& side = m_linsysmerger.ckLocalBranch()->side();

  // lambda to find out if we own the old (as in file) global node id. We
  // perform a linear search on the map value so we do not have to invert the
  // map and return the map key (the new global id, new as in producing
  // contiguous-row-id linear system contributions, see also Partitioner.h).
  auto own = [ this ]( std::size_t id ) -> std::pair< bool, std::size_t > {
    for (const auto& i : this->m_cid)
      if (i.second == id)
        return { true, i.first };
    return { false, 0 };
  };

  // lambda to find out if the user has specified a Dirichlet BC on the given
  // side set for any component of any of the PDEs integrated
  auto userbc = []( int sideset ) -> bool {
    for (const auto& eq : g_pdes) if (eq.anydirbc(sideset)) return true;
    return false;
  };

  // Collect list of owned node IDs on which a Dirichlet BC is given by the user
  std::vector< std::size_t > bc;
  for (const auto& s : side) {
    auto ubc = userbc( s.first );
    for (auto n : s.second) {
      auto o = own(n);
      if (o.first && ubc) bc.push_back( o.second );
    }
  }

  // Make BC node list unique. The node list extracted here can still contain
  // repeating IDs, because the nodes are originally extracted by the Exodus
  // reader by interrogating elements whose faces are adjacent to side sets,
  // which can contain repeating node IDs for those nodes that elements share.
  // See also the Exodus user manual.
  tk::unique(bc);

  return bc;
}

void
Performer::sendBCs( const std::vector< std::size_t >& bc )
// *****************************************************************************
// Send node list to our LinSysMerger branch which is then used to set BCs
//! \param[in] bc List of node IDs to send
//! \author J. Bakosi
// *****************************************************************************
{
  // Offer list of owned node IDs mapped to side sets to LinSysMerger
  m_linsysmerger.ckLocalBranch()->charebc( thisIndex, bc );
}

std::vector< std::size_t >
Performer::old( const std::vector< std::size_t >& newids )
// *****************************************************************************
// Query old node IDs for a list of new node IDs
//! \param[in] newids Vector of new node IDs
//! \details 'old' as in file, 'new' as in as in producing contiguous-row-id
//!   linear system contributions, see also Partitioner.h.
//! \author J. Bakosi
// *****************************************************************************
{
  std::vector< std::size_t > oldids;
  for (auto i : newids) oldids.push_back( tk::cref_find(m_cid,i) );
  return oldids;
}

void
Performer::requestBCs()
// *****************************************************************************
// Request owned node IDs on which a Dirichlet BC is set by the user
//! \details Called from host (Conductor), contributing the result back to host
//! \author J. Bakosi
// *****************************************************************************
{
   auto stream = tk::serialize( old( queryBCs() ) );
   CkCallback cb( CkIndex_Conductor::doverifybc(nullptr), m_conductor );
   contribute( stream.first, stream.second.get(), VerifyBCMerger, cb );
}

void
Performer::oldID( int frompe, const std::vector< std::size_t >& newids )
// *****************************************************************************
// Look up and return old node IDs for new ones
//! \param[in] newids Vector of new node IDs
//! \details Old (as in file), new (as in producing contiguous-row-id linear
//!   system contributions
//! \author J. Bakosi
// *****************************************************************************
{
  m_linsysmerger[ frompe ].oldID( thisIndex, old(newids) );
}

void
Performer::bcval( int frompe, const std::vector< std::size_t >& nodes )
// *****************************************************************************
// Look up boundary condition values at node IDs for all PDEs
//! \param[in] nodes Vector of node IDs at which to query BC values
//! \author J. Bakosi
// *****************************************************************************
{
  // Access all side sets from LinSysMerger
  auto& side = m_linsysmerger.ckLocalBranch()->side();

  // lambda to query the user-specified Dirichlet BCs on a given side set for
  // all components of all the PDEs integrated
  auto bc = []( int sideset ) -> std::vector< std::pair< bool, tk::real > > {
    std::vector< std::pair< bool, tk::real > > b;
    for (const auto& eq : g_pdes) {
      auto e = eq.dirbc( sideset );  // query BC values for all components of eq
      b.insert( end(b), begin(e), end(e) );
    }
    return b;
  };

  // lambda to find out whether 'new' node id is in the list of 'old' node ids
  // given in s (which contains the old node ids of a side set). Here 'old'
  // means as in file, while 'new' means as in producing contiguous-row-id
  // linear system contributions. See also Partitioner.h.
  auto inset = [ this ]( std::size_t id, const std::vector< std::size_t >& s )
  -> bool {
    for (auto n : s) if (tk::cref_find(this->m_cid,id) == n) return true;
    return false;
  };

  // Collect vector of pairs of bool and BC value, where the bool indicates
  // whether the BC value is set at the given node by the user. The size of the
  // vectors is the number of PDEs integrated times the number of scalar
  // components in all PDEs. The vector is associated to global node IDs at
  // which the boundary condition will be set. If a node belongs to multiple
  // side sets in the file, we keep it associated to the first side set given by
  // the user.
  std::unordered_map< std::size_t,
                      std::vector< std::pair< bool, tk::real > > > bcv;
  for (const auto& s : side) {
    auto b = bc( s.first );    // query BC values for all components of all PDEs
    for (auto n : nodes)
      if (inset( n, s.second )) {
        auto& v = bcv[n];
        if (v.empty()) v.insert( begin(v), begin(b), end(b) );
      }
  }

//  for (const auto& n : bcv) {
//     std::cout << n.first+1 << ':';
//     for (const auto& p : n.second) {
//       if (p.first) std::cout << p.second << ','; else std::cout << "*,";
//     }
//     std::cout << '\n';
//   }
//   std::cout << "----\n";

  m_linsysmerger[ frompe ].charebcval( bcv );
}

void
Performer::init( tk::real dt )
// *****************************************************************************
// Initialize linear system
//! \author J. Bakosi
// *****************************************************************************
{
  // Set initial conditions for all PDEs
  for (const auto& eq : g_pdes) eq.initialize( m_coord, m_u, m_t );

  // Output initial conditions to file (it = 0, time = 0.0)
  writeFields( m_it, m_t );

  // Send off initial conditions for assembly
  m_linsysmerger.ckLocalBranch()->charesol( thisIndex, m_gid, m_u );

  // Call back to Conductor::initcomplete(), signaling that the initialization
  // is complete and we are now starting time stepping
  contribute(
      CkCallback( CkReductionTarget( Conductor, initcomplete ), m_conductor ) );

  // Compute left-hand side of PDE
  lhs();
  // Start advancing PDE in time at time step stage 0
  advance( 0, dt, m_it, m_t );
}

void
Performer::lhs()
// *****************************************************************************
// Compute left-hand side of PDE
//! \author J. Bakosi
// *****************************************************************************
{
  // Compute left-hand side matrix for all equations
  for (const auto& eq : g_pdes)
    eq.lhs( m_coord, m_inpoel, m_psup, m_lhsd, m_lhso );

  // Send off left hand side for assembly
  m_linsysmerger.ckLocalBranch()->
    charelhs( thisIndex, m_gid, m_psup, m_lhsd, m_lhso );
}

void
Performer::rhs( tk::real mult,
                tk::real dt,
                const tk::MeshNodes& sol,
                tk::MeshNodes& rhs )
// *****************************************************************************
// Compute right-hand side of PDE
//! \param[in] mult Multiplier differentiating the different stages in
//!    multi-stage time stepping
//! \param[in] dt Size of time step
//! \param[in] sol Solution vector at current stage
//! \param[inout] rhs Right-hand side vector computed
//! \author J. Bakosi
// *****************************************************************************
{
  // Compute right-hand side vector for all equations
  for (const auto& eq : g_pdes)
    eq.rhs( mult, dt, m_coord, m_inpoel, sol, m_u, rhs );

  // Send off right-hand sides for assembly
  m_linsysmerger.ckLocalBranch()->charerhs( thisIndex, m_gid, rhs );
}

void
Performer::readCoords()
// *****************************************************************************
//  Read coordinates of mesh nodes from file
//! \author J. Bakosi
// *****************************************************************************
{
  tk::ExodusIIMeshReader
    er( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );

  auto& x = m_coord[0];
  auto& y = m_coord[1];
  auto& z = m_coord[2];
  for (auto p : m_gid) er.readNode( tk::cref_find(m_cid,p), x, y, z );
}

void
Performer::writeMesh()
// *****************************************************************************
// Output chare element blocks to file
//! \author J. Bakosi
// *****************************************************************************
{
  // Create ExodusII writer
  tk::ExodusIIMeshWriter ew( m_outFilename, tk::ExoWriter::CREATE );
  // Write chare mesh initializing element connectivity and point coords
  ew.writeMesh( tk::UnsMesh( m_inpoel, m_coord ) );
}

void
Performer::writeSolution( const tk::ExodusIIMeshWriter& ew,
                          uint64_t it,
                          const std::vector< std::vector< tk::real > >& u )
  const
// *****************************************************************************
// Output solution to file
//! \param[in] ew ExodusII mesh-based writer object
//! \param[in] it Iteration count
//! \param[in] varid Exodus variable ID
//! \param[in] u Vector of fields to write to file
//! \author J. Bakosi
// *****************************************************************************
{
  int varid = 0;
  for (const auto& f : u) ew.writeNodeScalar( it, ++varid, f );
}

void
Performer::writeMeta() const
// *****************************************************************************
// Output mesh-based fields metadata to file
//! \author J. Bakosi
// *****************************************************************************
{
  // Create ExodusII writer
  tk::ExodusIIMeshWriter ew( m_outFilename, tk::ExoWriter::OPEN );

  // Collect nodal field output names from all PDEs
  std::vector< std::string > names;
  for (const auto& eq : g_pdes) {
    auto n = eq.names();
    names.insert( end(names), begin(n), end(n) );
  }

  // Write node field names
  ew.writeNodeVarNames( names );
}

void
Performer::writeFields( uint64_t it, tk::real time )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] it Iteration count
//! \param[in] time Physical time
//! \author J. Bakosi
// *****************************************************************************
{
  // Increase field output iteration count
  ++m_itf;

  // Create ExodusII writer
  tk::ExodusIIMeshWriter ew( m_outFilename, tk::ExoWriter::OPEN );

  // Write time stamp
  ew.writeTimeStamp( m_itf, time );

  // Collect node fields output from all PDEs
  m_un = m_u;   // make a copy as eq::output() is allowed to overwrite its arg
  std::vector< std::vector< tk::real > > output;
  for (const auto& eq : g_pdes) {
    auto o = eq.output( time, m_coord, m_un );
    output.insert( end(output), begin(o), end(o) );
  }
  // Write node fields
  writeSolution( ew, m_itf, output );

  // Write particle fields
  if (g_inputdeck.get< tag::param, tag::compns, tag::npar >() > 0) {
    m_particlewriter.ckLocalBranch()->writeTimeStamp( it );
    m_particlewriter.ckLocalBranch()->writeCoords( m_particles.extract(0,0),
                                                   m_particles.extract(1,0),
                                                   m_particles.extract(2,0) );
  }
}

void
Performer::advance( uint8_t stage, tk::real dt, uint64_t it, tk::real t )
// *****************************************************************************
// Advance equations to next stage in multi-stage time stepping
//! \param[in] stage Stage in multi-stage time stepping
//! \param[in] dt Size of time step
//! \param[in] it Iteration count
//! \param[in] t Physical time
//! \author J. Bakosi
// *****************************************************************************
{
  // Update local copy of time step stage
  m_stage = stage;

  // Advance stage in multi-stage time stepping by updating the rhs
  if (m_stage < 1) {

    rhs( 0.5, dt, m_u, m_uf );

  } else {

    // Update local copy of physical time and iteration count at the final stage
    m_t = t;
    m_it = it;

    rhs( 1.0, dt, m_uf, m_un );
    
  }
}

void
Performer::genPar()
// *****************************************************************************
// Generate particles to each of our mesh cells
//! \author F.J. Gonzalez
// *****************************************************************************
{
  auto rng = tk::RNGSSE< gm19_state, unsigned, gm19_generate_ >
                       ( static_cast<unsigned>(m_nperf), gm19_init_sequence_ );

  // Create a reference of mesh point coordinates
  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  // Generate npar number of particles into each mesh cell
  auto npar = g_inputdeck.get< tag::param, tag::compns, tag::npar >();
  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    for (std::size_t p=0; p<npar; ++p) {
      std::array< tk::real, 4 > N;
      rng.uniform( thisIndex, 3, N.data() );
      N[3] = 1.0 - N[0] - N[1] - N[2];
      if ( std::min(N[0],1-N[0]) > 0 && std::min(N[1],1-N[1]) > 0 &&
           std::min(N[2],1-N[2]) > 0 && std::min(N[3],1-N[3]) > 0 ) {
        const auto A = m_inpoel[e*4+0];
        const auto B = m_inpoel[e*4+1];
        const auto C = m_inpoel[e*4+2];
        const auto D = m_inpoel[e*4+3];
        const auto i = e * npar + p;
        m_particles( i, 0, 0 ) = x[A]*N[0] + x[B]*N[1] + x[C]*N[2] + x[D]*N[3];
        m_particles( i, 1, 0 ) = y[A]*N[0] + y[B]*N[1] + y[C]*N[2] + y[D]*N[3];
        m_particles( i, 2, 0 ) = z[A]*N[0] + z[B]*N[1] + z[C]*N[2] + z[D]*N[3];
        // std::cout"p "<<i<<"in e "<<e<<std::endl;
      } else --p; // retry if particle was not generated into cell
    }
  }
}

void
Performer::track()
// *****************************************************************************
// Search particles in our chunk of the mesh
//! \author F.J. Gonzalez
// *****************************************************************************
{
  // Create a reference of mesh point coordinates
  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  bool found = false;
  // Loop over the number of particles and evaluate shapefunctions at each
  // particle's location
  for (std::size_t i=0; i<m_particles.nunk(); ++i) {
    // Create a reference of particle coordinates
    const auto& xp = m_particles( i, 0, 0 );
    const auto& yp = m_particles( i, 1, 0 );
    const auto& zp = m_particles( i, 2, 0 );
    for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
      const auto A = m_inpoel[e*4+0];
      const auto B = m_inpoel[e*4+1];
      const auto C = m_inpoel[e*4+2]; 
      const auto D = m_inpoel[e*4+3];

      // Evaluate shapefunctions at particle locations using Cramer's Rule
      //
      //    | xp |   | x1 x2 x3 x4 |   | N1 |
      //    | yp | = | y1 y2 y3 y4 | • | N2 |
      //    | zp |   | z1 z2 z3 z4 |   | N3 |
      //    | 1  |   | 1  1  1  1  |   | N4 |
      //
      
      tk::real DetX = (y[B]*z[C] - y[C]*z[B] - y[B]*z[D] + y[D]*z[B] + 
        y[C]*z[D] - y[D]*z[C])*x[A] + x[B]*y[C]*z[A] - x[B]*y[A]*z[C] +
        x[C]*y[A]*z[B] - x[C]*y[B]*z[A] + x[B]*y[A]*z[D] - x[B]*y[D]*z[A] -
        x[D]*y[A]*z[B] + x[D]*y[B]*z[A] - x[C]*y[A]*z[D] + x[C]*y[D]*z[A] +
        x[D]*y[A]*z[C] - x[D]*y[C]*z[A] - x[B]*y[C]*z[D] + x[B]*y[D]*z[C] +
        x[C]*y[B]*z[D] - x[C]*y[D]*z[B] - x[D]*y[B]*z[C] + x[D]*y[C]*z[B];

      tk::real DetX1 = (y[D]*z[C] - y[C]*z[D] + y[C]*zp - yp*z[C] -
        y[D]*zp + yp*z[D])*x[B] + x[C]*y[B]*z[D] - x[C]*y[D]*z[B] - 
        x[D]*y[B]*z[C] + x[D]*y[C]*z[B] - x[C]*y[B]*zp + x[C]*yp*z[B] +
        xp*y[B]*z[C] - xp*y[C]*z[B] + x[D]*y[B]*zp - x[D]*yp*z[B] - 
        xp*y[B]*z[D] + xp*y[D]*z[B] + x[C]*y[D]*zp - x[C]*yp*z[D] - 
        x[D]*y[C]*zp + x[D]*yp*z[C] + xp*y[C]*z[D] - xp*y[D]*z[C];

      tk::real DetX2 = (y[C]*z[D] - y[D]*z[C] - y[C]*zp + yp*z[C] +
        y[D]*zp - yp*z[D])*x[A] + x[C]*y[D]*z[A] - x[C]*y[A]*z[D] +
        x[D]*y[A]*z[C] - x[D]*y[C]*z[A] + x[C]*y[A]*zp - x[C]*yp*z[A] -
        xp*y[A]*z[C] + xp*y[C]*z[A] - x[D]*y[A]*zp + x[D]*yp*z[A] +
        xp*y[A]*z[D] - xp*y[D]*z[A] - x[C]*y[D]*zp + x[C]*yp*z[D] + 
        x[D]*y[C]*zp - x[D]*yp*z[C] - xp*y[C]*z[D] + xp*y[D]*z[C];

      tk::real DetX3 = (y[D]*z[B] - y[B]*z[D] + y[B]*zp - yp*z[B] -
        y[D]*zp + yp*z[D])*x[A] + x[B]*y[A]*z[D] - x[B]*y[D]*z[A] -
        x[D]*y[A]*z[B] + x[D]*y[B]*z[A] - x[B]*y[A]*zp + x[B]*yp*z[A] +
        xp*y[A]*z[B] - xp*y[B]*z[A] + x[D]*y[A]*zp - x[D]*yp*z[A] - 
        xp*y[A]*z[D] + xp*y[D]*z[A] + x[B]*y[D]*zp - x[B]*yp*z[D] - 
        x[D]*y[B]*zp + x[D]*yp*z[B] + xp*y[B]*z[D] - xp*y[D]*z[B];

      tk::real DetX4 = (y[B]*z[C] - y[C]*z[B] - y[B]*zp + yp*z[B] +
        y[C]*zp - yp*z[C])*x[A] + x[B]*y[C]*z[A] - x[B]*y[A]*z[C] +
        x[C]*y[A]*z[B] - x[C]*y[B]*z[A] + x[B]*y[A]*zp - x[B]*yp*z[A] -
        xp*y[A]*z[B] + xp*y[B]*z[A] - x[C]*y[A]*zp + x[C]*yp*z[A] +
        xp*y[A]*z[C] - xp*y[C]*z[A] - x[B]*y[C]*zp + x[B]*yp*z[C] +
        x[C]*y[B]*zp - x[C]*yp*z[B] - xp*y[B]*z[C] + xp*y[C]*z[B];

      std::array< tk::real, 4 > N {{ DetX1/DetX,
                                     DetX2/DetX,
                                     DetX3/DetX,
                                     DetX4/DetX }};

      // If particle is found, advance, and process next one
      if ( std::min(N[0],1-N[0]) > 0 && std::min(N[1],1-N[1]) > 0 &&
           std::min(N[2],1-N[2]) > 0 && std::min(N[3],1-N[3]) > 0 ) {
        advanceParticle( i, e, N );
        //std::cout<<"a: p "<<i<<"in e "<<e<<std::endl;
        found = true;
        e = m_inpoel.size()/4;  // search for next particle
      }
    }
    if ( found == false ) {
      std::cout<<"particle "<<i<<" not found"<<std::endl;
      std::cout<<"xp: "<<xp<<std::endl;
      std::cout<<"yp: "<<yp<<std::endl;
      std::cout<<"zp: "<<zp<<std::endl;
      CkExit();
    }
  }
}

void
Performer::advanceParticle( std::size_t i,
                            std::size_t e,
                            const std::array< tk::real, 4>& N)
// *****************************************************************************
// Advance particle based on velocity from mesh cell
//! \author F.J. Gonzalez
// *****************************************************************************
{
  auto dt = g_inputdeck.get< tag::discr, tag::dt >();
  const auto A = m_inpoel[e*4+0];
  const auto B = m_inpoel[e*4+1];
  const auto C = m_inpoel[e*4+2]; 
  const auto D = m_inpoel[e*4+3];
  
  // Update particle coordinates
  // If this condition is true, search for and interpolate the particle
  // velocity using the shape functions defined above. 
  tk::real dvx[4] = { m_u(A,1,0)/m_u(A,0,0) - m_up(A,1,0)/m_up(A,0,0),
                      m_u(B,1,0)/m_u(B,0,0) - m_up(B,1,0)/m_up(B,0,0),
                      m_u(C,1,0)/m_u(C,0,0) - m_up(C,1,0)/m_up(C,0,0),
                      m_u(D,1,0)/m_u(D,0,0) - m_up(D,1,0)/m_up(D,0,0) };
  tk::real dvy[4] = { m_u(A,2,0)/m_u(A,0,0) - m_up(A,2,0)/m_up(A,0,0),
                      m_u(B,2,0)/m_u(B,0,0) - m_up(B,2,0)/m_up(B,0,0),
                      m_u(C,2,0)/m_u(C,0,0) - m_up(C,2,0)/m_up(C,0,0),
                      m_u(D,2,0)/m_u(D,0,0) - m_up(D,2,0)/m_up(D,0,0) };
  tk::real dvz[4] = { m_u(A,3,0)/m_u(A,0,0) - m_up(A,3,0)/m_up(A,0,0),
                      m_u(B,3,0)/m_u(B,0,0) - m_up(B,3,0)/m_up(B,0,0),
                      m_u(C,3,0)/m_u(C,0,0) - m_up(C,3,0)/m_up(C,0,0),
                      m_u(D,3,0)/m_u(D,0,0) - m_up(D,3,0)/m_up(D,0,0) };
        
  m_particles( i, 0, 0) +=
    dt*(N[0]*dvx[0] + N[1]*dvx[1] + N[2]*dvx[2] + N[3]*dvx[3]);
  m_particles( i, 1, 0) +=
    dt*(N[0]*dvy[0] + N[1]*dvy[1] + N[2]*dvy[2] + N[3]*dvy[3]);
  m_particles( i, 2, 0) +=
    dt*(N[0]*dvz[0] + N[1]*dvz[1] + N[2]*dvz[2] + N[3]*dvz[3]);
  
  applyParBC( i );
}

void
Performer::applyParBC( std::size_t i )
// *****************************************************************************
// Apply boundary conditions to particles
//! \author F.J. Gonzalez
// *****************************************************************************
{
  auto& x = m_particles(i,0,0);
  auto& y = m_particles(i,1,0);
  auto& z = m_particles(i,2,0);

  if ( z > 0 ) z -= 40;
  if ( z < -40 ) z += 40;
  if ( y < -7.5 ) y = -7.5 - (y + 7.5); 
  if ( y > 7.5 ) y = 7.5 - (y - 7.5);
  if ( x < -0.125 ) y = -0.125 - (y + 0.125); 
  if ( x > 0.125 ) y = 0.125 - (y - 0.125);
  // Cylinder boundary conditions
  if ( sqrt( y*y + (z+10.5)*(x+10.5) ) <= 0.5 ) z -= 0.5;
  
}

void
Performer::updateSolution( const std::vector< std::size_t >& gid,
                           const std::vector< tk::real >& u )
// *****************************************************************************
// Update solution vector
//! \param[in] gid Global row indices of the vector updated
//! \param[in] u Portion of the unknown/solution vector updated
//! \author J. Bakosi
// *****************************************************************************
{
  auto ncomp = g_inputdeck.get< tag::component >().nprop();
  Assert( gid.size() * ncomp == u.size(),
          "Size of row ID vector times the number of scalar components and the "
          "size of the solution vector must equal" );

  // Receive update of solution vector
  for (std::size_t i=0; i<gid.size(); ++i) {
    auto id = tk::cref_find( m_lid, gid[i] );
    for (ncomp_t c=0; c<ncomp; ++c) m_un( id, c, 0 ) = u[ i*ncomp+c ];
  }

  // Count number of solution nodes updated
  m_nsol += gid.size();

  // If all contributions we own have been received, continue by updating a
  // different solution vector depending on time step stage
  if (m_nsol == m_gid.size()) {

    if (m_stage < 1) {

      m_uf = m_un;

    } else {

      m_up = m_u;
      m_u = m_un;
      
      track();

      // Optionally output field data
      if (!((m_it+1) % g_inputdeck.get< tag::interval, tag::field >()))
        writeFields( m_it+1, m_t + g_inputdeck.get< tag::discr, tag::dt >() );
      else // if no fields at this time, still signal back to host
        if (g_inputdeck.get< tag::param, tag::compns, tag::npar >() > 0)
          contribute(
            CkCallback( CkReductionTarget( Conductor, parcomplete ),
                        m_conductor ) );

      // Optionally contribute diagnostics, e.g., residuals, back to host
      if (!((m_it+1) % g_inputdeck.get< tag::interval, tag::diag >()))
        m_linsysmerger.ckLocalBranch()->diagnostics();
      else // if no diagnostics at this time, still signal back to host
        contribute(
          CkCallback( CkReductionTarget( Conductor, diagcomplete ),
                      m_conductor ) );
    }

    // Prepare for next time step stage
    m_nsol = 0;

    // Tell the Charm++ runtime system to call back to Conductor::evaluateTime()
    // once all Performer chares have received the update. The reduction is done
    // via creating a callback that invokes the typed reduction client, where
    // m_conductor is the proxy on which the reduction target method,
    // evaluateTime(), is called upon completion of the reduction.
    contribute(
      CkCallback( CkReductionTarget( Conductor, evaluateTime ), m_conductor ) );

//     // TEST FEATURE: Manually migrate this chare by using migrateMe to see if
//     // all relevant state variables are being PUPed correctly.
//     //CkPrintf("I'm performer chare %d on PE %d\n",thisIndex,CkMyPe());
//     if (thisIndex == 2 && CkMyPe() == 2) {
//       /*int j;
//       for (int i; i < 50*std::pow(thisIndex,4); i++) {
//         j = i*thisIndex;
//       }*/
//       CkPrintf("I'm performer chare %d on PE %d\n",thisIndex,CkMyPe());
//       migrateMe(1);
//    }
//    if (thisIndex == 2 && CkMyPe() == 1) {
//      CkPrintf("I'm performer chare %d on PE %d\n",thisIndex,CkMyPe());
//      migrateMe(2);
//    }

  }
}

#include "NoWarning/performer.def.h"
