// *****************************************************************************
/*!
  \file      src/Inciter/Performer.C
  \author    J. Bakosi
  \date      Thu 21 Jul 2016 02:49:02 PM MDT
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

#include "Performer.h"
#include "Vector.h"
#include "Reader.h"
#include "ContainerUtil.h"
#include "UnsMesh.h"
#include "Reorder.h"
#include "ExodusIIMeshReader.h"
#include "ExodusIIMeshWriter.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "DerivedData.h"
#include "PDE.h"
#include "Tracker.h"
#include "LinSysMerger.h"

// Force the compiler to not instantiate the template below as it is
// instantiated in LinSys/LinSysMerger.C (only required on mac)
extern template class tk::LinSysMerger< inciter::CProxy_Conductor,
                                        inciter::CProxy_Performer >;

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< PDE > g_pdes;

//! \brief Charm++ BC merger reducer for verifying BCs
//! \details This variable is defined here in the .C file and declared as extern
//!   in Performer.h. If instead one defines it in the header (as static),
//!   a new version of the variable is created any time the header file is
//!   included, yielding no compilation nor linking errors. However, that leads
//!   to runtime errors, since Performer::registerVerifyBCMerger(), a Charm++
//!   "initnode" entry method, *may* fill one while contribute() may use the
//!   other (unregistered) one. Result: undefined behavior, segfault, and
//!   formatting the internet ...
CkReduction::reducerType VerifyBCMerger;

} // inciter::

using inciter::Performer;

Performer::Performer(
  const ConductorProxy& conductor,
  const LinSysMergerProxy& lsm,
  const TrackerProxy& tracker,
  const std::vector< std::size_t >& conn,
  const std::unordered_map< std::size_t, std::size_t >& cid )
:
  m_it( 0 ),
  m_itf( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_stage( 0 ),
  m_nsol( 0 ),
  m_outFilename( g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "." +
                 std::to_string( thisIndex ) ),
  m_conductor( conductor ),
  m_linsysmerger( lsm ),
  m_tracker( tracker ),
  m_cid( cid ),
  m_el( tk::global2local( conn ) ),     // fills m_inpoel and m_gid
  m_lid(),
  m_coord(),
  m_psup( tk::genPsup( m_inpoel, 4, tk::genEsup(m_inpoel,4) ) ),
  m_u( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_uf( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_un( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_lhsd( m_psup.second.size()-1, g_inputdeck.get< tag::component >().nprop() ),
  m_lhso( m_psup.first.size(), g_inputdeck.get< tag::component >().nprop() )
// *****************************************************************************
//  Constructor
//! \param[in] conductor Host (Conductor) proxy
//! \param[in] lsm Linear system merger (LinSysMerger) proxy
//! \param[in] conn Vector of mesh element connectivity owned (global IDs)
//! \param[in] cid Map associating old node IDs (as in file) to new node IDs (as
//!   in producing contiguous-row-id linear system contributions)
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( m_psup.second.size()-1 == m_gid.size(),
          "Number of mesh points and number of global IDs unequal" );

  // Register ourselves with the linear system merger
  m_linsysmerger.ckLocalBranch()->checkin();
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
  genPar( 10 );
  // Output chare mesh to file
  writeMesh();
  // Output mesh-based fields metadata to file
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

  // Output initial conditions to file (it = 1, time = 0.0)
  writeFields( m_t );

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
Performer::recPartLoc( const std::vector< tk::real >& x,
                       const std::vector< tk::real >& y,
                       const std::vector< tk::real >& z )
// *****************************************************************************
// Receive particles from the Tracker 
//! \param[in] pcoord Particle coordinates
//! \author F.J. Gonzalez
// *****************************************************************************
{
  std::cout << "performer id " << thisIndex << ": recPartLoc\n";
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
// Output chare mesh to file
//! \author J. Bakosi
// *****************************************************************************
{
  // Create mesh object initializing element connectivity and point coords
  tk::UnsMesh mesh( m_inpoel, m_coord );

  // Create ExodusII writer
  tk::ExodusIIMeshWriter ew( m_outFilename, tk::ExoWriter::CREATE );

  // Write chare mesh
  ew.writeMesh( mesh );
}

void
Performer::writeChareId( const tk::ExodusIIMeshWriter& ew,
                         uint64_t it ) const
// *****************************************************************************
// Output chare id field to file
//! \param[in] ew ExodusII mesh-based writer object
//! \param[in] it Iteration count
//! \author J. Bakosi
// *****************************************************************************
{
  // Write elem chare id field to mesh
  std::vector< tk::real > chid( m_inpoel.size()/4,
                                static_cast<tk::real>(thisIndex) );
  ew.writeElemScalar( it, 1, chid );
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

  ew.writeElemVarNames( { "Chare Id" } );

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
Performer::writeFields( tk::real time )
// *****************************************************************************
// Output mesh-based fields to file
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

  // Write element fields
  writeChareId( ew, m_itf );

  // Collect node fields output from all PDEs
  m_un = m_u;   // make a copy as eq::output() is allowed to overwrite its arg
  std::vector< std::vector< tk::real > > output;
  for (const auto& eq : g_pdes) {
    auto o = eq.output( time, m_coord, m_un );
    output.insert( end(output), begin(o), end(o) );
  }
  // Write node fields
  writeSolution( ew, m_itf, output );
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
Performer::genPar( std::size_t npar )
// *****************************************************************************
// Update solution vector
//! \param[in] npar Number of particles to generate
// *****************************************************************************
{
  std::vector< tk::real > xp(npar), yp(npar), zp(npar);
  const auto& x = m_coord[0];                                                          const auto& y = m_coord[1];
  const auto& z = m_coord[2];
  for (std::size_t i=0; i<npar; ++i) { 
    tk::real NA=0.1, NB=0.2, NC=0.3, ND = 1-NA-NB-NC;
    if ( std::min(NA,1-NA) > 0 && 
         std::min(NB,1-NB) > 0 && 
         std::min(NC,1-NC) > 0 && 
         std::min(ND,1-ND) > 0 ) {
      for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
        const auto A = m_inpoel[e*4+0];
        const auto B = m_inpoel[e*4+1];
        const auto C = m_inpoel[e*4+2]; 
        const auto D = m_inpoel[e*4+3];
        xp[i] = x[A]*NA + x[B]*NB + x[C]*NC + x[D]*ND;
        yp[i] = y[A]*NA + y[B]*NB + y[C]*NC + y[D]*ND;
        zp[i] = z[A]*NA + z[B]*NB + z[C]*NC + z[D]*ND;
      }
    }
  }

  m_tracker[ thisIndex ].advance( 0.0, m_it, m_t );
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
      m_u = m_un;
      if (!((m_it+1) % g_inputdeck.get< tag::interval, tag::field >()))
        writeFields( m_t + g_inputdeck.get< tag::discr, tag::dt >() );
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
