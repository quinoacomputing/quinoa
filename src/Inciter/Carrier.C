// *****************************************************************************
/*!
  \file      src/Inciter/Carrier.C
  \author    J. Bakosi
  \date      Fri 02 Sep 2016 03:49:19 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Carrier advances a system of transport equations
  \details   Carrier advances a system of transport equations. There are a
    potentially large number of Carrier Charm++ chares created by Transporter.
    Each carrier gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances a PDE in time.
*/
// *****************************************************************************

#include <string>
#include <cmath>
#include <array>
#include <set>

#include <gm19.h>

#include "Carrier.h"
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
extern template class tk::LinSysMerger< inciter::CProxy_Transporter,
                                        inciter::CProxy_Carrier >;

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< PDE > g_pdes;

//! \brief Charm++ reducers used by Carrier
//! \details These variables are defined here in the .C file and declared as
//!   extern in Carrier.h. If instead one defines it in the header (as static),
//!   a new version of the variable is created any time the header file is
//!   included, yielding no compilation nor linking errors. However, that leads
//!   to runtime errors, since Carrier::registerReducers(), a Charm++
//!   "initnode" entry method, *may* fill one while contribute() may use the
//!   other (unregistered) one. Result: undefined behavior, segfault, and
//!   formatting the internet ...
CkReduction::reducerType VerifyBCMerger;
CkReduction::reducerType MeshNodeMerger;

} // inciter::

using inciter::Carrier;

Carrier::Carrier( const TransporterProxy& transporter,
                  const LinSysMergerProxy& lsm,
                  const ParticleWriterProxy& pw,
                  const std::vector< std::size_t >& conn,
                  const std::unordered_map< std::size_t, std::size_t >& cid,
                  int ncarr ) :
  m_it( 0 ),
  m_itf( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_stage( 0 ),
  m_nsol( 0 ),
  m_nchpar( 0 ),
  m_ncarr( static_cast< std::size_t >( ncarr ) ),
  m_outFilename( g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "." +
                 std::to_string( thisIndex ) ),
  m_transporter( transporter ),
  m_linsysmerger( lsm ),
  m_particlewriter( pw ),
  m_cid( cid ),
  m_el( tk::global2local( conn ) ),     // fills m_inpoel and m_gid
  m_lid(),
  m_coord(),
  m_fluxcorrector( tk::genEsup(m_inpoel,4) ),
  m_psup( tk::genPsup( m_inpoel, 4, m_fluxcorrector.esup() ) ),
  m_esupel( tk::genEsupel( m_inpoel, 4, m_fluxcorrector.esup() ) ),
  m_u( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_uf( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_un( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_up( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_p( m_u.nunk(), m_u.nprop()*2 ),
  m_lhsd( m_psup.second.size()-1, g_inputdeck.get< tag::component >().nprop() ),
  m_lhso( m_psup.first.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_particles( 0 * m_inpoel.size()/4, 3 ),
  m_elp( m_particles.nunk() ),
  m_msum(),
  m_parmiss(),
  m_parelse()
// *****************************************************************************
//  Constructor
//! \param[in] transporter Host (Transporter) proxy
//! \param[in] lsm Linear system merger (LinSysMerger) proxy
//! \param[in] conn Vector of mesh element connectivity owned (global IDs)
//! \param[in] cid Map associating old node IDs (as in file) to new node IDs (as
//!   in producing contiguous-row-id linear system contributions)
//! \param[in] nper Total number of Carrier chares
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( m_psup.second.size()-1 == m_gid.size(),
          "Number of mesh points and number of global IDs unequal" );

  // Register ourselves with the linear system merger
  m_linsysmerger.ckLocalBranch()->checkin();

  // Collect fellow chare IDs of our mesh neighbors. The goal here is to end up
  // with a list of chare IDs on each chare that the chare shares at least a
  // single mesh node ID, i.e., its mesh chunk borders another chare's mesh
  // chunk. Note that this is not necessarily the most efficient way of doing
  // this, since we do an all-to-all reduction from and to all Carrier chares,
  // sending all the global node indices of the mess chunk each Carrier chare
  // holds. This contribute() is the start of this, and msum() is the final
  // step. A better way of doing this is to only send the chare-boundary node
  // indices. That information is available in Partitioner (the host chare
  // group we are created by), if not directly, it can be computed from the
  // output of the partitioner. This data structure (m_msum) is used in two
  // different ways: (1) the list of neighboring chare IDs are used for
  // searching for particles traveling across the mesh that left our mesh chunk,
  // and (2) the list of neighboring mesh node IDs (that lie on the
  // chare-boundary) are used by the FCT algorithm in FluxCorrector.
  std::vector< std::pair< int, std::unordered_set< std::size_t > > >
    meshnodes{ { thisIndex, { begin(m_gid), end(m_gid) } } };
  auto stream = serialize( meshnodes );
  CkCallback cb( CkIndex_Carrier::msum(nullptr), thisProxy );
  contribute( stream.first, stream.second.get(), MeshNodeMerger, cb );
}

void
Carrier::msum( CkReductionMsg* msg )
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
        m_msum[ c.first ].push_back( i );
        break;
      }

  contribute(
     CkCallback( CkReductionTarget(Transporter,msumcomplete), m_transporter ) );
}

void
Carrier::setup()
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
  genpar();
  // Output chare mesh to file
  writeMesh();
  // Output fields metadata to output file
  writeMeta();
}

void
Carrier::setupIds()
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
Carrier::queryBCs()
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
Carrier::sendBCs( const std::vector< std::size_t >& bc )
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
Carrier::old( const std::vector< std::size_t >& newids )
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
Carrier::requestBCs()
// *****************************************************************************
// Request owned node IDs on which a Dirichlet BC is set by the user
//! \details Called from host (Transporter), contributing the result back to host
//! \author J. Bakosi
// *****************************************************************************
{
   auto stream = tk::serialize( old( queryBCs() ) );
   CkCallback cb( CkIndex_Transporter::doverifybc(nullptr), m_transporter );
   contribute( stream.first, stream.second.get(), VerifyBCMerger, cb );
}

void
Carrier::oldID( int frompe, const std::vector< std::size_t >& newids )
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
Carrier::bcval( int frompe, const std::vector< std::size_t >& nodes )
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
    Assert( b.size() == g_inputdeck.get< tag::component >().nprop(),
            "The total number of scalar components of all configured PDEs (" +
            std::to_string( g_inputdeck.get< tag::component >().nprop() ) + ") "
            "and the sum of the lengths of the Dirichlet BC vectors returned "
            "from all configured PDE::dirbc() calls (" +
            std::to_string( b.size() ) + ") does not match" );
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
  // which the boundary condition will be set. If a node gets boundary
  // conditions from multiple side sets, true values (the fact that a BC is set)
  // overwrite existing false ones, i.e., the union of the boundary conditions
  // will be applied at the node.
  std::unordered_map< std::size_t,
                      std::vector< std::pair< bool, tk::real > > > bcv;
  for (const auto& s : side) {
    auto b = bc( s.first );    // query BC values for all components of all PDEs
    for (auto n : nodes)
      if (inset( n, s.second )) {
        auto& v = bcv[n];
        v.resize( b.size() );
        for (std::size_t i=0; i<b.size(); ++i) if (b[i].first) v[i] = b[i];
      }
  }

  m_linsysmerger[ frompe ].charebcval( bcv );
}

void
Carrier::init( tk::real dt )
// *****************************************************************************
// Initialize linear system
//! \author J. Bakosi
// *****************************************************************************
{
  // Set initial conditions for all PDEs
  for (const auto& eq : g_pdes) eq.initialize( m_coord, m_u, m_t );

  // Output initial conditions to file (time = 0.0)
  //out();

  // Send off initial conditions for assembly
  m_linsysmerger.ckLocalBranch()->charesol( thisIndex, m_gid, m_u );

  // Call back to Transporter::initcomplete(), signaling that the initialization
  // is complete and we are now starting time stepping
  contribute(
      CkCallback(CkReductionTarget(Transporter,initcomplete), m_transporter));

  // Compute left-hand side of PDE
  lhs();
  // Start advancing PDE in time at time step stage 0
  advance( 0, dt, m_it, m_t );
}

void
Carrier::lhs()
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
Carrier::rhs( tk::real mult,
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

  // Compute lhs and rhs required for the low order solution
  auto low = m_fluxcorrector.low( m_coord, m_inpoel, sol );
  // Send off lumped mass lhs for assembly
  m_linsysmerger.ckLocalBranch()->charelump( thisIndex, m_gid, low.first );
  // Send off mass diffusion rhs addendum for assembly
  m_linsysmerger.ckLocalBranch()->charediff( thisIndex, m_gid, low.second );
}

void
Carrier::readCoords()
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
Carrier::writeMesh()
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
Carrier::writeSolution( const tk::ExodusIIMeshWriter& ew,
                        uint64_t it,
                        const std::vector< std::vector< tk::real > >& u ) const
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
Carrier::writeMeta() const
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
Carrier::writeFields( tk::real time )
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
Carrier::writeParticles()
// *****************************************************************************
// Output number of particles we will write to file in this step
//! \author J. Bakosi
// *****************************************************************************
{
  // Send number of partciles we will contribute to particle writer
  m_particlewriter.ckLocalBranch()->npar( m_particles.nunk() );
  // Signal back to host that we are done with sending our number of particles
  contribute(
      CkCallback(CkReductionTarget(Transporter,nparcomplete), m_transporter) );
}

void
Carrier::doWriteParticles()
// *****************************************************************************
// Output particles fields to file
//! \author J. Bakosi
// *****************************************************************************
{
  m_particlewriter.ckLocalBranch()->writeCoords( m_it,
                                                 m_particles.extract(0,0),
                                                 m_particles.extract(1,0),
                                                 m_particles.extract(2,0) );
}

void
Carrier::aec()
// *****************************************************************************
//  Compute and sum antidiffusive element contributions to mesh nodes
//! \author J. Bakosi
// *****************************************************************************
{
  // Compute and sum antidiffusive element contributions to mesh nodes. Note
  // that the sums are complete on nodes that are not shared with other chares
  // and only partial sums on chare-boundary nodes.
  m_p = m_fluxcorrector.aec( m_coord, m_inpoel, m_u, m_un );

  // Send contributions to chare-boundary nodes to fellow chares
  for (const auto& n : m_msum) {
    std::vector< std::size_t > gid;
    std::vector< std::vector< tk::real > > p;
    for (auto i : n.second) {
      gid.push_back( i );
      p.push_back( m_p.extract( tk::cref_find(m_lid,i) ) );
    }
    thisProxy[ n.first ].sumaec( thisIndex, gid, p );
  }
}

void
Carrier::sumaec( int fromch,
                 const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& P )
// *****************************************************************************
//  Finish summing antidiffusive element contributions on chare-boundaries
//! \author J. Bakosi
// *****************************************************************************
{
  for (std::size_t i=0; i<gid.size(); ++i) {
    auto lid = tk::cref_find( m_lid, gid[i] );
    Assert( lid > m_p.nunk(),
            "Indexing out of unknowns summing AEC on chare-boundaries" );
    const auto& p = P[i];
    Assert( p.size() == m_p.nprop(),
            "Indexing out of components summing AEC on chare-boundaries" );
    for (ncomp_t c=0; c<m_p.nprop(); ++c) m_p(lid,c,0) += p[c];
  }
}

void
Carrier::advance( uint8_t stage, tk::real dt, uint64_t it, tk::real t )
// *****************************************************************************
// Advance equations to next stage in multi-stage time stepping
//! \param[in] stage Stage in multi-stage time stepping
//! \param[in] dt Size of time step
//! \param[in] it Iteration count
//! \param[in] t Physical time
//! \author J. Bakosi
// *****************************************************************************
{
  // Update local copy of time step stage, hysical time, and iteration count
  m_stage = stage;
  m_t = t;
  m_it = it;

  // Advance stage in multi-stage time stepping by updating the rhs
  if (m_stage < 1)
    rhs( 0.5, dt, m_u, m_uf );
  else
    rhs( 1.0, dt, m_uf, m_un );
}

void
Carrier::genpar()
// *****************************************************************************
// Generate particles to each of our mesh cells
//! \author F.J. Gonzalez
// *****************************************************************************
{
  auto rng = tk::RNGSSE< gm19_state, unsigned, gm19_generate_ >
                       ( static_cast<unsigned>(m_ncarr), gm19_init_sequence_ );

  // Create a reference of mesh point coordinates
  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  Assert( m_elp.size() >= m_particles.nunk(),
          "Element-of-particle array not large enough" );

  // Generate npar number of particles into each mesh cell
  auto npar = m_particles.nunk() / (m_inpoel.size()/4);
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
        m_particles(i,0,0) = x[A]*N[0] + x[B]*N[1] + x[C]*N[2] + x[D]*N[3];
        m_particles(i,1,0) = y[A]*N[0] + y[B]*N[1] + y[C]*N[2] + y[D]*N[3];
        m_particles(i,2,0) = z[A]*N[0] + z[B]*N[1] + z[C]*N[2] + z[D]*N[3];
        m_elp[i] = e;
      } else --p; // retry if particle was not generated into cell
    }
  }
}

bool
Carrier::parinel( std::size_t p, std::size_t e, std::array< tk::real, 4 >& N )
// *****************************************************************************
// Search particle in a single mesh cell
//! \param[in] p Particle index
//! \param[in] e Mesh cell index
//! \param[inout] N Shapefunctions evaluated at the particle position
//! \return True if particle is in mesh cell
//! \author F.J. Gonzalez
// *****************************************************************************
{
  // Tetrahedron node indices
  const auto A = m_inpoel[e*4+0];
  const auto B = m_inpoel[e*4+1];
  const auto C = m_inpoel[e*4+2]; 
  const auto D = m_inpoel[e*4+3];

  // Tetrahedron node coordinates
  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  // Particle coordinates
  const auto& xp = m_particles( p, 0, 0 );
  const auto& yp = m_particles( p, 1, 0 );
  const auto& zp = m_particles( p, 2, 0 );

  // Evaluate linear shapefunctions at particle locations using Cramer's Rule
  //    | xp |   | x1 x2 x3 x4 |   | N1 |
  //    | yp | = | y1 y2 y3 y4 | • | N2 |
  //    | zp |   | z1 z2 z3 z4 |   | N3 |
  //    | 1  |   | 1  1  1  1  |   | N4 |

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

  // Shape functions evaluated at particle location
  N[0] = DetX1/DetX;
  N[1] = DetX2/DetX;
  N[2] = DetX3/DetX;
  N[3] = DetX4/DetX;

  // if min( N^i, 1-N^i ) > 0 for all i, particle is in cell
  if ( std::min(N[0],1.0-N[0]) > 0 && std::min(N[1],1.0-N[1]) > 0 &&
       std::min(N[2],1.0-N[2]) > 0 && std::min(N[3],1.0-N[3]) > 0 )
  {
    m_elp.resize( p+1 );
    m_elp[ p ] = e;
    return true;
  } else {
    return false;
  }
}

std::vector< std::size_t >
Carrier::addpar( const std::vector< std::size_t >& miss,
                 const tk::Particles& ps )
// *****************************************************************************
// Try to find particles and add those found to the list of ours
//! \param[in] miss Indices of particles to find
//! \param[in] ps Particle data associated to those particle indices to find
//! \return Particle indices found
//! \author J. Bakosi
// *****************************************************************************
{
  std::vector< std::size_t > found;     // will store indices of particles found

  // try to find particles received
  for (std::size_t i=0; i<ps.nunk(); ++i)
    for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
      std::array< tk::real, 4 > N;
      auto last = m_particles.nunk();
      m_particles.add( ps[i] );
      if (parinel( last, e, N )) {
        found.push_back( miss[i] );
        m_elp.resize( last+1 );
        m_elp[ last ] = e;
      } else {
        m_particles.rm( { last } );
      }
    }

  return found;
}

void
Carrier::track()
// *****************************************************************************
// Advance our particles and search for their new mesh cells
//! \author F.J. Gonzalez
// *****************************************************************************
{
  // Only advance particles in the finel time step stage
  if (m_stage < 1) {
    contribute( CkCallback( CkReductionTarget( Transporter, parcomcomplete ),
                m_transporter));
    return;
  }

  // Lambda to attempt to find and advance particle i in element e. Returns true
  // if the particle was found (and advanced), false if it was not found.
  std::array< tk::real, 4 > N;
  auto adv = [ this, &N ]( std::size_t i, std::size_t e ) -> bool {
    if (this->parinel( i, e, N )) {
      advanceParticle( i, e, N );
      return true;
    }
    return false;
  };

  // Search cells of our mesh chunk for all particles
  for (std::size_t i=0; i<m_particles.nunk(); ++i) {
    // First try to find particle where it has last been found
    bool found = adv( i, m_elp[i] );
    // Next search in the elements surroundings the points of the element where
    // the particle has last been found
    if (!found) {
      auto last = m_esupel.second[m_elp[i]+1];
      for (auto j=m_esupel.second[m_elp[i]]+1; j<=last; ++j) {
        found = adv( i, m_esupel.first[j] );
        if (found) j = last;  // search for next particle
      }
    }
    // Next search all cells in our chunk of the mesh
    if (!found) {
      for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
        found = adv( i, e );
        if (found) e = m_inpoel.size()/4;  // search for next particle
      }
      // If the particle still has not been found, it left our chunk of the
      // mesh, mark as missing (will initiate communication to find it)
      if (!found) m_parmiss.insert( i );
    }
  }

  // If we have no missing particles, we are done, if we do, send out requests
  // to find them to those Carrier chares with which we neighbor mesh cells
  if (m_parmiss.empty()) {
    m_nchpar = 0;
    m_parmiss.clear();
    m_parelse.clear();
    contribute( CkCallback( CkReductionTarget( Transporter, parcomcomplete ),
                m_transporter));
  } else {
    decltype(m_particles) pexp( m_parmiss.size(), 3 );
    std::size_t p = 0;
    for (auto i : m_parmiss) {
      pexp( p, 0, 0 ) = m_particles( i, 0, 0 );
      pexp( p, 1, 0 ) = m_particles( i, 1, 0 );
      pexp( p, 2, 0 ) = m_particles( i, 2, 0 );
      ++p;
    }
    std::vector< std::size_t > miss( begin(m_parmiss), end(m_parmiss) );
    for (const auto& n : m_msum)
      thisProxy[ n.first ].findpar( thisIndex, miss, pexp );
  }
}

void
Carrier::findpar( int fromch,
                  const std::vector< std::size_t >& miss,
                  const tk::Particles& ps )
// *****************************************************************************
// Find particles missing by the requestor and make those found ours
//! \param[in] fromch Chare ID the request originates from
//! \param[in] miss Indices of particles to find
//! \param[in] ps Particle data associated to those particle indices to find
//! \author J. Bakosi
// *****************************************************************************
{
  // Try to find particles missing by the requestor and make those found ours
  auto found = addpar( miss, ps );

  // Send the particle indices we found back to the requestor
  thisProxy[ fromch ].foundpar( found );
}

void
Carrier::foundpar( const std::vector< std::size_t >& found )
// *****************************************************************************
// Receive particle indices found elsewhere (by fellow neighbors)
//! \param[in] found Indices of particles found
//! \author J. Bakosi
// *****************************************************************************
{
  m_parelse.insert( begin(found), end(found) );

  if (++m_nchpar == m_msum.size()) {    // if we have heard from all neighbors

    // find the particle indices that are still have not been found (by close
    // neighbors)
    decltype(m_parelse) far;
    std::set_difference( begin(m_parmiss), end(m_parmiss),
                         begin(m_parelse), end(m_parelse), 
                         std::inserter( far, begin(far) ) );
    m_parmiss = far;

    // if there are still missing particles (not found by close neighbors we
    // share mesh nodes with), we resort to requesting them to be searched by
    // all Carriers
    if (m_parmiss.empty()) {
      m_nchpar = 0;
      m_parmiss.clear();
      m_parelse.clear();
      contribute( CkCallback( CkReductionTarget( Transporter, parcomcomplete ),
                  m_transporter));
    } else {
      decltype(m_particles) pexp( m_parmiss.size(), 3 );
      std::size_t p = 0;
      for (auto i : m_parmiss) {
        pexp( p, 0, 0 ) = m_particles( i, 0, 0 );
        pexp( p, 1, 0 ) = m_particles( i, 1, 0 );
        pexp( p, 2, 0 ) = m_particles( i, 2, 0 );
        ++p;
      }
      m_particles.rm( m_parmiss );
      m_nchpar = 0;
      m_parelse.clear();
      std::vector< std::size_t > miss( begin(m_parmiss), end(m_parmiss) );
      thisProxy.collectpar( thisIndex, miss, pexp );  // broadcast to everyone
    }
  }
}

void
Carrier::collectpar( int fromch,
                     const std::vector< std::size_t >& miss,
                     const tk::Particles& ps )
// *****************************************************************************
// Find particles missing by the requestor and make those found ours
//! \param[in] fromch Chare ID the request originates from
//! \param[in] miss Indices of particles to find
//! \param[in] ps Particle data associated to those particle indices to find
//! \author J. Bakosi
// *****************************************************************************
{
  // Try to find particles missing by the requestor and make those found ours
  auto found = addpar( miss, ps );

  // Send the particle indices we found back to the requestor
  thisProxy[ fromch ].collectedpar( found );
}

void
Carrier::collectedpar( const std::vector< std::size_t >& found )
// *****************************************************************************
// Collect particle indices found elsewhere (by far fellows)
//! \param[in] found Indices of particles found
//! \author J. Bakosi
// *****************************************************************************
{
  // Collect particle indices found elsewhere (by distant neighbors)
  m_parelse.insert( begin(found), end(found) );

  if (++m_nchpar == m_ncarr) {  // if we have heard from everyone
    Assert( m_parmiss == m_parelse, "Not all particles have been found on "
            "chare " + std::to_string(thisIndex) + '\n' );
    m_nchpar = 0;
    m_parmiss.clear();
    m_parelse.clear();
    contribute( CkCallback( CkReductionTarget( Transporter, parcomcomplete ),
                m_transporter));
  }
}

void
Carrier::advanceParticle( std::size_t i,
                          std::size_t e,
                          const std::array< tk::real, 4 >& Np )
// *****************************************************************************
// Advance particle based on velocity from mesh cell
//! \author F.J. Gonzalez
// *****************************************************************************
{
  auto dt = g_inputdeck.get< tag::discr, tag::dt >();

  // Extract velocity at the four cell nodes at current and previous time step.
  // To keep the code below as general as possible, we interrogate all PDEs
  // configured and use the last nonzero velocity vector.
  std::vector< std::array< tk::real, 4 > > c, p;
  for (const auto& eq : g_pdes) {
    const std::array< std::size_t, 4 > N{{ m_inpoel[e*4+0], m_inpoel[e*4+1],
                                           m_inpoel[e*4+2], m_inpoel[e*4+3] }}; 
    auto v = eq.velocity( m_u, m_coord, N );
    if (!v.empty()) c = std::move(v);
    auto w = eq.velocity( m_up, m_coord, N );
    if (!w.empty()) p = std::move(w);
  }

  Assert(c.size() == 3 && p.size() == 3, "PDE velocity must have 3 components");

  // Compute velocity differnce between current and previous time steps
  const auto& cx = c[0];
  const auto& px = p[0];
  const auto& cy = c[1];
  const auto& py = p[1];
  const auto& cz = c[2];
  const auto& pz = p[2];
  tk::real dvx[4] = { cx[0]-px[0], cx[1]-px[1], cx[2]-px[2], cx[3]-px[3] };
  tk::real dvy[4] = { cy[0]-py[0], cy[1]-py[1], cy[2]-py[2], cy[3]-py[3] };
  tk::real dvz[4] = { cz[0]-pz[0], cz[1]-pz[1], cz[2]-pz[2], cz[3]-pz[3] };

  // Advance particle coordinates using the interpolated velocity
  m_particles( i, 0, 0 ) +=
    dt*(Np[0]*dvx[0] + Np[1]*dvx[1] + Np[2]*dvx[2] + Np[3]*dvx[3]);
  m_particles( i, 1, 0 ) +=
    dt*(Np[0]*dvy[0] + Np[1]*dvy[1] + Np[2]*dvy[2] + Np[3]*dvy[3]);
  m_particles( i, 2, 0 ) +=
    dt*(Np[0]*dvz[0] + Np[1]*dvz[1] + Np[2]*dvz[2] + Np[3]*dvz[3]);

  // Apply boundary conditions to particle
  applyParBC( i );
}

void
Carrier::applyParBC( std::size_t i )
// *****************************************************************************
// Apply boundary conditions to particles
//! \author F.J. Gonzalez
// *****************************************************************************
{
  auto& x = m_particles(i,0,0);
  auto& y = m_particles(i,1,0);
  auto& z = m_particles(i,2,0);

//   if ( z > 0 ) z -= 40;
//   if ( z < -40 ) z += 40;
//   if ( y < -5.0 ) y = -5.0 - (y + 5.0); 
//   if ( y > 5.0 ) y = 5.0 - (y - 5.0);
//   if ( x < -0.25 ) y = -0.25 - (y + 0.25); 
//   if ( x > 0.25 ) y = 0.25 - (y - 0.25);
//   // Cylinder boundary conditions
//   if ( sqrt( y*y + (z+10.5)*(x+10.5) ) <= 0.5 ) z -= 0.5;

  if (z > 1) z = 0.9;
  if (z < 0) z = 0.1;
  if (y > 1) y = 0.9;
  if (y < 0) y = 0.1;
  if (x > 1) x = 0.9;
  if (x < 0) x = 0.1;
}

void
Carrier::out()
// *****************************************************************************
// Output mesh and particle fields
//! \author J. Bakosi
// *****************************************************************************
{
  // Optionally output field and particle data
  if (m_stage == 1 &&
      !(m_it % g_inputdeck.get< tag::interval, tag::field >()))
  {
    writeFields( m_t + g_inputdeck.get< tag::discr, tag::dt >() );
    writeParticles();
  } else
    contribute(
       CkCallback(CkReductionTarget(Transporter,outcomplete), m_transporter) );
}

void
Carrier::updateSolution( const std::vector< std::size_t >& gid,
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

    //aec();

    if (m_stage < 1) {
      m_uf = m_un;
    } else {
      m_up = m_u;       // save time n solution for particle update
      m_u = m_un;
    }

    m_nsol = 0;

    // Advance particles
    track();

    // Optionally contribute diagnostics, e.g., residuals, back to host
    if (m_stage == 1 && !(m_it % g_inputdeck.get< tag::interval, tag::diag >()))
      m_linsysmerger.ckLocalBranch()->diagnostics();
    else
      contribute(
        CkCallback(CkReductionTarget(Transporter,diagcomplete), m_transporter));

//     // TEST FEATURE: Manually migrate this chare by using migrateMe to see if
//     // all relevant state variables are being PUPed correctly.
//     //CkPrintf("I'm carrier chare %d on PE %d\n",thisIndex,CkMyPe());
//     if (thisIndex == 2 && CkMyPe() == 2) {
//       /*int j;
//       for (int i; i < 50*std::pow(thisIndex,4); i++) {
//         j = i*thisIndex;
//       }*/
//       CkPrintf("I'm carrier chare %d on PE %d\n",thisIndex,CkMyPe());
//       migrateMe(1);
//    }
//    if (thisIndex == 2 && CkMyPe() == 1) {
//      CkPrintf("I'm carrier chare %d on PE %d\n",thisIndex,CkMyPe());
//      migrateMe(2);
//    }

  }
}

#include "NoWarning/carrier.def.h"
