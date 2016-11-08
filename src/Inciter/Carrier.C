// *****************************************************************************
/*!
  \file      src/Inciter/Carrier.C
  \author    J. Bakosi
  \date      Tue 08 Nov 2016 09:05:24 AM MST
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
#include <algorithm>

#include "Carrier.h"
#include "LinSysMerger.h"
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

// Force the compiler to not instantiate the template below as it is
// instantiated in LinSys/LinSysMerger.C (only required on mac)
extern template class tk::LinSysMerger< inciter::CProxy_Transporter,
                                        inciter::CProxy_Carrier,
                                        inciter::AuxSolverLumpMassDiff >;

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
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

} // inciter::

using inciter::Carrier;

Carrier::Carrier( const TransporterProxy& transporter,
                  const LinSysMergerProxy& lsm,
                  const ParticleWriterProxy& pw,
                  const std::vector< std::size_t >& conn,
                  const std::unordered_map< int,
                          std::vector< std::size_t > >& msum,
                  const std::unordered_map< std::size_t, std::size_t >& cid,
                  int ncarr ) :
  __dep(),
  m_it( 0 ),
  m_itf( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_dt( g_inputdeck.get< tag::discr, tag::dt >() ),
  m_stage( 0 ),
  m_nvol( 0 ),
  m_nhsol( 0 ),
  m_nlsol( 0 ),
  m_naec( 0 ),
  m_nalw( 0 ),
  m_nlim( 0 ),
  m_ncarr( static_cast< std::size_t >( ncarr ) ),
  m_outFilename( g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "." +
                 std::to_string( thisIndex ) ),
  m_transporter( transporter ),
  m_linsysmerger( lsm ),
  m_particlewriter( pw ),
  m_cid( cid ),
  m_el( tk::global2local( conn ) ),     // fills m_inpoel, m_gid, m_lid
  m_coord(),
  m_fluxcorrector( m_inpoel.size() ),
  m_psup( tk::genPsup( m_inpoel, 4, tk::genEsup(m_inpoel,4) ) ),
  m_u( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_ul( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_uf( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_ulf( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_du( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_dul( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_p( m_u.nunk(), m_u.nprop()*2 ),
  m_q( m_u.nunk(), m_u.nprop()*2 ),
  m_a( m_gid.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_lhsd( m_psup.second.size()-1, g_inputdeck.get< tag::component >().nprop() ),
  m_lhso( m_psup.first.size(), g_inputdeck.get< tag::component >().nprop() ),
  m_msum( msum ),
  m_vol( m_gid.size(), 0.0 ),
  m_bid(),
  m_pc(),
  m_qc(),
  m_ac(),
                                                        // 0 = no particles
  m_tracker( g_inputdeck.get< tag::cmd, tag::feedback >(), 0, m_inpoel ),
  m_bc()
// *****************************************************************************
//  Constructor
//! \param[in] transporter Host (Transporter) proxy
//! \param[in] lsm Linear system merger (LinSysMerger) proxy
//! \param[in] pw Particle writer proxy
//! \param[in] conn Vector of mesh element connectivity owned (global IDs)
//! \param[in] msum Global mesh node IDs associated to chare IDs bordering the
//!   mesh chunk we operate on
//! \param[in] cid Map associating old node IDs (as in file) to new node IDs (as
//!   in producing contiguous-row-id linear system contributions)
//! \param[in] ncarr Total number of Carrier chares
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( m_psup.second.size()-1 == m_gid.size(),
          "Number of mesh points and number of global IDs unequal" );

  // Register ourselves with the linear system merger
  m_linsysmerger.ckLocalBranch()->checkin();

  // Count the number of mesh nodes at which we receive data from other chares
  // and compute map associating boundary-chare node ID associated to global ID
  std::vector< std::size_t > c;
  for (const auto& n : m_msum) for (auto i : n.second) c.push_back( i );
  tk::unique( c );
  m_bid = tk::assignLid( c );

  // Allocate receive buffers for FCT
  m_pc.resize( m_bid.size() );
  for (auto& b : m_pc) b.resize( m_u.nprop()*2 );
  m_qc.resize( m_bid.size() );
  for (auto& b : m_qc) b.resize( m_u.nprop()*2 );
  m_ac.resize( m_bid.size() );
  for (auto& b : m_ac) b.resize( m_u.nprop() );
}

void
Carrier::vol()
// *****************************************************************************
//  Read mesh node coordinates, sum mesh volumes to nodes, and start
//  communicating them on chare-boundaries
//! \author J. Bakosi
// *****************************************************************************
{
  // Read coordinates of nodes of the mesh chunk we operate on
  readCoords();

  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  // Compute nodal volumes on our chunk of the mesh
  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ m_inpoel[e*4+0], m_inpoel[e*4+1],
                                           m_inpoel[e*4+2], m_inpoel[e*4+3] }};
    // compute element Jacobi determinant * 5/120 = element volume / 4
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da ) * 5.0 / 120.0;
    // scatter add V/4 to nodes
    for (std::size_t j=0; j<4; ++j) m_vol[N[j]] += J;
  }

  // Send our nodal volume contributions to neighbor chares
  if (m_msum.empty())
    contribute(
       CkCallback(CkReductionTarget(Transporter,volcomplete), m_transporter) );
  else
    for (const auto& n : m_msum) {
      std::vector< tk::real > v;
      for (auto i : n.second) v.push_back( m_vol[ tk::cref_find(m_lid,i) ] );
      thisProxy[ n.first ].comvol( n.second, v );
    }
}

void
Carrier::comvol( const std::vector< std::size_t >& gid,
                 const std::vector< tk::real >& V )
// *****************************************************************************
//  Receive nodal volumes on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive volume contributions
//! \param[in] V Partial sums of nodal volume contributions to chare-boundary
//!   nodes
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( V.size() == gid.size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto lid = tk::cref_find( m_lid, gid[i] );
    Assert( lid < m_vol.size(), "Indexing out of bounds" );
    m_vol[ lid ] += V[i];
  }

  if (++m_nvol == m_msum.size()) {
    m_nvol = 0;
    contribute(
       CkCallback(CkReductionTarget(Transporter,volcomplete), m_transporter) );
  }
}

void
Carrier::setup()
// *****************************************************************************
// Setup rows, query boundary conditions, generate particles, output mesh, etc.
//! \author J. Bakosi
// *****************************************************************************
{
  // Send off global row IDs to linear system merger
  m_linsysmerger.ckLocalBranch()->charerow( thisIndex, m_gid );
  // Send node IDs from element side sets matched to BC set IDs
  bc();
  // Generate particles
  m_tracker.genpar( m_coord, m_inpoel, m_ncarr, thisIndex );
  // Output chare mesh to file
  writeMesh();
  // Output fields metadata to output file
  writeMeta();
}

void
Carrier::bc()
// *****************************************************************************
//  Extract node IDs from side set node lists and match to user-specified
//  boundary conditions
//! \details Boundary conditions (BC), mathematically speaking, are applied on
//!   finite surfaces. These finite surfaces are given by element sets. This
//!   function queries the node lists associated to side set IDs as read in from
//!   file (old mesh node IDs as in file associated to all side sets in file).
//!   Then the user-specified boundary conditions, their values and which side
//!   set they are assigned to, are interrogated and only those nodes and their
//!   BC values are extracted that we operate on. This BC data structure is then
//!   sent to the linear system merger which needs to know about this to apply
//!   BCs before a linear solve.
//! \author J. Bakosi
// *****************************************************************************
{
  // Access all side sets from LinSysMerger
  auto& side = m_linsysmerger.ckLocalBranch()->side();

  // Invert m_cid, a map associating old node IDs (as in file) to new node IDs
  // (as in producing contiguous-row-id linear system contributions), so we can
  // search more efficiently for old node IDs.
  decltype(m_cid) rcid;
  for (const auto& i : m_cid) rcid[ i.second ] = i.first;

  // lambda to find out if we own the old (as in file) global node id
  auto own = [ &rcid ]( std::size_t id ) -> std::pair< bool, std::size_t > {
    auto it = rcid.find( id );
    if (it != end(rcid)) return { true, it->second }; else return { false, 0 };
  };

  // lambda to query the Dirichlet BCs on a side set for all components of all
  // the PDEs integrated
  auto bcval = []( int sideset ) {
    // storage for whether the a BC is set and its value for all components of
    // all PDEs configured to be solved
    std::vector< std::pair< bool, tk::real > > b;
    for (const auto& eq : g_pdes) {
      auto e = eq.dirbc( sideset );
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

  // Collect mesh node IDs we contribute to at which a Dirichlet BC is set. The
  // data structure in which this information is stored associates a vector of
  // pairs of bool and BC value to new global mesh node IDs. 'New' as in
  // producing contiguous-row-id linear system contributions, see also
  // Partitioner.h. The bool indicates whether the BC value is set at the given
  // node by the user. The size of the vectors is the number of PDEs integrated
  // times the number of scalar components in all PDEs. The vector is associated
  // to global node IDs at which the boundary condition will be set. If a node
  // gets boundary conditions from multiple side sets, true values (the fact
  // that a BC is set) overwrite existing false ones, i.e., the union of the
  // boundary conditions will be applied at the node.
  for (const auto& s : side) {
    // get BC values for side set
    auto b = bcval( s.first );
    // query if any BC is to be set on the side set for any component of any of
    // the PDEs integrated
    bool ubc = std::any_of( b.cbegin(), b.cend(),
                            [](const std::pair< bool, tk::real >& p)
                            { return p.first; } );
    // for all node IDs we contribute to in the side set associate BC values and
    // whether they are set for a component (for all PDE components)
    for (auto o : s.second) {
      auto n = own(o);
      if (n.first && inset(n.second,s.second) && ubc) {
        auto& v = m_bc[ n.second ];
        v.resize( b.size() );
        for (std::size_t i=0; i<b.size(); ++i) if (b[i].first) v[i] = b[i];
      }
    }
  }

  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() )
    m_transporter.chbcmatched();

  // Send off list of owned node IDs mapped to side sets to LinSysMerger
  m_linsysmerger.ckLocalBranch()->charebc( m_bc );
}

void
Carrier::init()
// *****************************************************************************
// Set ICs, compute initial time step size, output initial field data, compute
// left-hand-side matrix
//! \author J. Bakosi
// *****************************************************************************
{
  // zero initial solution vector
  m_du.fill( 0.0 );

  // Send off initial guess for assembly
  m_linsysmerger.ckLocalBranch()->charesol( thisIndex, m_gid, m_du );

  // Set initial and boundary conditions for all PDEs
  for (const auto& eq : g_pdes) eq.initialize( m_coord, m_u, m_t, m_gid, m_bc );

  // Compute initial time step size
  dt();

  // Output initial conditions to file (time = 0.0)
  if ( !g_inputdeck.get< tag::cmd, tag::benchmark >() ) writeFields( 0.0 );

  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() )
    m_transporter.chic();

  // Compute left-hand side of PDE
  lhs();
}

void
Carrier::dt()
// *****************************************************************************
// Start computing minimum time step size
//! \author J. Bakosi
// *****************************************************************************
{
  tk::real mindt = std::numeric_limits< tk::real >::max();

  auto const_dt = g_inputdeck.get< tag::discr, tag::dt >();
  auto def_const_dt = g_inputdeck_defaults.get< tag::discr, tag::dt >();
  auto eps = std::numeric_limits< tk::real >::epsilon();

  // use constant dt if configured
  if (std::abs(const_dt - def_const_dt) > eps) {

    mindt = const_dt;

  } else {      // compute dt based on CFL

    // find the minimum dt across all PDEs integrated
    for (const auto& eq : g_pdes) {
      auto eqdt = eq.dt( m_coord, m_inpoel, m_u );
      if (eqdt < mindt) mindt = eqdt;
    }

    // Scale smallest dt with CFL coefficient
    mindt *= g_inputdeck.get< tag::discr, tag::cfl >();

  }

  // Contribute to mindt across all Carrier chares
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
    CkCallback(CkReductionTarget(Transporter,dt), m_transporter) );
}

void
Carrier::lhs()
// *****************************************************************************
// Compute left-hand side of transport equations
//! \author J. Bakosi
// *****************************************************************************
{
  // Compute left-hand side matrix for all equations
  for (const auto& eq : g_pdes)
    eq.lhs( m_coord, m_inpoel, m_psup, m_lhsd, m_lhso );

  // Send off left hand side for assembly
  m_linsysmerger.ckLocalBranch()->
    charelhs( thisIndex, m_gid, m_psup, m_lhsd, m_lhso );

  // Compute lumped mass lhs required for the low order solution
  auto lump = m_fluxcorrector.lump( m_coord, m_inpoel );
  // Send off lumped mass lhs for assembly
  m_linsysmerger.ckLocalBranch()->chareauxlhs( thisIndex, m_gid, lump );

  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() )
    m_transporter.chlhs();

  // Call back to Transporter::initcomplete(), signaling that the initialization
  // is complete and we are now starting time stepping
  contribute(
      CkCallback(CkReductionTarget(Transporter,initcomplete), m_transporter));
}

void
Carrier::rhs( tk::real mult, const tk::Fields& sol )
// *****************************************************************************
// Compute right-hand side of transport equations
//! \param[in] mult Multiplier differentiating the different stages in
//!    multi-stage time stepping
//! \param[in] sol Solution vector at current stage
//! \author J. Bakosi
// *****************************************************************************
{
  // Initialize FCT data structures for new time step stage
  m_p.fill( 0.0 );
  m_a.fill( 0.0 );
  for (std::size_t p=0; p<m_u.nunk(); ++p)
    for (ncomp_t c=0; c<g_inputdeck.get<tag::component>().nprop(); ++c) {
      m_q(p,c*2+0,0) = -std::numeric_limits< tk::real >::max();
      m_q(p,c*2+1,0) = std::numeric_limits< tk::real >::max();
    }

  for (auto& b : m_pc) std::fill( begin(b), end(b), 0.0 );
  for (auto& b : m_ac) std::fill( begin(b), end(b), 0.0 );
  for (auto& b : m_qc)
    for (ncomp_t c=0; c<m_u.nprop(); ++c) {
      b[c*2+0] = -std::numeric_limits< tk::real >::max();
      b[c*2+1] = std::numeric_limits< tk::real >::max();
    }

  // Compute right-hand side vector for all equations
  tk::Fields r( m_gid.size(), g_inputdeck.get< tag::component >().nprop() );
  for (const auto& eq : g_pdes) eq.rhs( mult, m_dt, m_coord, m_inpoel, sol, r );
  // Send off right-hand sides for assembly
  m_linsysmerger.ckLocalBranch()->charerhs( thisIndex, m_gid, r );

  // Compute mass diffusion rhs contribution required for the low order solution
  auto diff = m_fluxcorrector.diff( m_coord, m_inpoel, sol );
  // Send off mass diffusion rhs contribution for assembly
  m_linsysmerger.ckLocalBranch()->chareauxrhs( thisIndex, m_gid, diff );

  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() )
    m_transporter.chrhs();
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
  if (!g_inputdeck.get< tag::cmd, tag::benchmark >()) {
    // Create ExodusII writer
    tk::ExodusIIMeshWriter ew( m_outFilename, tk::ExoWriter::CREATE );
    // Write chare mesh initializing element connectivity and point coords
    ew.writeMesh( tk::UnsMesh( m_inpoel, m_coord ) );
  }
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
  if (!g_inputdeck.get< tag::cmd, tag::benchmark >()) {

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
  auto u = m_u;   // make a copy as eq::output() is allowed to overwrite its arg
  std::vector< std::vector< tk::real > > output;
  for (const auto& eq : g_pdes) {
    auto o = eq.output( time, m_coord, u );
    output.insert( end(output), begin(o), end(o) );
  }
  // Write node fields
  writeSolution( ew, m_itf, output );
}

void
Carrier::doWriteParticles()
// *****************************************************************************
// Output particles fields to file
//! \author J. Bakosi
// *****************************************************************************
{
  if (!g_inputdeck.get< tag::cmd, tag::benchmark >())
    m_tracker.doWriteParticles( m_particlewriter, m_it );
}

void
Carrier::aec( const tk::Fields& Un )
// *****************************************************************************
//  Compute and sum antidiffusive element contributions (AEC) to mesh nodes
//! \details This function computes and starts communicating m_p, which stores
//!    the sum of all positive (negative) antidiffusive element contributions to
//!    nodes (Lohner: P^{+,-}_i), see also FluxCorrector::aec().
//! \author J. Bakosi
// *****************************************************************************
{
  // Compute and sum antidiffusive element contributions to mesh nodes. Note
  // that the sums are complete on nodes that are not shared with other chares
  // and only partial sums on chare-boundary nodes.
  m_fluxcorrector.aec( m_coord, m_inpoel, m_vol, m_bc, m_du, Un, m_p );
  ownaec_complete();
  #ifndef NDEBUG
  ownaec_complete();
  #endif

  if (m_msum.empty())
    comaec_complete();
  else // send contributions to chare-boundary nodes to fellow chares
    for (const auto& n : m_msum) {
      std::vector< std::vector< tk::real > > p;
      for (auto i : n.second) p.push_back( m_p[ tk::cref_find(m_lid,i) ] );
      thisProxy[ n.first ].comaec( n.second, p );
    }
}

void
Carrier::comaec( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& P )
// *****************************************************************************
//  Receive sums of antidiffusive element contributions on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive AEC contributions
//! \param[in] P Partial sums of positive (negative) antidiffusive element
//!   contributions to chare-boundary nodes
//! \details This function receives contributions to m_p, which stores the
//!   sum of all positive (negative) antidiffusive element contributions to
//!   nodes (Lohner: P^{+,-}_i), see also FluxCorrector::aec(). While m_p stores
//!   own contributions, m_pc collects the neighbor chare contributions during
//!   communication. This way work on m_p and m_pc is overlapped. The two are
//!   combined in lim().
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( P.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto bid = tk::cref_find( m_bid, gid[i] );
    Assert( bid < m_pc.size(), "Indexing out of bounds" );
    m_pc[ bid ] += P[i];
  }

  if (++m_naec == m_msum.size()) {
    m_naec = 0;
    comaec_complete();
  }
}

void
Carrier::alw( const tk::Fields& Un, const tk::Fields& Ul )
// *****************************************************************************
//  Compute the maximum and minimum unknowns of elements surrounding nodes
//! \param[in] Un Solution at previous time step stage
//! \param[in] Ul Low order solution
//! \details This function computes and starts communicating m_q, which stores
//!    the maximum and mimimum unknowns of all elements surrounding each node
//!    (Lohner: u^{max,min}_i), see also FluxCorrector::alw().
//! \author J. Bakosi
// *****************************************************************************
{
  // Compute the maximum and minimum unknowns of all elements surrounding nodes
  // Note that the maximum and minimum unknowns are complete on nodes that are
  // not shared with other chares and only partially complete on chare-boundary
  // nodes.
  m_fluxcorrector.alw( m_inpoel, Un, Ul, m_q );
  ownalw_complete();
  #ifndef NDEBUG
  ownalw_complete();
  #endif

  if (m_msum.empty())
    comalw_complete();
  else // send contributions at chare-boundary nodes to fellow chares
    for (const auto& n : m_msum) {
      std::vector< std::vector< tk::real > > q;
      for (auto i : n.second) q.push_back( m_q[ tk::cref_find(m_lid,i) ] );
      thisProxy[ n.first ].comalw( n.second, q );
    }
}

void
Carrier::comalw( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& Q )
// *****************************************************************************
// Receive contributions to the maxima and minima of unknowns of all elements
// surrounding mesh nodes on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] Q Partial contributions to maximum and minimum unknowns of all
//!   elements surrounding nodes to chare-boundary nodes
//! \details This function receives contributions to m_q, which stores the
//!   maximum and mimimum unknowns of all elements surrounding each node
//!   (Lohner: u^{max,min}_i), see also FluxCorrector::alw(). While m_q stores
//!   own contributions, m_qc collects the neighbor chare contributions during
//!   communication. This way work on m_q and m_qc is overlapped. The two are
//!   combined in lim().
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( Q.size() == gid.size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto bid = tk::cref_find( m_bid, gid[i] );
    Assert( bid < m_qc.size(), "Indexing out of bounds" );
    auto& o = m_qc[ bid ];
    const auto& q = Q[i];
    for (ncomp_t c=0; c<m_u.nprop(); ++c) {
      if (q[c*2+0] > o[c*2+0]) o[c*2+0] = q[c*2+0];
      if (q[c*2+1] < o[c*2+1]) o[c*2+1] = q[c*2+1];
    }
  }

  if (++m_nalw == m_msum.size()) {
    m_nalw = 0;
    comalw_complete();
  }
}

void
Carrier::lim()
// *****************************************************************************
//  Compute the limited antidiffusive element contributions
//! \details This function computes and starts communicating m_a, which stores
//!   the limited antidiffusive element contributions assembled to nodes
//!   (Lohner: AEC^c), see also FluxCorrector::limit().
//! \author J. Bakosi
// *****************************************************************************
{
  // Combine own and communicated contributions to P and Q
  for (const auto& b : m_bid) {
    auto lid = tk::cref_find( m_lid, b.first );
    const auto& bpc = m_pc[ b.second ];
    const auto& bqc = m_qc[ b.second ];
    for (ncomp_t c=0; c<m_p.nprop()/2; ++c) {
      m_p(lid,c*2+0,0) += bpc[c*2+0];
      m_p(lid,c*2+1,0) += bpc[c*2+1];
      if (bqc[c*2+0] > m_q(lid,c*2+0,0)) m_q(lid,c*2+0,0) = bqc[c*2+0];
      if (bqc[c*2+1] < m_q(lid,c*2+1,0)) m_q(lid,c*2+1,0) = bqc[c*2+1];
    }
  }

  m_fluxcorrector.lim( m_inpoel, m_bc, m_p, (m_stage<1?m_ulf:m_ul), m_q, m_a );
  ownlim_complete();

  if (m_msum.empty())
    comlim_complete();
  else // send contributions to chare-boundary nodes to fellow chares
    for (const auto& n : m_msum) {
      std::vector< std::vector< tk::real > > a;
      for (auto i : n.second) a.push_back( m_a[ tk::cref_find(m_lid,i) ] );
      thisProxy[ n.first ].comlim( n.second, a );
    }
}

void
Carrier::comlim( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& A )
// *****************************************************************************
//  Receive contributions of limited antidiffusive element contributions on
//  chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] A Partial contributions to antidiffusive element contributions to
//!   chare-boundary nodes
//! \details This function receives contributions to m_a, which stores the
//!   limited antidiffusive element contributions assembled to nodes (Lohner:
//!   AEC^c), see also FluxCorrector::limit(). While m_a stores own
//!   contributions, m_ac collects the neighbor chare contributions during
//!   communication. This way work on m_a and m_ac is overlapped. The two are
//!   combined in apply().
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( A.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto bid = tk::cref_find( m_bid, gid[i] );
    Assert( bid < m_ac.size(), "Indexing out of bounds" );
    m_ac[ bid ] += A[i];
  }
 
  if (++m_nlim == m_msum.size()) {
    m_nlim = 0;
    comlim_complete();
  }
}

void
Carrier::advance( uint8_t stage, tk::real newdt, uint64_t it, tk::real t )
// *****************************************************************************
// Advance equations to next stage in multi-stage time stepping
//! \param[in] stage Stage in multi-stage time stepping
//! \param[in] newdt Size of this new time step
//! \param[in] it Iteration count
//! \param[in] t Physical time
//! \author J. Bakosi
// *****************************************************************************
{
  // Update local copy of time step stage, physical time and time step size, and
  // iteration count
  m_stage = stage;
  m_t = t;
  m_dt = newdt;
  m_it = it;

  // Activate SDAG-waits
  #ifndef NDEBUG
  wait4ver();
  #endif
  wait4fct();
  wait4app();

  // Advance stage in multi-stage time stepping by updating the rhs
  if (m_stage < 1)
    rhs( 0.5, m_u );
  else
    rhs( 1.0, m_uf );
}

void
Carrier::out()
// *****************************************************************************
// Output mesh and particle fields
//! \author J. Bakosi
// *****************************************************************************
{
  // Optionally output field and particle data
  if ( m_stage == 1 &&
       !(m_it % g_inputdeck.get< tag::interval, tag::field >()) &&
       !g_inputdeck.get< tag::cmd, tag::benchmark >() )
  {
    writeFields( m_t+m_dt );
    m_tracker.writeParticles( m_transporter, m_particlewriter, this );
  } else
    contribute(
       CkCallback(CkReductionTarget(Transporter,outcomplete), m_transporter) );
}

void
Carrier::updateLowSol( const std::vector< std::size_t >& gid,
                       const std::vector< tk::real >& du )
// *****************************************************************************
// Update low order solution vector
//! \param[in] gid Global row indices of the vector updated
//! \param[in] du Portion of the unknown/solution vector update
//! \author J. Bakosi
// *****************************************************************************
{
  auto ncomp = g_inputdeck.get< tag::component >().nprop();
  Assert( gid.size() * ncomp == du.size(),
          "Size of row ID vector times the number of scalar components and the "
          "size of the low order solution vector must equal" );

  // Receive update of solution vector
  for (std::size_t i=0; i<gid.size(); ++i) {
    auto id = tk::cref_find( m_lid, gid[i] );
    for (ncomp_t c=0; c<ncomp; ++c) m_dul( id, c, 0 ) = du[ i*ncomp+c ];
  }

  // Count number of solution nodes updated
  m_nlsol += gid.size();

  // If all contributions we own have been received, continue by updating a
  // different solution vector depending on time step stage
  if (m_nlsol == m_gid.size()) {
    m_nlsol = 0;
    if (m_stage < 1) {
      m_ulf = m_u + m_dul;
      alw( m_u, m_ulf );
     } else {
      m_ul = m_u + m_dul;
      alw( m_uf, m_ul );
    }
  }
}

void
Carrier::updateSol( const std::vector< std::size_t >& gid,
                    const std::vector< tk::real >& du )
// *****************************************************************************
// Update high order solution vector
//! \param[in] gid Global row indices of the vector updated
//! \param[in] du Portion of the unknown/solution vector update
//! \author J. Bakosi
// *****************************************************************************
{
  auto ncomp = g_inputdeck.get< tag::component >().nprop();
  Assert( gid.size() * ncomp == du.size(),
          "Size of row ID vector times the number of scalar components and the "
          "size of the high order solution vector must equal" );

  // Receive update of solution vector
  for (std::size_t i=0; i<gid.size(); ++i) {
    auto id = tk::cref_find( m_lid, gid[i] );
    for (ncomp_t c=0; c<ncomp; ++c) m_du( id, c, 0 ) = du[ i*ncomp+c ];
  }

  // Count number of solution nodes updated
  m_nhsol += gid.size();

  // If all contributions we own have been received, continue by updating a
  // different solution vector depending on time step stage
  if (m_nhsol == m_gid.size()) {
    m_nhsol = 0;
    aec( m_stage < 1 ? m_u : m_uf );
  }
}

void
Carrier::verify()
// *****************************************************************************
// Verify antidiffusive element contributions up to linear solver convergence
//! \author J. Bakosi
// *****************************************************************************
{
  if (m_fluxcorrector.verify( m_ncarr, m_inpoel, m_du, m_dul ))
    contribute(
      CkCallback( CkReductionTarget(Transporter,verified), m_transporter) );
}

void
Carrier::diagnostics()
// *****************************************************************************
// Compute diagnostics, e.g., residuals
//! \author J. Bakosi
// *****************************************************************************
{
  // Optionally send latest solution to LinSysMerger for computing diagnostics
  if (m_stage == 1 && !(m_it % g_inputdeck.get< tag::interval, tag::diag >()))
    m_linsysmerger.ckLocalBranch()->charediag( thisIndex, m_gid, m_u );
  else
    contribute(
      CkCallback(CkReductionTarget(Transporter,diagcomplete), m_transporter));
}

std::array< std::array< tk::real, 4 >, 3 >
Carrier::velocity( std::size_t e )
// *****************************************************************************
// Extract velocity at the four cell nodes of a mesh element
//! \param[in] e Element id
//! \return Array of 3 arrays of 4 nodal velocity differences
//! \details This funcion extracts the velocity field, defined by each PDE
//!   integrated. All PDEs configured are interrogated and the last nonzero
//!   velocity vector is returned at all four nodes of the mesh element.
//! \author J. Bakosi
// *****************************************************************************
{
  std::array< std::array< tk::real, 4 >, 3 > c;
  for (const auto& eq : g_pdes) {
    const std::array< std::size_t, 4 > N{{ m_inpoel[e*4+0], m_inpoel[e*4+1],
                                           m_inpoel[e*4+2], m_inpoel[e*4+3] }}; 
    auto v = eq.velocity( m_u, m_coord, N );
    if (!v.empty()) c = std::move(v);
  }

  return c;
}

void
Carrier::apply()
// *****************************************************************************
// Apply limited antidiffusive element contributions
//! \author J. Bakosi
// *****************************************************************************
{
  // Combine own and communicated contributions to A
  for (const auto& b : m_bid) {
    auto lid = tk::cref_find( m_lid, b.first );
    const auto& bac = m_ac[ b.second ];
    for (ncomp_t c=0; c<m_a.nprop(); ++c) m_a(lid,c,0) += bac[c];
  }

  // Apply limited antidiffusive element contributions to low order solution
  if (m_stage < 1)
    m_uf = m_ulf + m_a;   // update half-time solution with limited solution
  else
    m_u = m_ul + m_a;     // update solution with new (limited) solution

  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() )
    m_transporter.chlim();

  // Advance particles
  m_tracker.track( m_transporter, thisProxy, m_coord, m_inpoel, m_msum,
                   thisIndex, this, m_stage, m_dt );

  // Compute diagnostics, e.g., residuals
  diagnostics();

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

#include "NoWarning/carrier.def.h"
