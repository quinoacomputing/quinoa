// *****************************************************************************
/*!
  \file      src/Inciter/Carrier.C
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
CkReduction::reducerType PDFMerger;

} // inciter::

using inciter::Carrier;

Carrier::Carrier( const TransporterProxy& transporter,
                  const LinSysMergerProxy& lsm,
                  const ParticleWriterProxy& pw,
                  const std::vector< std::size_t >& conn,
                  const std::unordered_map< int,
                          std::unordered_set< std::size_t > >& msum,
                  const std::unordered_map< std::size_t, std::size_t >&
                          filenodes,
                  const tk::UnsMesh::EdgeNodes& edgenodes,
                  int ncarr ) :
  __dep(),
  m_it( 0 ),
  m_itf( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_dt( g_inputdeck.get< tag::discr, tag::dt >() ),
  m_lastFieldWriteTime( -1.0 ),
  m_stage( 0 ),
  m_nvol( 0 ),
  m_nhsol( 0 ),
  m_nlsol( 0 ),
  m_naec( 0 ),
  m_nalw( 0 ),
  m_nlim( 0 ),
  m_V( 0.0 ),
  m_ncarr( static_cast< std::size_t >( ncarr ) ),
  m_outFilename( g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "." +
                 std::to_string( thisIndex ) ),
  m_transporter( transporter ),
  m_linsysmerger( lsm ),
  m_particlewriter( pw ),
  m_filenodes( filenodes ),
  m_edgenodes( edgenodes ),
  m_el( tk::global2local( conn ) ),     // fills m_inpoel, m_gid, m_lid
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
  m_v( m_gid.size(), 0.0 ),
  m_vol( m_gid.size(), 0.0 ),
  m_bid(),
  m_pc(),
  m_qc(),
  m_ac(),                                               // 0 = no particles
  m_tracker( g_inputdeck.get< tag::cmd, tag::feedback >(), 0, m_inpoel )
// *****************************************************************************
//  Constructor
//! \param[in] transporter Host (Transporter) proxy
//! \param[in] lsm Linear system merger (LinSysMerger) proxy
//! \param[in] pw Particle writer proxy
//! \param[in] conn Vector of mesh element connectivity owned (global IDs)
//! \param[in] msum Global mesh node IDs associated to chare IDs bordering the
//!   mesh chunk we operate on
//! \param[in] filenodes Map associating old node IDs (as in file) to new node
//!   IDs (as in producing contiguous-row-id linear system contributions)
//! \param[in] edgenodes Map associating new node IDs ('new' as in producing
//!   contiguous-row-id linear system contributions) to edges (a pair of old
//!   node IDs ('old' as in file). These 'new' node IDs are the ones newly
//!   added during inital uniform mesh refinement.
//! \param[in] ncarr Total number of Carrier chares
//! \details "Contiguous-row-id" here means that the numbering of the mesh nodes
//!   (which corresponds to rows in the linear system) are (approximately)
//!   contiguous (as much as this can be done with an unstructured mesh) as the
//!   problem is distirbuted across PEs, held by LinSysMerger objects. This
//!   ordering is in start contrast with "as-in-file" ordering, which is the
//!   ordering of the mesh nodes as it is stored in the file from which the mesh
//!   is read in. The as-in-file ordering is highly non-contiguous across the
//!   distributed problem.
// *****************************************************************************
{
  Assert( m_psup.second.size()-1 == m_gid.size(),
          "Number of mesh points and number of global IDs unequal" );

  // Convert neighbor nodes to vectors from sets
  for (const auto& n : msum) {
    auto& v = m_msum[ n.first ];
    v.insert( end(v), begin(n.second), end(n.second) );
  }

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

  // Invert file-node map, a map associating old node IDs (as in file) to new
  // node IDs (as in producing contiguous-row-id linear system contributions),
  // so we can search more efficiently for old node IDs.
  decltype(m_filenodes) linnodes;
  for (const auto& i : m_filenodes) {
    auto n = tk::cref_find(m_lid,i.first);
    Assert( n < m_gid.size(),
            "Local IDs must be lower than the local number of grid points" );
    linnodes[ i.second ] = n;
  }

  // Access all side sets and their old node IDs (as in file) from LinSysMerger
  auto& oldside = m_linsysmerger.ckLocalBranch()->side();

  // Create map that assigns the local mesh node IDs mapped to side set ids,
  // storing only those nodes for a given side set that are part of our chunk of
  // the mesh.
  for (const auto& s : oldside) {
    auto& n = m_side[ s.first ];
    for (auto o : s.second) {
      auto it = linnodes.find( o );
      if (it != end(linnodes))
        n.push_back( it->second );
    }
  }
}

void
Carrier::coord()
// *****************************************************************************
//  Read mesh node coordinates and optionally add new edge-nodes in case of
//  initial uniform refinement
// *****************************************************************************
{
  // Read coordinates of nodes of the mesh chunk we operate on
  readCoords();
  // Add coordinates of mesh nodes newly generated to edge-mid points during
  // initial refinement
  addEdgeNodeCoords();
  // Compute mesh cell volumes
  vol();
  // Compute mesh cell statistics
  stat();
}

void
Carrier::vol()
// *****************************************************************************
// Sum mesh volumes to nodes, start communicating them on chare-boundaries
// *****************************************************************************
{
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
    ErrChk( J > 0, "Element Jacobian non-positive: PE:" +
                   std::to_string(CkMyPe()) + ", node IDs: " +
                   std::to_string(m_gid[N[0]]) + ',' +
                   std::to_string(m_gid[N[1]]) + ',' +
                   std::to_string(m_gid[N[2]]) + ',' +
                   std::to_string(m_gid[N[3]]) + ", coords: (" +
                   std::to_string(x[N[0]]) + ", " +
                   std::to_string(y[N[0]]) + ", " +
                   std::to_string(z[N[0]]) + "), (" +
                   std::to_string(x[N[1]]) + ", " +
                   std::to_string(y[N[1]]) + ", " +
                   std::to_string(z[N[1]]) + "), (" +
                   std::to_string(x[N[2]]) + ", " +
                   std::to_string(y[N[2]]) + ", " +
                   std::to_string(z[N[2]]) + "), (" +
                   std::to_string(x[N[3]]) + ", " +
                   std::to_string(y[N[3]]) + ", " +
                   std::to_string(z[N[3]]) + ')' );
    // scatter add V/4 to nodes
    for (std::size_t j=0; j<4; ++j) { m_v[N[j]] = m_vol[N[j]] += J; }
  }

  // Sum mesh volume to host
  tk::real V = 0.0;
  for (auto v : m_v) V += v;
  contribute( sizeof(tk::real), &V, CkReduction::sum_double,
    CkCallback(CkReductionTarget(Transporter,vol), m_transporter) );

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
Carrier::stat()
// *****************************************************************************
// Compute mesh volume statistics
// *****************************************************************************
{
  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  auto MIN = -std::numeric_limits< tk::real >::max();
  auto MAX = std::numeric_limits< tk::real >::max();
  std::array< tk::real, 2 > min = {{ MAX, MAX }};
  std::array< tk::real, 2 > max = {{ MIN, MIN }};
  std::array< tk::real, 4 > sum{{ 0.0, 0.0, 0.0, 0.0 }};
  tk::UniPDF edgePDF( 1e-4 );
  tk::UniPDF volPDF( 1e-4 );

  // Compute edge length statistics
  // Note that while the min and max edge lengths are independent of the number
  // of CPUs (by the time they are aggregated across all chares), the sum of
  // the edge lengths and the edge length PDF are not. This is because the
  // edges on the chare-boundary are counted multiple times and we
  // conscientiously do not make an effort to precisely compute this, because
  // that would require communication and more complex logic. Since these
  // statistics are intended as simple average diagnostics, we ignore these
  // small differences. For reproducible average edge lengths and edge length
  // PDFs, run the mesh in serial.
  for (std::size_t p=0; p<m_gid.size(); ++p)
    for (auto i=m_psup.second[p]+1; i<=m_psup.second[p+1]; ++i) {
       const auto dx = x[ m_psup.first[i] ] - x[ p ];
       const auto dy = y[ m_psup.first[i] ] - y[ p ];
       const auto dz = z[ m_psup.first[i] ] - z[ p ];
       const auto length = std::sqrt( dx*dx + dy*dy + dz*dz );
       if (length < min[0]) min[0] = length;
       if (length > max[0]) max[0] = length;
       sum[0] += 1.0;
       sum[1] += length;
       edgePDF.add( length );
    }

  // Compute mesh cell volume statistics
  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ m_inpoel[e*4+0], m_inpoel[e*4+1],
                                           m_inpoel[e*4+2], m_inpoel[e*4+3] }};
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto L = std::cbrt( tk::triple( ba, ca, da ) / 6.0 );
    if (L < min[1]) min[1] = L;
    if (L > max[1]) max[1] = L;
    sum[2] += 1.0;
    sum[3] += L;
    volPDF.add( L );
  }

  // Contribute to mesh statistics across all Carrier chares
  contribute( min.size()*sizeof(tk::real), min.data(), CkReduction::min_double,
    CkCallback(CkReductionTarget(Transporter,minstat), m_transporter) );
  contribute( max.size()*sizeof(tk::real), max.data(), CkReduction::max_double,
    CkCallback(CkReductionTarget(Transporter,maxstat), m_transporter) );
  contribute( sum.size()*sizeof(tk::real), sum.data(), CkReduction::sum_double,
    CkCallback(CkReductionTarget(Transporter,sumstat), m_transporter) );

  // Serialize PDFs to raw stream
  auto stream = tk::serialize( { edgePDF, volPDF } );
  // Create Charm++ callback function for reduction of PDFs with
  // Transporter::pdfstat() as the final target where the results will appear.
  CkCallback cb( CkIndex_Transporter::pdfstat(nullptr), m_transporter );
  // Contribute serialized PDF of partial sums to host via Charm++ reduction
  contribute( stream.first, stream.second.get(), PDFMerger, cb );
}

void
Carrier::comvol( const std::vector< std::size_t >& gid,
                 const std::vector< tk::real >& V )
// *****************************************************************************
//  Receive nodal volumes on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive volume contributions
//! \param[in] V Partial sums of nodal volume contributions to chare-boundary
//!   nodes
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
Carrier::setup( tk::real v )
// *****************************************************************************
// Setup rows, query boundary conditions, generate particles, output mesh, etc.
//! \param[in] v Total mesh volume
// *****************************************************************************
{
  // Store total mesh volume
  m_V = v;
  // Output chare mesh to file
  writeMesh();
  // Send off global row IDs to linear system merger
  m_linsysmerger.ckLocalBranch()->charerow( thisIndex, m_gid );
  // Generate particles
  m_tracker.genpar( m_coord, m_inpoel, m_ncarr, thisIndex );
  // Output fields metadata to output file
  writeMeta();
}

void
Carrier::bc()
// *****************************************************************************
//  Extract node IDs from side set node lists and match to user-specified
//  boundary conditions
//! \details Boundary conditions (BC), mathematically speaking, are applied on
//!   finite surfaces. These finite surfaces are given by element sets (i.e., a
//!   list of elements). This function queries Dirichlet boundary condition
//!   values from all PDEs in the system of PDEs integrated at the node lists
//!   associated to side set IDs (previously read from file). As a response to
//!   this query, each PDE system returns a BC data structure which is then
//!   sent to the linear system merger which needs to know about this to apply
//!   BCs before a linear solve. Note that the BC mesh nodes that this function
//!   results in, stored in dirbc and sent to the linear system merger, only
//!   contains those nodes that this chare contributes to, i.e., it does not
//!   contain those BC nodes at which other chares enforce Dirichlet BCs. The
//!   linear system merger then collects these and communicates to other PEs so
//!   that BC data held in LinSysMerger::m_bc are the same on all PEs. That is
// *****************************************************************************
{
  // Vector of pairs of bool and boundary condition value associated to mesh
  // node IDs at which the user has set Dirichlet boundary conditions for all
  // PDEs integrated
  std::unordered_map< std::size_t, NodeBC > dirbc;

  // Query Dirichlet boundary conditions for all PDEs integrated and assign to
  // nodes. This is where the individual system of PDEs are queried for boundary
  // conditions. The outer loop goes through all sides sets that exists in the
  // input file and passes the map's value_type (a pair of the side set id and a
  // vector of local node IDs) to PDE::dirbc(). PDE::dirbc() returns a new map
  // that associates a vector of pairs associated to local node IDs. (The pair
  // is a pair of bool and real value, the former is the fact that the BC is to
  // be set while the latter is the value if it is to be set). The length of
  // this NodeBC vector, returning from each system of PDEs equals to the number
  // of scalar components the given PDE integrates. Here then we contatenate
  // this map for all PDEs integrated. If there are multiple BCs set at a mesh
  // node (dirbc::key), either because (1) in the same PDE system the user
  // prescribed BCs on side sets that share nodes or (2) because more than a
  // single PDE system assigns BCs to a given node (on different variables), the
  // NodeBC vector must be correctly stored. "Correctly" here means that the
  // size of the NodeBC vectors must all be the same and qual to the sum of all
  // scalar components integrated by all PDE systems integrated. Example:
  // single-phase compressible flow (density, momentum, energy = 5) +
  // transported scalars of 10 variables -> NodeBC vector length = 15. Note that
  // in case (1) above a new node encountered must "overwrite" the already
  // existing space for the NodeBC vector. "Overwrite" here means that it should
  // keep the existing BCs and add the new ones yielding the union the two
  // prescription for BCs but in the same space that already exist in the NodeBC
  // vector. In case (2), however, the NodeBC pairs must go to the location in
  // the vector assigned to the given PDE system, i.e., using the above example
  // BCs for the 10 (or less) scalars should go in the positions starting at 5,
  // leaving the first 5 false, indicating no BCs for the flow variables.
  //
  // TODO: Note that the logic described above is only partially implemented at
  // this point. What works is the correct insertion of multiple BCs for nodes
  // shared among multiple side sets, e.g., corners, originating from the same
  // PDE system. What is not yet implemented is the case when there are no BCs
  // set for flow variables but there are BCs for transport, the else branch
  // below will incorrectly NOT skip the space for the flow variables. In other
  // words, this only works for a single PDE system and a sytem of systems. This
  // machinery is only tested with a single system of PDEs at this point.
  for (const auto& s : m_side) {
    std::size_t c = 0;
    for (std::size_t eq=0; eq<g_pdes.size(); ++eq) {
      auto eqbc = g_pdes[eq].dirbc( m_t, m_dt, s, m_coord );
      for (const auto& n : eqbc) {
        auto id = n.first;                      // BC node ID
        const auto& bc = n.second;              // BCs
        auto& nodebc = dirbc[ m_gid[id] ];      // BCs to be set for node
        if (nodebc.size() > c) {        // node already has BCs from this PDE
          Assert( nodebc.size() == c+bc.size(), "Size mismatch" );
          for (std::size_t i=0; i<bc.size(); i++)
            if (bc[i].first)
              nodebc[c+i] = bc[i];
        } else {        // node does not yet have BCs from this PDE
          // This branch needs to be completed for system of systems of PDEs.
          // See note above.
          nodebc.insert( end(nodebc), begin(bc), end(bc) );
        }
      }
      if (!eqbc.empty()) c += eqbc.cbegin()->second.size();
    }
  }

  // Verify the size of each NodeBC vectors. They must have the same lengths and
  // equal to the total number of scalar components for all systems of PDEs
  // integrated. This is intentional, because this way the linear system merger
  // does not have to (and does not) know about individual equation systems.
  for (const auto& n : dirbc) {
    IGNORE(n);
    Assert( n.second.size() == m_u.nprop(), "Size of NodeBC vector incorrect" );
  }

  // Send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() )
    m_transporter.chbcmatched();

  // Send off list of owned node IDs mapped to side sets to LinSysMerger
  m_linsysmerger.ckLocalBranch()->charebc( dirbc );
}

void
Carrier::init()
// *****************************************************************************
// Set ICs, compute initial time step size, output initial field data, compute
// left-hand-side matrix
// *****************************************************************************
{
  // zero initial solution vector
  m_du.fill( 0.0 );

  // Send off initial guess for assembly
  m_linsysmerger.ckLocalBranch()->charesol( thisIndex, m_gid, m_du );

  // Set initial conditions for all PDEs
  for (const auto& eq : g_pdes) eq.initialize( m_coord, m_u, m_t, m_gid );

  // Compute initial time step size
  dt();

  // Output initial conditions to file (regardless of whether it was requested)
  if ( !g_inputdeck.get< tag::cmd, tag::benchmark >() ) writeFields( m_t );

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
Carrier::rhs( const tk::Fields& sol )
// *****************************************************************************
// Compute right-hand side of transport equations
//! \param[in] sol Solution vector at current stage
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

  // Compute right-hand side and query Dirichlet BCs for all equations
  tk::Fields r( m_gid.size(), g_inputdeck.get< tag::component >().nprop() );
  // Scale dt by 0.5 for first stage in 2-stage time stepping
  if (m_stage < 1) m_dt *= 0.5;
  for (const auto& eq : g_pdes) eq.rhs( m_t, m_dt, m_coord, m_inpoel, sol, r );
  // Query Dirichlet BCs and send to linear system merger
  bc();
  // Revert dt so that our copy is consistent with that of the host
  if (m_stage < 1) m_dt /= 0.5;
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
// *****************************************************************************
{
  tk::ExodusIIMeshReader
    er( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );

  auto nnode = er.readHeader();

  auto& x = m_coord[0];
  auto& y = m_coord[1];
  auto& z = m_coord[2];

  auto nn = m_lid.size();
  x.resize( nn );
  y.resize( nn );
  z.resize( nn );

  for (auto p : m_gid) {
    auto n = m_filenodes.find(p);
    if (n != end(m_filenodes) && n->second < nnode)
      er.readNode( n->second, tk::cref_find(m_lid,n->first), x, y, z );
  }
}

void
Carrier::addEdgeNodeCoords()
// *****************************************************************************
//  Add coordinates of mesh nodes newly generated to edge-mid points during
//  initial refinement
// *****************************************************************************
{
  if (m_edgenodes.empty()) return;

  auto& x = m_coord[0];
  auto& y = m_coord[1];
  auto& z = m_coord[2];
  Assert( x.size() == y.size() && x.size() == z.size(), "Size mismatch" );

  // Lambda to add coordinates for a single new node on an edge
  auto addnode = [ this, &x, &y, &z ]
                 ( const decltype(m_edgenodes)::value_type& e )
  {
    auto p = tk::cref_find( m_lid, e.first[0] );
    auto q = tk::cref_find( m_lid, e.first[1] );
    auto i = tk::cref_find( m_lid, e.second );
    Assert( p < x.size(), "Carrier chare " + std::to_string(thisIndex) +
                          " indexing out of bounds: " + std::to_string(p)
                          + " must be lower than " + std::to_string(x.size()) );
    Assert( q < x.size(), "Carrier chare " + std::to_string(thisIndex) +
                          " indexing out of bounds: " + std::to_string(q)
                          + " must be lower than " + std::to_string(x.size()) );
    Assert( i < x.size(), "Carrier chare " + std::to_string(thisIndex) +
                          " indexing out of bounds: " + std::to_string(i)
                          + " must be lower than " + std::to_string(x.size()) );
    x[i] = (x[p]+x[q])/2.0;
    y[i] = (y[p]+y[q])/2.0;
    z[i] = (z[p]+z[q])/2.0;
  };

  // add new nodes
  for (const auto& e : m_edgenodes) addnode( e );
}

void
Carrier::writeMesh()
// *****************************************************************************
// Output chare element blocks to file
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
//! \param[in] u Vector of fields to write to file
// *****************************************************************************
{
  int varid = 0;
  for (const auto& f : u) ew.writeNodeScalar( it, ++varid, f );
}

void
Carrier::writeMeta() const
// *****************************************************************************
// Output mesh-based fields metadata to file
// *****************************************************************************
{
  if (!g_inputdeck.get< tag::cmd, tag::benchmark >()) {

    // Create ExodusII writer
    tk::ExodusIIMeshWriter ew( m_outFilename, tk::ExoWriter::OPEN );

    // Collect nodal field output names from all PDEs
    std::vector< std::string > names;
    for (const auto& eq : g_pdes) {
      auto n = eq.fieldNames();
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
// *****************************************************************************
{
  // Only write if the last time is different than the current one
  if (std::abs(m_lastFieldWriteTime - time) <
      std::numeric_limits< tk::real >::epsilon() )
    return;

  // Save time stamp at which the last field write happened
  m_lastFieldWriteTime = time;

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
    auto o = eq.fieldOutput( time, m_V, m_coord, m_v, u );
    output.insert( end(output), begin(o), end(o) );
  }
  // Write node fields
  writeSolution( ew, m_itf, output );
}

void
Carrier::doWriteParticles()
// *****************************************************************************
// Output particles fields to file
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
// *****************************************************************************
{
  // Compute and sum antidiffusive element contributions to mesh nodes. Note
  // that the sums are complete on nodes that are not shared with other chares
  // and only partial sums on chare-boundary nodes.
  auto& dbc = m_linsysmerger.ckLocalBranch()->dirbc();
  m_fluxcorrector.aec( m_coord, m_inpoel, m_vol, dbc, m_gid, m_du, Un, m_p );
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

  m_fluxcorrector.lim( m_inpoel, m_p, (m_stage<1?m_ulf:m_ul), m_q, m_a );
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
  if (m_stage < 1) rhs( m_u ); else rhs( m_uf );
}

void
Carrier::out()
// *****************************************************************************
// Output mesh and particle fields
// *****************************************************************************
{
  // Optionally output field and particle data
  if ( m_stage == 1 &&
       !((m_it+1) % g_inputdeck.get< tag::interval, tag::field >()) &&
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
// *****************************************************************************
{
  // Optionally: collect analytical solutions and send both the latest
  // analytical and numerical solutions to LinSysMerger for computing and
  // outputing diagnostics
  if (m_stage == 1 && !(m_it % g_inputdeck.get< tag::interval, tag::diag >())) {

    // Collect analytical solutions (if available) from all PDEs. Note that
    // calling the polymorphic PDE::initialize() is assumed to evaluate the
    // analytical solution for a PDE. For those PDE problems that have
    // analytical solutions, this is the same as used for setting the initial
    // conditions, since if the analytical solution is available for a problem
    // it is (so far anyway) always initialized to its analytical solution and
    // that is done by calling PDE::initialize(). If the analytical solution is
    // a function of time, that is already incorporated in setting initial
    // conditions. For those PDEs where the analytical solution is not
    // available, initialize() returns the initial conditions (obviously), and
    // thus the "error", defined between this "analytical" and numerical
    // solution will be a measure of the "distance" between the initial
    // condition and the current numerical solution. This is not necessarily
    // useful, but simplies the logic because all PDEs can be treated as baing
    // able to compute an error based on some "analytical" solution, which
    // really the initial condition.
    for (const auto& eq : g_pdes)
      eq.initialize( m_coord, m_ul, m_t+m_dt, m_gid );

    // Prepare for computing diagnostics. Diagnostics are defined as the L2
    // norm of a quantity, computed in mesh nodes, A, as || A ||_2 = sqrt[
    // sum_i ( A_i )^2 V_i ], where the sum is taken over all mesh nodes and
    // V_i is the nodal volume. We send two sets of quantities to the host for
    // aggregation across the whole mesh: (1) the numerical solutions of all
    // components of all PDEs, and their error, defined as A_i = (a_i - n_i),
    // where a_i and n_i are the analytical and numerical solutions at node i,
    // respectively. We send these to LinSysMerger and the final aggregated
    // solution will end up in Transporter::diagnostics(). Note that we
    // weigh/multiply all data here by sqrt(V_i), so that the nodal volumes do
    // not have to be communicated separately. In LinSysMerger::diagnostics(),
    // where we collect all contributions from chares on a PE, all data is
    // squared. LinSysMerger::diagnostics() is where the sums are computed,
    // then the sums are summed accross the whole problem in
    // Transporter::diagnostics(), where the final square-root of the L2 norm,
    // defined above, is taken.
    for (std::size_t p=0; p<m_u.nunk(); ++p)
      for (ncomp_t c=0; c<g_inputdeck.get<tag::component>().nprop(); ++c) {
        m_ul(p,c,0) = m_ul(p,c,0);
        m_ulf(p,c,0) = m_u(p,c,0);
      }

    // Send both numerical and analytical solutions to linsysmerger
    m_linsysmerger.ckLocalBranch()->
      charediag( thisIndex, m_gid, m_ulf, m_ul, m_v );

  } else
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

bool
Carrier::correctBC()
// *****************************************************************************
//  Verify that the change in the solution at those nodes where Dirichlet
//  boundary conditions are set is exactly the amount the BCs prescribe
//! \return True if the solution is correct at Dirichlet boundary condition
//!   nodes
// *****************************************************************************
{
  auto& dirbc = m_linsysmerger.ckLocalBranch()->dirbc();

  if (dirbc.empty()) return true;

  for (const auto& s : m_side)
    for (auto i : s.second) {
      auto u = dirbc.find( m_gid[i] );
      if (u != end(dirbc)) {
        const auto& b = u->second;
        Assert( b.size() == m_u.nprop(), "Size mismatch" );
        for (std::size_t c=0; c<b.size(); ++c)
          if ( b[c].first &&
               std::abs( m_dul(i,c,0) + m_a(i,c,0) - b[c].second ) >
                 std::numeric_limits< tk::real >::epsilon() ) {
             return false;
          }
      }
  }

  return true;
}

void
Carrier::apply()
// *****************************************************************************
// Apply limited antidiffusive element contributions
// *****************************************************************************
{
  // Combine own and communicated contributions to A
  for (const auto& b : m_bid) {
    auto lid = tk::cref_find( m_lid, b.first );
    const auto& bac = m_ac[ b.second ];
    for (ncomp_t c=0; c<m_a.nprop(); ++c) m_a(lid,c,0) += bac[c];
  }

  // Verify that solution values do not change at Dirichlet BC nodes
  Assert( correctBC(), "Dirichlet boundary condition incorrect" );

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

  // Output final field data to file (regardless of whether it was requested)
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  if ( (m_stage == 1) &&
       (std::fabs(m_t+m_dt-term) < eps || (m_it+1) >= nstep) &&
       (!g_inputdeck.get< tag::cmd, tag::benchmark >()) )
    writeFields( m_t+m_dt );

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
