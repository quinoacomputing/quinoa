// *****************************************************************************
/*!
  \file      src/Inciter/Carrier.C
  \author    J. Bakosi
  \date      Thu 13 Oct 2016 10:07:20 AM MDT
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
CkReduction::reducerType FieldsMerger;

} // inciter::

using inciter::Carrier;

Carrier::Carrier( const TransporterProxy& transporter,
                  const LinSysMergerProxy& lsm,
                  const ParticleWriterProxy& pw,
                  const std::vector< std::size_t >& conn,
                  const std::unordered_map< std::size_t, std::size_t >& cid,
                  int ncarr ) :
  __dep(),
  m_it( 0 ),
  m_itf( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_dt( g_inputdeck.get< tag::discr, tag::dt >() ),
  m_stage( 0 ),
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
  m_msum(),
  m_vol( m_gid.size(), 0.0 ),
  m_bid(),
  m_pc(),
  m_qc(),
  m_ac(),
  m_tracker( 0, m_inpoel ) // 0 = no particles
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
}

void
Carrier::comm()
// *****************************************************************************
// Starts collecting global mesh node IDs bordering the mesh chunk held by
// fellow Carrier chares associated to their chare IDs
//! \author J. Bakosi
// *****************************************************************************
{
  // Read coordinates of nodes of the mesh chunk we operate on
  readCoords();

  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

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

  // Collect fellow chare IDs (as well as the global mesh node IDs and the
  // volume associated to them) to/from our mesh neighbors. The goal here is to
  // end up with a list of chare IDs (as the keys in m_msum) on each chare that
  // the chare shares at least a single mesh node ID with, i.e., its mesh chunk
  // borders another chare's mesh chunk. Note that this is not necessarily the
  // most efficient way of doing this, since we do an all-to-all reduction from
  // and to all Carrier chares, sending all the global node indices of the mesh
  // chunk each Carrier chare holds. This contribute() is the start of this, and
  // msum() is the final step. A better way of doing this is to only send the
  // chare-boundary node indices. That information might already be available in
  // Partitioner (the host chare group we are created by), or if not directly,
  // it might be possible to be computed from the output of the partitioner.
  // This data structure (m_msum) is used in various ways: (1) the list of
  // neighboring chare IDs are used for searching for particles traveling across
  // the mesh that left our mesh chunk, (2) the list of neighboring mesh node
  // IDs (that lie on the chare-boundary) and their associated volumes are used
  // by the FCT algorithm.
  std::vector< std::pair< int, comm_t > >
    meshnodes{ { thisIndex, std::make_tuple( m_gid, m_vol ) } };
  auto stream = serialize( meshnodes );
  CkCallback cb( CkIndex_Carrier::msum(nullptr), thisProxy );
  contribute( stream.first, stream.second.get(), FieldsMerger, cb );
}

void
Carrier::msum( CkReductionMsg* msg )
// *****************************************************************************
// Reduction target finishing collecting global mesh node IDs bordering the
// mesh chunk held by fellow Carrier chares associated to their chare IDs
//! \param[in] msg Serialized aggregated mesh node IDs categorized by chares
//! \author J. Bakosi
// *****************************************************************************
{
  std::vector< std::pair< int, comm_t > > meshnodes;
  PUP::fromMem creator( msg->getData() );
  creator | meshnodes;
  delete msg;

  // Store nodes that our mesh chunk neighbors associated to chare IDs
  for (const auto& c : meshnodes)
    if (c.first != thisIndex) {
      if (m_ncarr == 1) Throw( "Should not be called with a single chare" );
      const auto& gid = std::get< 0 >( c.second );
      const auto& vol = std::get< 1 >( c.second );
      for (std::size_t i=0; i<m_gid.size(); ++i)
        for (std::size_t n=0; n<gid.size(); ++n)
          if (m_gid[i] == gid[n]) {
            m_msum[ c.first ].push_back( m_gid[i] );
            m_vol[ tk::cref_find(m_lid,m_gid[i]) ] += vol[n];
            n = gid.size();
          }
    }

//   std::cout << thisIndex << "> ";
//   for (const auto& n : m_msum) {
//     std::cout << n.first << ": ";
//     for (auto i : n.second) std::cout << i << ' ';
//     std::cout << ' ';
//   }
//   std::cout << '\n';

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
  // Generate particles
  m_tracker.genpar( m_coord, m_inpoel, m_ncarr, thisIndex );
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
Carrier::init()
// *****************************************************************************
// Initialize linear system
//! \author J. Bakosi
// *****************************************************************************
{
  // zero initial solution vector
  m_du.fill( 0.0 );

  // Send off initial guess for assembly
  m_linsysmerger.ckLocalBranch()->charesol( thisIndex, m_gid, m_du );

  // Set initial conditions for all PDEs
  for (const auto& eq : g_pdes) eq.initialize( m_coord, m_u, m_t );

  // Compute initial time step size
  dt();

  // Output initial conditions to file (time = 0.0)
  writeFields( 0.0 );

  // Call back to Transporter::initcomplete(), signaling that the initialization
  // is complete and we are now starting time stepping
  contribute(
      CkCallback(CkReductionTarget(Transporter,initcomplete), m_transporter));

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

  if (std::abs(const_dt - def_const_dt) > eps) {        // use constant dt

    mindt = const_dt;

  } else if (m_stage == 1) {      // compute dt based on CFL only in final stage

    const auto& x = m_coord[0];
    const auto& y = m_coord[1];
    const auto& z = m_coord[2];

    for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
      // Find largest of element nodal velocities
      auto ev = evel( e );
      const auto& evx = ev[0], evy = ev[1], evz = ev[2];
      const std::array< tk::real, 4 > vel = {{
        std::sqrt(evx[0]*evx[0] + evy[0]*evy[0] + evz[0]*evz[0]),
        std::sqrt(evx[1]*evx[1] + evy[1]*evy[1] + evz[1]*evz[1]),
        std::sqrt(evx[2]*evx[2] + evy[2]*evy[2] + evz[2]*evz[2]),
        std::sqrt(evx[3]*evx[3] + evy[3]*evy[3] + evz[3]*evz[3]) }};
      auto maxvel = *std::max_element( begin(vel), end(vel) );

      if (maxvel < std::numeric_limits<tk::real>::epsilon()) maxvel = 1.0;

      // Compute smallest edge length
      const std::array< std::size_t, 4 > N{{ m_inpoel[e*4+0], m_inpoel[e*4+1],
                                             m_inpoel[e*4+2], m_inpoel[e*4+3] }};
      const std::array< std::array< tk::real, 3 >, 6 > evec = {{
        {{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
        {{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
        {{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }},
        {{ x[N[2]]-x[N[1]], y[N[2]]-y[N[1]], z[N[2]]-z[N[1]] }},
        {{ x[N[3]]-x[N[1]], y[N[3]]-y[N[1]], z[N[3]]-z[N[1]] }},
        {{ x[N[3]]-x[N[2]], y[N[3]]-y[N[2]], z[N[3]]-z[N[2]] }} }};
      const std::array< tk::real, 6 > edge = {{
        std::sqrt( tk::dot( evec[0], evec[0] ) ),
        std::sqrt( tk::dot( evec[1], evec[1] ) ),
        std::sqrt( tk::dot( evec[2], evec[2] ) ),
        std::sqrt( tk::dot( evec[3], evec[3] ) ),
        std::sqrt( tk::dot( evec[4], evec[4] ) ),
        std::sqrt( tk::dot( evec[5], evec[5] ) ) }};
      auto minedge = *std::min_element( begin(edge), end(edge) );

      // Find smallest dt
      tk::real dt = g_inputdeck.get< tag::discr, tag::cfl >() * minedge / maxvel;
      if (dt < mindt) mindt = dt;
    }

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
  m_fluxcorrector.aec( m_coord, m_inpoel, m_vol, m_du, Un, m_p );
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
  // Update local copy of time step stage, physical time and time step size, and
  // iteration count
  m_stage = stage;
  m_t = t;
  m_dt = dt;
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
  if (m_stage == 1 &&
      !(m_it % g_inputdeck.get< tag::interval, tag::field >()))
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

std::array< std::array< tk::real, 4 >, 3 >
Carrier::evel( std::size_t e )
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

  // Advance particles
  m_tracker.track( m_transporter, thisProxy, m_coord, m_inpoel, m_msum,
                   thisIndex, this, m_stage, m_dt );

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

#include "NoWarning/carrier.def.h"
