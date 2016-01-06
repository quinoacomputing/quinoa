//******************************************************************************
/*!
  \file      src/Inciter/Performer.C
  \author    J. Bakosi
  \date      Wed 06 Jan 2016 09:46:57 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Performer advances a PDE
  \details   Performer advances a PDE. There are a potentially
    large number of Performer Charm++ chares created by Conductor. Each
    performer gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances a PDE in time.
*/
//******************************************************************************

#include <string>

#include "Performer.h"
#include "Vector.h"
#include "Reader.h"
#include "ContainerUtil.h"
#include "UnsMesh.h"
#include "Reorder.h"
#include "ExodusIIMeshReader.h"
#include "ExodusIIMeshWriter.h"
#include "LinSysMerger.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "DerivedData.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::Performer;

Performer::Performer(
  int id,
  ConductorProxy& conductor,
  LinSysMergerProxy& linsysmerger,
  SpawnerProxy& spawner,
  const std::vector< std::size_t >& conn,
  const std::unordered_map< std::size_t, std::size_t >& cid )
:
  m_id( id ),
  m_it( 0 ),
  m_itf( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_stage( 0 ),
  m_nsol( 0 ),
  m_conductor( conductor ),
  m_linsysmerger( linsysmerger ),
  m_spanwer( spawner ),
  m_conn( conn ),
  m_cid( cid )
//******************************************************************************
//  Constructor
//! \param[in] id Charm++ global array id
//! \param[in] host Host proxy
//! \param[in] lsm Linear system merger (LinSysMerger) proxy
//! \param[in] conn Vector of mesh element connectivity owned (global IDs)
//! \param[in] cid Map associating old node IDs (as in file) to new node IDs (as
//!   in producing contiguous-row-id linear system contributions)
//! \details Since a Performer chare array is created separately on each PE, the
//!   chare array index, thisIndex, is a local index. The global index, unknown
//!   to Charm, is unique across all PEs, stored in m_id.
//! \author J. Bakosi
//******************************************************************************
{
  // Register ourselves with the linear system merger
  m_linsysmerger.ckLocalBranch()->checkin();
}

void
Performer::setup()
//******************************************************************************
// Initialize mesh IDs, element connectivity, coordinates
//! \author J. Bakosi
//******************************************************************************
{
  // Initialize local->global, global->local node ids, element connectivity
  setupIds( m_conn );
  // Read coordinates of owned and received mesh nodes
  readCoords();
  // Output chare mesh to file
  writeMesh();
  // Output mesh-based fields metadata to file
  writeMeta();
}

void
Performer::setupIds( const std::vector< std::size_t >& gelem )
//******************************************************************************
//! Initialize local->global, global->local node ids, element connectivity
//! \param[in] gelem Set of unique owned global ids
//! \author J. Bakosi
//******************************************************************************
{
  // Generate connectivity graph storing local node ids
  std::tie( m_inpoel, m_gid ) = tk::global2local( gelem );

  // Send off number of columns per row to linear system merger
  m_linsysmerger.ckLocalBranch()->charerow( thisProxy, m_id, thisIndex, m_gid );
}

void
Performer::init( tk::real dt )
//******************************************************************************
// Initialize linear system
//! \author J. Bakosi
//******************************************************************************
{
  // Set initial conditions
  ic();

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
Performer::ic()
//******************************************************************************
// Set initial conditions
//! \author J. Bakosi
//******************************************************************************
{
  m_u.resize( m_gid.size() );
  m_uf.resize( m_gid.size() );
  m_un.resize( m_gid.size() );

  for (std::size_t i=0; i<m_gid.size(); ++i)
    m_u[i] = ansol_shear( i, m_t );
    //m_u[i] = ansol_gauss( i );

  // Output initial conditions to file (it = 1, time = 0.0)
  writeFields( m_t );

  m_linsysmerger.ckLocalBranch()->charesol( m_id, m_gid, m_u );
}

void
Performer::lhs()
//******************************************************************************
// Compute left-hand side of PDE
//! \author J. Bakosi
//******************************************************************************
{
  // Generate points surrounding points
  const auto psup = tk::genPsup( m_inpoel, 4, tk::genEsup(m_inpoel,4) );

  Assert( psup.second.size()-1 == m_gid.size(),
          "Number of mesh points and number of global IDs unequal" );

  // Sparse matrix storing the nonzero matrix values at rows and columns given
  // by psup. The format is similar to compressed row storage, but the diagonal
  // and off-diagonal data are stored in separate vectors. For the off-diagonal
  // data the local row and column indices, at which values are nonzero, are
  // stored by psup (psup1 and _psup2, where psup2 holds the indices at which
  // psup1 holds the point ids surrounding points, see also tk::genPsup()). Note
  // that the number of mesh points (our chunk) npoin = psup.second.size()-1.
  std::vector< tk::real > lhsd( psup.second.size()-1, 0.0 );
  std::vector< tk::real > lhso( psup.first.size(), 0.0 );

  // Lambda to compute the sparse matrix vector index for row and column
  // indices. Used only for off-diagonal entries.
  auto spidx = [ &psup ]( std::size_t r, std::size_t c ) -> std::size_t {
    Assert( r != c, "Only for computing the off-diagonal indices" );
    for (auto i=psup.second[r]+1; i<=psup.second[r+1]; ++i)
      if (c == psup.first[i]) return i;
    Throw( "Cannot find row, column: " + std::to_string(r) + ',' +
           std::to_string(c) + " in sparse matrix" );
  };

  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    const auto a = m_inpoel[e*4+0];
    const auto b = m_inpoel[e*4+1];
    const auto c = m_inpoel[e*4+2];
    const auto d = m_inpoel[e*4+3];
    std::array< tk::real, 3 > ba{{ x[b]-x[a], y[b]-y[a], z[b]-z[a] }},
                              ca{{ x[c]-x[a], y[c]-y[a], z[c]-z[a] }},
                              da{{ x[d]-x[a], y[d]-y[a], z[d]-z[a] }};
    const auto J = tk::triple( ba, ca, da ) / 120.0;

    lhsd[ a ] += 2.0*J;
    lhsd[ b ] += 2.0*J;
    lhsd[ c ] += 2.0*J;
    lhsd[ d ] += 2.0*J;

    lhso[ spidx(a,b) ] += J;
    lhso[ spidx(a,c) ] += J;
    lhso[ spidx(a,d) ] += J;

    lhso[ spidx(b,a) ] += J;
    lhso[ spidx(b,c) ] += J;
    lhso[ spidx(b,d) ] += J;

    lhso[ spidx(c,a) ] += J;
    lhso[ spidx(c,b) ] += J;
    lhso[ spidx(c,d) ] += J;

    lhso[ spidx(d,a) ] += J;
    lhso[ spidx(d,b) ] += J;
    lhso[ spidx(d,c) ] += J;
  }

  m_linsysmerger.ckLocalBranch()->charelhs( m_id, m_gid, psup, lhsd, lhso );
}

void
Performer::rhs( tk::real mult,
                tk::real dt,
                const std::vector< tk::real >& sol,
                std::vector< tk::real >& rhs )
//******************************************************************************
// Compute right-hand side of PDE
//! \param[in] mult Multiplier differentiating the different stages in
//!    multi-stage time stepping
//! \param[in] dt Size of time step
//! \param[in] sol Solution vector at current stage
//! \param[inout] rhs Right-hand side vector computed
//! \author J. Bakosi
//******************************************************************************
{
  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  const tk::real U0 = 0.5;
  const tk::real LAMBDA = 5.0e-4;
  const tk::real D = 10.0;

  std::fill( begin(rhs), end(rhs), 0.0 );

  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    const auto a = m_inpoel[e*4+0];
    const auto b = m_inpoel[e*4+1];
    const auto c = m_inpoel[e*4+2];
    const auto d = m_inpoel[e*4+3];

    // compute element Jacobi determinant
    const std::array< tk::real, 3 > ba{{ x[b]-x[a], y[b]-y[a], z[b]-z[a] }},
                                    ca{{ x[c]-x[a], y[c]-y[a], z[c]-z[a] }},
                                    da{{ x[d]-x[a], y[d]-y[a], z[d]-z[a] }};
    const auto J = tk::triple( ba, ca, da );

    // construct tetrahedron element-level matrices

    // consistent mass
    std::array< std::array< tk::real, 4 >, 4 > mass;  // nnode*nnode [4][4]
    // diagonal
    mass[0][0] = mass[1][1] = mass[2][2] = mass[3][3] = J/60.0;
    // off-diagonal
    mass[0][1] = mass[0][2] = mass[0][3] =
    mass[1][0] = mass[1][2] = mass[1][3] =
    mass[2][0] = mass[2][1] = mass[2][3] =
    mass[3][0] = mass[3][1] = mass[3][2] = J/120.0;

    // prescribed shear velocity
    std::array< std::array< tk::real, 4 >, 3 > vel;  // ndim*nnode [3][4]
    vel[0][0] = U0 + LAMBDA*y[a];  vel[1][0] = 0.0;  vel[2][0] = 0.0;
    vel[0][1] = U0 + LAMBDA*y[b];  vel[1][1] = 0.0;  vel[2][1] = 0.0;
    vel[0][2] = U0 + LAMBDA*y[c];  vel[1][2] = 0.0;  vel[2][2] = 0.0;
    vel[0][3] = U0 + LAMBDA*y[d];  vel[1][3] = 0.0;  vel[2][3] = 0.0;

    // shape function derivatives
    std::array< std::array< tk::real, 3 >, 4 > grad;  // nnode*ndim [4][3]
    grad[1] = tk::crossdiv( ca, da, J );
    grad[2] = tk::crossdiv( da, ba, J );
    grad[3] = tk::crossdiv( ba, ca, J );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

    // solution at nodes at time n
    const std::array< tk::real, 4 > u{{ m_u[a], m_u[b], m_u[c], m_u[d] }};
    // solution at nodes at time n (at stage 0) and n+1/2 (at stage 1)
    const std::array< tk::real, 4 > s{{ sol[a], sol[b], sol[c], sol[d] }};
    // pointers to rhs at nodes
    std::array< tk::real*, 4 > r{{ &rhs[a], &rhs[b], &rhs[c], &rhs[d] }};

    // add mass contribution to rhs
    for (std::size_t i=0; i<4; ++i)
      for (std::size_t j=0; j<4; ++j)
        *r[i] += mass[i][j] * u[j];

    // add advection contribution to rhs
    for (std::size_t i=0; i<4; ++i)
      for (std::size_t j=0; j<4; ++j)
        for (std::size_t k=0; k<3; ++k)
          for (std::size_t l=0; l<4; ++l)
            *r[i] -= mult * dt * mass[i][j] * vel[k][j] * grad[l][k] * s[l];

    // add diffusion contribution to rhs
    for (std::size_t i=0; i<4; ++i)
      for (std::size_t j=0; j<4; ++j)
        for (std::size_t k=0; k<3; ++k)
          *r[i] -= mult * dt * D * J/6.0 * grad[i][k] * grad[j][k] * s[j];
  }

  m_linsysmerger.ckLocalBranch()->charerhs( m_id, m_gid, rhs );
}

void
Performer::readCoords()
//******************************************************************************
//  Read coordinates of mesh nodes from file
//! \author J. Bakosi
//******************************************************************************
{
  tk::ExodusIIMeshReader
    er( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );

  auto& x = m_coord[0];
  auto& y = m_coord[1];
  auto& z = m_coord[2];
  for (auto p : m_gid) er.readNode( tk::val_find(m_cid,p), x, y, z );
}

void
Performer::writeMesh()
//******************************************************************************
// Output chare mesh to file
//! \author J. Bakosi
//******************************************************************************
{
  // Create mesh object initializing element connectivity and point coords
  tk::UnsMesh mesh( m_inpoel, m_coord );

  // Create ExodusII writer
  tk::ExodusIIMeshWriter
    ew( g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "_" +
          std::to_string( m_id ),
        tk::ExoWriter::CREATE );

  // Write chare mesh
  ew.writeMesh( mesh );
}

void
Performer::writeChareId( const tk::ExodusIIMeshWriter& ew,
                         uint64_t it ) const
//******************************************************************************
// Output chare id field to file
//! \param[in] ew ExodusII mesh-based writer object
//! \param[in] it Iteration count
//! \author J. Bakosi
//******************************************************************************
{
  // Write elem chare id field to mesh
  std::vector< tk::real > chid( m_inpoel.size()/4, static_cast<tk::real>(m_id) );
  ew.writeElemScalar( it, 1, chid );
}

void
Performer::writeSolution( const tk::ExodusIIMeshWriter& ew,
                          uint64_t it,
                          int varid,
                          const std::vector< tk::real >& u ) const
//******************************************************************************
// Output solution to file
//! \param[in] ew ExodusII mesh-based writer object
//! \param[in] it Iteration count
//! \param[in] varid Exodus variable ID
//! \param[in] u Field to write
//! \author J. Bakosi
//******************************************************************************
{
  ew.writeNodeScalar( it, varid, u );
}

void
Performer::writeMeta() const
//******************************************************************************
// Output mesh-based fields metadata to file
//! \author J. Bakosi
//******************************************************************************
{
  // Create ExodusII writer
  tk::ExodusIIMeshWriter
    ew( g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "_" +
          std::to_string( m_id ),
        tk::ExoWriter::OPEN );

  ew.writeElemVarNames( { "Chare Id" } );
  ew.writeNodeVarNames( { "NumSol", "AnSol" } );
}

void
Performer::writeFields( tk::real time )
//******************************************************************************
// Output mesh-based fields to file
//! \param[in] time Physical time
//! \author J. Bakosi
//******************************************************************************
{
  // Increase field output iteration count
  ++m_itf;

  // Create ExodusII writer
  tk::ExodusIIMeshWriter
    ew( g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "_" +
          std::to_string( m_id ),
        tk::ExoWriter::OPEN );

  // Write time stamp
  ew.writeTimeStamp( m_itf, time );

  // Write mesh-based fields
  writeChareId( ew, m_itf );
  writeSolution( ew, m_itf, 1, m_u );

  // Analytical solution for this time
  for (std::size_t i=0; i<m_gid.size(); ++i)
    m_un[i] = ansol_shear( i, time );
  writeSolution( ew, m_itf, 2, m_un );
}

void
Performer::advance( uint8_t stage, tk::real dt, uint64_t it, tk::real t )
//******************************************************************************
// Advance equations to next stage in multi-stage time stepping
//! \param[in] stage Stage in multi-stage time stepping
//! \param[in] dt Size of time step
//! \param[in] it Iteration count
//! \param[in] t Physical time
//! \author J. Bakosi
//******************************************************************************
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
Performer::updateSolution( const std::vector< std::size_t >& gid,
                           const std::vector< tk::real >& sol )
//******************************************************************************
// Update solution vector
//! \param[in] gid Global row indices of the vector updated
//! \param[in] sol Portion of the unknown/solution vector updated
//! \author J. Bakosi
//******************************************************************************
{
  Assert( gid.size() == sol.size(),
          "Size of solution and row ID vectors must equal" );

  // Receive update of solution vector
  for (std::size_t i=0; i<gid.size(); ++i) {
    auto it = std::find( cbegin(m_gid), cend(m_gid), gid[i] );
    if (it != cend(m_gid)) {
      auto p = static_cast< std::size_t >( std::distance( cbegin(m_gid), it ) );
      m_un[ p ] = sol[i];
    } else Throw( "Cannot find global node ID for solution received" );
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

  }
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "performer.def.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
