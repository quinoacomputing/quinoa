//******************************************************************************
/*!
  \file      src/Inciter/Performer.C
  \author    J. Bakosi
  \date      Fri 20 Nov 2015 06:19:01 AM MST
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

Performer::Performer( int id,
                      ConductorProxy& conductor,
                      LinSysMergerProxy& linsysmerger,
                      SpawnerProxy& spawner,
                      const std::vector< std::size_t >& element ) :
  m_id( static_cast< std::size_t >( id ) ),
  m_it( 0 ),
  m_itf( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_stage( 0 ),
  m_conductor( conductor ),
  m_linsysmerger( linsysmerger ),
  m_spanwer( spawner ),
  m_elem( element )
//******************************************************************************
//  Constructor
//! \param[in] id Charm++ global array id
//! \param[in] host Host proxy
//! \param[in] lsm Linear system merger (LinSysMerger) proxy
//! \param[in] element Vector global mesh element IDs owned
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
//std::cout << "setup: " << CkMyPe() << ", " << m_id << ", " << m_elem.size() << '\n';
  // Initialize local->global, global->local node ids, element connectivity
  initIds( m_elem );
  // Read coordinates of owned and received mesh nodes
  initCoords();
  // Output chare mesh to file
  writeMesh();
  // Output mesh-based fields metadata to file
  writeMeta();
}

void
Performer::init( tk::real dt )
//******************************************************************************
// Initialize linear system
//! \author J. Bakosi
//******************************************************************************
{
//std::cout << "init: " << CkMyPe() << ", " << m_id << '\n';
  // Set initial conditions
  ic();

  // If the desired max number of time steps is zero or the desired max time is
  // less than the initial time step size, finish righ away.
  if ( g_inputdeck.get< tag::discr, tag::nstep >() == 0 ||
       g_inputdeck.get< tag::discr, tag::t0 >() >
         g_inputdeck.get< tag::discr, tag::term >() ||
       g_inputdeck.get< tag::discr, tag::term >() < dt ) {

    // Send time stamps to the host
    m_conductor.arrTimestamp( m_timestamp );

    // Tell the Charm++ runtime system to call back to Conductor::finish(). The
    // reduction is done via creating a callback that invokes the typed
    // reduction client, where m_conductor is the proxy on which the reduction
    // target method, finish(), is called upon completion of the reduction.
    contribute(
      CkCallback( CkReductionTarget( Conductor, finish ), m_conductor ) );

  } else {

    // Compute left-hand side of PDE
    lhs();
    // Advance PDE in time (start at stage 0)
    advance( 0, dt, m_it, m_t );
    // Send time stamps to the host
    m_conductor.arrTimestamp( m_timestamp );

  }
}

void
Performer::ic()
//******************************************************************************
// Set initial conditions
//! \author J. Bakosi
//******************************************************************************
{
  for (std::size_t i=0; i<m_gid.size(); ++i)
    m_u[ m_gid[i] ] = ansol_shear( i, m_t );
    //m_u[ m_gid[i] ] = ansol_gauss( i );

  // Output initial conditions to file (it = 1, time = 0.0)
  writeFields( m_t );

  m_linsysmerger.ckLocalBranch()->charesol( static_cast<int>(m_id), m_u );
}

void
Performer::lhs()
//******************************************************************************
// Compute left-hand side of PDE
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  // Sparse matrix: global mesh point row and column ids, and nonzero value
  std::map< std::size_t, std::map< std::size_t, tk::real > > lhs;

  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    const auto a = m_inpoel[e*4+0];
    const auto b = m_inpoel[e*4+1];
    const auto c = m_inpoel[e*4+2];
    const auto d = m_inpoel[e*4+3];
    std::array< tk::real, 3 > ba{{ x[b]-x[a], y[b]-y[a], z[b]-z[a] }},
                              ca{{ x[c]-x[a], y[c]-y[a], z[c]-z[a] }},
                              da{{ x[d]-x[a], y[d]-y[a], z[d]-z[a] }};
    const auto J = tk::triple( ba, ca, da ) / 120.0;

    const auto A = m_gid[a];
    const auto B = m_gid[b];
    const auto C = m_gid[c];
    const auto D = m_gid[d];
    auto& lA = lhs[A];
    lA[A] += 2.0*J;
    lA[B] += J;
    lA[C] += J;
    lA[D] += J;
    auto& lB = lhs[B];
    lB[A] += J;
    lB[B] += 2.0*J;
    lB[C] += J;
    lB[D] += J;
    auto& lC = lhs[C];
    lC[A] += J;
    lC[B] += J;
    lC[C] += 2.0*J;
    lC[D] += J;
    auto& lD = lhs[D];
    lD[A] += J;
    lD[B] += J;
    lD[C] += J;
    lD[D] += 2.0*J;
  }

  m_timestamp.emplace_back( "Compute left-hand side matrix", t.dsec() );

  m_linsysmerger.ckLocalBranch()->charelhs( static_cast<int>(m_id), lhs );
}

void
Performer::rhs( tk::real mult,
                tk::real dt,
                const std::map< std::size_t, tk::real >& sol,
                std::map< std::size_t, tk::real >& newrhs )
//******************************************************************************
// Compute right-hand side of PDE
//! \param[in] mult Multiplier differentiating the different stages in
//!    multi-stage time stepping
//! \param[in] dt Size of time step
//! \param[in] sol Solution vector at current stage
//! \param[inout] newrhs Right-hand side vector computed
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  const tk::real U0 = 0.5;
  const tk::real LAMBDA = 5.0e-4;
  const tk::real D = 10.0;

  newrhs.clear();

  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    const auto a = m_inpoel[e*4+0];
    const auto b = m_inpoel[e*4+1];
    const auto c = m_inpoel[e*4+2];
    const auto d = m_inpoel[e*4+3];

    // compute element Jacobi determinant
    std::array< tk::real, 3 > ba{{ x[b]-x[a], y[b]-y[a], z[b]-z[a] }},
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
    std::array< tk::real, 4 > u {{ tk::cref_find( m_u, m_gid[a] ),
                                   tk::cref_find( m_u, m_gid[b] ),
                                   tk::cref_find( m_u, m_gid[c] ),
                                   tk::cref_find( m_u, m_gid[d] ) }};
    // solution at nodes at time n (at stage 0) and n+1/2 (at stage 1)
    std::array< tk::real, 4 > s {{ tk::cref_find( sol, m_gid[a] ),
                                   tk::cref_find( sol, m_gid[b] ),
                                   tk::cref_find( sol, m_gid[c] ),
                                   tk::cref_find( sol, m_gid[d] ) }};
    // pointers to rhs at nodes
    std::array< tk::real*, 4 > r {{ &newrhs[ m_gid[a] ],
                                    &newrhs[ m_gid[b] ],
                                    &newrhs[ m_gid[c] ],
                                    &newrhs[ m_gid[d] ] }};

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

  m_timestamp.emplace_back( "Compute right-hand side vector", t.dsec() );

//std::cout << "sendrhs: " << CkMyPe() << ", " << m_id << '\n';
  m_linsysmerger.ckLocalBranch()->charerhs( static_cast<int>(m_id), newrhs );
}

void
Performer::initIds( const std::vector< std::size_t >& gelem )
//******************************************************************************
//! Initialize local->global, global->local node ids, element connectivity
//! \param[in] gelem Set of unique owned global element ids
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  tk::ExodusIIMeshReader
    er( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );

  // Read element block IDs from ExodusII file
  er.readElemBlockIDs();

  std::vector< std::size_t > gtetinpoel;

  // Read global element connectivity of owned tetrahedron elements
  for (auto e : gelem) er.readElement( e, tk::ExoElemType::TET, gtetinpoel );

  m_timestamp.emplace_back( "Read mesh element connectivity from file",
                            t.dsec() );

  tk::Timer t2;

  // Generate connectivity graph storing local node ids
  std::tie( m_inpoel, m_gid ) = tk::global2local( gtetinpoel );

  // Send off number of columns per row to linear system merger
  m_linsysmerger.ckLocalBranch()->charerow( thisProxy,
                                            static_cast<int>(m_id),
                                            m_gid );

  m_timestamp.emplace_back( "Initialize mesh point ids, element connectivity",
                            t2.dsec() );
}

void
Performer::initCoords()
//******************************************************************************
//  Read coordinates of mesh nodes from file
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  tk::ExodusIIMeshReader
    er( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );

  auto& x = m_coord[0];
  auto& y = m_coord[1];
  auto& z = m_coord[2];
  for (auto p : m_gid) er.readNode( p, x, y, z );

  m_timestamp.emplace_back( "Read mesh point coordinates from file", t.dsec() );
}

void
Performer::writeMesh()
//******************************************************************************
// Output chare mesh to file
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  // Create mesh object initializing element connectivity and point coords
  tk::UnsMesh mesh( m_inpoel, m_coord );

  // Create ExodusII writer
  tk::ExodusIIMeshWriter
    ew( g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "_" +
          std::to_string( m_id ),
        tk::ExoWriter::CREATE );

  // Write chare mesh
  ew.writeMesh( mesh );

  m_timestamp.emplace_back( "Write chare mesh to file", t.dsec() );
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
                          const std::map< std::size_t, tk::real >& u ) const
//******************************************************************************
// Output solution to file
//! \param[in] ew ExodusII mesh-based writer object
//! \param[in] it Iteration count
//! \param[in] varid Exodus variable ID
//! \param[in] u Field to write
//! \author J. Bakosi
//******************************************************************************
{
  std::vector< tk::real > sol;
  for (const auto& p : u) sol.push_back( p.second );
  ew.writeNodeScalar( it, varid, sol );
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
  tk::Timer t;

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
  m_un.clear();
  for (std::size_t i=0; i<m_gid.size(); ++i)
    m_un[ m_gid[i] ] = ansol_shear( i, time );
  writeSolution( ew, m_itf, 2, m_un );

  m_timestamp.emplace_back( "Write mesh-based fields to file", t.dsec() );
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
Performer::updateSolution( const std::map< std::size_t, tk::real >& sol )
//******************************************************************************
// Update solution vector
//! \author J. Bakosi
//******************************************************************************
{
  // Receive update of solution vector
  for (const auto& s : sol) m_ur[ s.first ] = s.second;

//std::cout << "UpdateRecvd " << CkMyPe() << ", " << m_id << ", ursize now: " << m_ur.size() << ", when complete: " << m_gid.size() << '\n';

  // If all contributions we own have been received, continue by updating a
  // different solution vector depending on time step stage
  if (m_ur.size() == m_gid.size()) {

//std::cout << "allUpdateRecvd " << CkMyPe() << ", " << m_id << '\n';

    if (m_stage < 1) {

      m_uf = std::move( m_ur );

    } else {

      m_u = std::move( m_ur );
      if (!((m_it+1) % g_inputdeck.get< tag::interval, tag::field >()))
        writeFields( m_t + g_inputdeck.get< tag::discr, tag::dt >() );

    }

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
