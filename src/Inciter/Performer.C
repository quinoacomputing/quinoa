//******************************************************************************
/*!
  \file      src/Inciter/Performer.C
  \author    J. Bakosi
  \date      Thu 20 Aug 2015 11:39:03 AM MDT
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
#include "ContainerUtil.h"
#include "UnsMesh.h"
#include "ExodusIIMeshReader.h"
#include "ExodusIIMeshWriter.h"
#include "LinSysMerger.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< std::size_t > g_tetinpoel;
extern std::vector< std::size_t > g_meshfilemap;
extern std::vector< std::vector< std::size_t > > g_point;
extern std::vector< std::vector< std::size_t > > g_element;
extern std::vector< std::map< std::size_t, std::vector<std::size_t> > > g_ecomm;

} // inciter::

using inciter::Performer;

Performer::Performer( CProxy_Conductor& hostproxy, LinSysMergerProxy& lsmproxy )
: m_id( static_cast< std::size_t >( thisIndex ) ),
  m_it( 0 ),
  m_t( 0.0 ),
  m_hostproxy( hostproxy ),
  m_lsmproxy( lsmproxy ),
  m_point( g_point[ m_id ] )
//******************************************************************************
// Constructor
//! \param[in] hostproxy Host proxy
//! \param[in] lsmproxy Linear system merger (LinSysMerger) proxy
//! \author J. Bakosi
//******************************************************************************
{
  // Register ourselves with the linear system merger
  m_lsmproxy.ckLocalBranch()->checkin();
  // Tell the Charm++ runtime system to call back to Conductor::registered()
  // once all Performer chares have registered themselves, i.e., checked in,
  // with their local branch of the linear system merger group, LinSysMerger.
  // The reduction is done via creating a callback that invokes the typed
  // reduction client, where m_hostproxy is the proxy on which the reduction
  // target method, registered(), is called upon completion of the reduction.
  contribute(
    CkCallback( CkReductionTarget( Conductor, registered ), m_hostproxy ) );
}

void
Performer::setup()
//******************************************************************************
// Setup
//! \author J. Bakosi
//******************************************************************************
{
  // Initialize import maps
  initImports();
  // Initialize local->global, global->local node ids, element connectivity
  initIds( g_element[ m_id ] );
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
  // Set initial conditions
  ic();
  // Compute left-hand side of PDE
  lhs();
  // Compute righ-hand side of PDE
  rhs( 0.5, dt, m_u );
  // Send some time stamps to the host
  m_hostproxy.arrTimestamp( m_timestamp );
}

void
Performer::initImports()
//******************************************************************************
// Initialize import maps by inverting export maps
//! \author J. Bakosi
//******************************************************************************
{
  std::size_t h = 0;
  for (const auto& m : g_ecomm) {
    for (const auto& x : m)
      if (m_id == x.first)
        for (auto p : x.second)
          m_toimport[ h ].push_back( p );
    ++h;
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
    //m_u[ m_gid[i] ] = ansol_shear( i, 2400.0 );
    m_u[ m_gid[i] ] = ansol_gauss( i );

  m_lsmproxy.ckLocalBranch()->charesol( thisIndex, m_u );
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
    auto& lA = m_lhs[A];
    lA[A] += 2.0*J;
    lA[B] += J;
    lA[C] += J;
    lA[D] += J;
    auto& lB = m_lhs[B];
    lB[A] += J;
    lB[B] += 2.0*J;
    lB[C] += J;
    lB[D] += J;
    auto& lC = m_lhs[C];
    lC[A] += J;
    lC[B] += J;
    lC[C] += 2.0*J;
    lC[D] += J;
    auto& lD = m_lhs[D];
    lD[A] += J;
    lD[B] += J;
    lD[C] += J;
    lD[D] += 2.0*J;
  }

  m_timestamp.emplace_back( "Compute left-hand side matrix", t.dsec() );

  m_lsmproxy.ckLocalBranch()->charelhs( thisIndex, m_lhs );
}

void
Performer::rhs( tk::real mult,
                tk::real dt,
                std::map< std::size_t, tk::real >& unk )
//******************************************************************************
// Compute right-hand side of PDE
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  const tk::real U0 = 0.5;
  const tk::real LAMBDA = 5.0e-4;

  m_rhs.clear();

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
    vel[0][0] = U0 + LAMBDA * y[a];  vel[1][0] = 0.0;  vel[2][0] = 0.0;
    vel[0][1] = U0 + LAMBDA * y[b];  vel[1][1] = 0.0;  vel[2][1] = 0.0;
    vel[0][2] = U0 + LAMBDA * y[c];  vel[1][2] = 0.0;  vel[2][2] = 0.0;
    vel[0][3] = U0 + LAMBDA * y[d];  vel[1][3] = 0.0;  vel[2][3] = 0.0;

    // shape function derivatives
    std::array< std::array< tk::real, 3 >, 4 > grad;  // nnode*ndim [4][3]
    grad[1] = tk::crossdiv( ca, da, J );
    grad[2] = tk::crossdiv( ba, da, J );
    grad[3] = tk::crossdiv( ba, ca, J );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

    // add mass contribution to rhs
    for (std::size_t i=0; i<4; ++i) {
      const auto g = m_gid[ m_inpoel[e*4+i] ];  // global node id
      const auto& u = unk[g];                   // ref to unknown at g
      auto& r = m_rhs[g];                       // ref to rhs at g
      for (std::size_t j=0; j<4; ++j)           // add contribution to rhs
        r += mass[i][j] * 1.5;//u;
    }

//     // add advection contribution to rhs
//     for (std::size_t i=0; i<4; ++i) {
//       const auto g = m_gid[ m_inpoel[e*4+i] ];  // global node id
//       auto& r = m_rhs[g];                       // ref to rhs at g
//       for (std::size_t j=0; j<4; ++j)           // add contribution to rhs
//         for (std::size_t k=0; k<3; ++k)
//           for (std::size_t l=0; l<4; ++l)
//             r -= mult*dt*( mass[i][j] * vel[k][j] * grad[l][k] *
//                            unk[ m_gid[ m_inpoel[e*4+l] ] ] );
//     }
// 
//     // add diffusion contribution to rhs
//     for (std::size_t i=0; i<4; ++i) {
//       const auto g = m_gid[ m_inpoel[e*4+i] ];  // global node id
//       auto& r = m_rhs[g];                       // ref to rhs at g
//       for (std::size_t j=0; j<4; ++j) {         // add contribution to rhs
//         const auto& u = unk[ m_gid[ m_inpoel[e*4+j] ] ];  // ref to unknown
//         for (std::size_t k=0; k<3; ++k)
//           r -= mult*dt*( grad[i][k] * grad[j][k] * u );
//       }
//     }
  }

  m_timestamp.emplace_back( "Compute right-hand side vector", t.dsec() );

  m_lsmproxy.ckLocalBranch()->charerhs( thisIndex, m_rhs );
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

  // Build unique global node ids of owned elements
  for (auto e : gelem) {
    m_gid.push_back( g_tetinpoel[e*4] );
    m_gid.push_back( g_tetinpoel[e*4+1] );
    m_gid.push_back( g_tetinpoel[e*4+2] );
    m_gid.push_back( g_tetinpoel[e*4+3] );
  }
  tk::unique( m_gid );

  // Assign local node ids to global node ids
  const auto lnode = assignLid( m_gid );

  // Generate element connectivity for owned elements using local point ids
  for (auto e : gelem) {
    m_inpoel.push_back( tk::lid( lnode, g_tetinpoel[e*4] ) );
    m_inpoel.push_back( tk::lid( lnode, g_tetinpoel[e*4+1] ) );
    m_inpoel.push_back( tk::lid( lnode, g_tetinpoel[e*4+2] ) );
    m_inpoel.push_back( tk::lid( lnode, g_tetinpoel[e*4+3] ) );
  }

  // Send off number of columns per row to linear system merger
  m_lsmproxy.ckLocalBranch()->charerow( thisProxy, thisIndex, m_gid );

  m_timestamp.emplace_back( "Initialize mesh point ids, element connectivity",
                            t.dsec() );
}

std::map< std::size_t, std::size_t >
Performer::assignLid( const std::vector< std::size_t >& gid ) const
//******************************************************************************
//! Assign local ids to global ids
//! \param[in] gid
//! \return Map associating global ids to local ids
//! \author J. Bakosi
//******************************************************************************
{
  std::map< std::size_t, std::size_t > lid;
  std::size_t l = 0;
  for (auto p : gid) lid[p] = l++;
  return lid;
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
  if (!g_meshfilemap.empty())
    for (auto p : m_gid) er.readNode( g_meshfilemap[p], x, y, z );
  else
    for (auto p : m_gid) er.readNode( p, x, y, z );

  m_timestamp.emplace_back( "Read mesh point coordinates from file", t.dsec() );
}

void
Performer::updateSolution( const std::map< std::size_t, tk::real >& sol )
//******************************************************************************
// Update solution vector
//! \author J. Bakosi
//******************************************************************************
{
  for (const auto& r : sol) m_uf[ r.first ] = r.second;

  // If all contributions we own have been received, advance time step
  if (m_uf.size() == m_gid.size()) {

    writeFields( m_uf );
    m_uf.clear();

    // Tell the Charm++ runtime system to call back to Conductor::registered()
    // once all Performer chares have registered themselves, i.e., checked in,
    // with their local branch of the linear system merger group, LinSysMerger.
    // The reduction is done via creating a callback that invokes the typed
    // reduction client, where m_hostproxy is the proxy on which the reduction
    // target method, registered(), is called upon completion of the reduction.
    contribute(
      CkCallback( CkReductionTarget( Conductor, evaluateTime ), m_hostproxy ) );
  }
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
Performer::writeChareId( const tk::ExodusIIMeshWriter& ew ) const
//******************************************************************************
// Output chare id field to file
//! \param[in] ew ExodusII mesh-based writer object
//! \author J. Bakosi
//******************************************************************************
{
  // Write elem chare id field to mesh
  std::vector< tk::real > chid( m_inpoel.size()/4, static_cast<tk::real>(m_id) );
  ew.writeElemScalar( m_it+1, 1, chid );
}

void
Performer::writeSolution( const tk::ExodusIIMeshWriter& ew,
                          const std::map< std::size_t, tk::real >& u ) const
//******************************************************************************
// Output solution to file
//! \param[in] ew ExodusII mesh-based writer object
//! \author J. Bakosi
//******************************************************************************
{
  std::vector< tk::real > sol;
  for (const auto& p : u) sol.push_back( p.second );
  ew.writeNodeScalar( m_it+1, 1, sol );
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
  ew.writeNodeVarNames( { "Scalar" } );
}

void
Performer::writeFields( const std::map< std::size_t, tk::real >& u )
//******************************************************************************
// Output mesh-based fields to file
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  // Create ExodusII writer
  tk::ExodusIIMeshWriter
    ew( g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "_" +
          std::to_string( m_id ),
        tk::ExoWriter::OPEN );

  // Write time stamp
  ew.writeTimeStamp( m_it+1, m_t );

  // Write mesh-based fields
  writeChareId( ew );
  writeSolution( ew, u );

  m_timestamp.emplace_back( "Write mesh-based fields to file", t.dsec() );
}

void
Performer::advance( tk::real dt, uint64_t it, tk::real t )
//******************************************************************************
// Advance equations in time
//! \param[in] dt Size of time step
//! \param[in] it Iteration count
//! \param[in] t Physical time
//! \author J. Bakosi
//******************************************************************************
{
  // Update physical time and iteration count
  m_t = t;
  m_it = it;

  // Advance equations one step in time

  // Compute righ-hand side of PDE
  m_lsmproxy.ckLocalBranch()->enable_wait4rhs();

  rhs( 0.5, dt, m_u );

}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "performer.def.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
