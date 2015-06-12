//******************************************************************************
/*!
  \file      src/Inciter/Performer.C
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 09:59:45 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Performer advances the Euler equations
  \details   Performer advances the Euler equations. There are a potentially
    large number of Performer Charm++ chares created by Conductor. Each
    performer gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances the Euler equations in time.
*/
//******************************************************************************

#include <string>
#include <iterator>
#include <cmath>

#include "Constants.h"
#include "Exception.h"
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
extern std::pair< std::vector<std::size_t>, std::vector<std::size_t> > g_esup;
extern std::vector< std::size_t > g_tetinpoel;
extern std::vector< std::size_t > g_meshfilemap;
extern std::vector< std::vector< std::size_t > > g_point;
extern std::vector< std::vector< std::size_t > > g_element;
extern std::vector< std::map< std::size_t, std::vector<std::size_t> > > g_ecomm;

} // inciter::

using inciter::Performer;

Performer::Performer( CProxy_Conductor& hostproxy, LinSysMergerProxy& lsmproxy )
: m_id( static_cast< std::size_t >( thisIndex ) ),
  m_hostproxy( hostproxy ),
  m_lsmproxy( lsmproxy ),
  m_it( 0 ),
  m_t( 0.0 ),
  m_point( g_point[ m_id ] )
//******************************************************************************
// Constructor
//! \param[in] hostproxy Host proxy
//! \param[in] lsmproxy Linear system merger (LinSysMerger) proxy
//! \author J. Bakosi
//******************************************************************************
{
  // Activate SDAG waits
  wait4mass();
  wait4rhs();
  // Initialize import maps
  initImports();
  // Take over global mesh element ids of owned elements
  std::vector< std::size_t > gelem( g_element[ m_id ] );
  // Initialize local->global, global->local node ids, element connectivity
  std::vector< std::size_t > gnode, inpoel;
  std::tie( gnode, inpoel ) = initIds( gelem );
  // Read coordinates of owned and received mesh nodes
  auto coord = initCoords( gnode );
  // Output chare mesh to file
  writeMesh( inpoel, coord );
  // Output mesh-based fields metadata to file
  writeMeta();
  // Set initial conditions
  ic( gnode, coord );
  // Compute left-hand side of PDE
  lhs( gnode, inpoel, coord );
  // Compute righ-hand side of PDE
  rhs( gnode, inpoel, coord );
  // Output field data to file
  writeFields();
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

//   for (const auto& m : m_toimport) {
//     std::cout << m_id << " will import from " << m.first << ": ";
//     for (auto p : m.second)
//       std::cout << p << " ";
//     std::cout << '\n';
//   }
}

void
Performer::ic( const std::vector< std::size_t >& gnode,
               const std::array< std::vector< tk::real >, 3 >& coord )
//******************************************************************************
// Set initial conditions
//! \param[in] gnode Global node ids of owned elements
//! \param[in] coord Mesh point coordinates
//! \author J. Bakosi
//******************************************************************************
{
  // Lambda to compute initial conditions for dispersion in simple shear flow
  auto ansol_shear = [ &coord ]( std::size_t i, tk::real t ){
    const tk::real X0 = 0.5;           // x position of source
    const tk::real Y0 = 0.5;            // y position of source
    const tk::real U0 = 0.5;            // velocity in x direction
    const tk::real t0 = 2.4;            // initial time
    const tk::real LAMBDA = 5.0e-4;     // amount of shear, (lambda = du/dy)
    const tk::real D = 10.0;            // scalar diffusivity
    const tk::real MASS = 4 * tk::pi<tk::real>() * t0 *
                          std::sqrt( 1.0 + LAMBDA * LAMBDA * t0 * t0 / 12.0 );

    const auto& x = coord[0];
    const auto& y = coord[1];

    tk::real a = x[i] - X0 - U0*t - 0.5*LAMBDA*y[i]*t;
    tk::real b = 1.0 + LAMBDA*LAMBDA*t*t/12.0;
    return MASS * exp( -a*a/(4*D*t*b) - y[i]*y[i]/(4*D*t) )
                / ( 4.0 * tk::pi<tk::real>() * t * std::sqrt(b) );
  };

  // Lambda to compute initial conditions representing a joint Gaussian
  auto ansol_gauss = [ &coord ]( std::size_t i ){
    const tk::real X0 = 0.5, Y0 = 0.5, Z0 = 0.5;
    const tk::real v1 = 0.1, v2 = 0.05, v3 = 0.05;

    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];

    const tk::real rx = x[i]-X0;
    const tk::real ry = y[i]-Y0;
    const tk::real rz = z[i]-Z0;

    return std::exp( -0.5 * (rx*rx/v1 + ry*ry/v2 + rz*rz/v3) ) /
           2.0 / tk::pi<tk::real>() / v1 / v2 / v3;
  };

  for (std::size_t i=0; i<gnode.size(); ++i)
    //m_x[ gnode[i] ] = ansol_shear( i, 2.4 );
    m_x[ gnode[i] ] = ansol_gauss( i );
}

void
Performer::lhs(
  const std::vector< std::size_t >& gnode,
  const std::vector< std::size_t >& inpoel,
  const std::array< std::vector< tk::real >, 3 >& coord )
//******************************************************************************
// Compute left-hand side of PDE
//! \param[in] gnode Global node ids of owned elements
//! \param[in] inpoel Tetrahedron element connectivity of owned elements
//! \param[in] coord Mesh point coordinates
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const auto a = inpoel[e*4+0];
    const auto b = inpoel[e*4+1];
    const auto c = inpoel[e*4+2];
    const auto d = inpoel[e*4+3];
    std::array< tk::real, 3 > ba{{ x[b]-x[a], y[b]-y[a], z[b]-z[a] }},
                              ca{{ x[c]-x[a], y[c]-y[a], z[c]-z[a] }},
                              da{{ x[d]-x[a], y[d]-y[a], z[d]-z[a] }};
    const auto J = tk::triple( ba, ca, da ) / 120.0;

    const auto A = gnode[a];
    const auto B = gnode[b];
    const auto C = gnode[c];
    const auto D = gnode[d];
    auto& lshA = m_lhs[A];
    lshA[A] += 2.0*J;
    lshA[B] += J;
    lshA[C] += J;
    lshA[D] += J;
    auto& lhsB = m_lhs[B];
    lhsB[A] += J;
    lhsB[B] += 2.0*J;
    lhsB[C] += J;
    lhsB[D] += J;
    auto& lhsC = m_lhs[C];
    lhsC[A] += J;
    lhsC[B] += J;
    lhsC[C] += 2.0*J;
    lhsC[D] += J;
    auto& lhsD = m_lhs[D];
    lhsD[A] += J;
    lhsD[B] += J;
    lhsD[C] += J;
    lhsD[D] += 2.0*J;
  }

  m_timestamp.emplace_back( "Compute left-hand side matrix", t.dsec() );

  // Update left-hand side matrix of PDE by communicating with fellow chares,
  // and send off our contribution to linear system merger
  commLhs( g_ecomm.size() > m_id ? g_ecomm[ m_id ] :
           std::map< std::size_t, std::vector< std::size_t > >() );
}

void
Performer::commLhs(
  const std::map< std::size_t, std::vector< std::size_t > >& exp )
//******************************************************************************
//  Perform the necessary communication among fellow Performers to update the
//  chare-boundaries for left-hand side matrix of PDE
//! \param[in] exp Chare export map, see tk::elemCommMaps()
//! \author J. Bakosi
//******************************************************************************
{
  m_timer[ TimerTag::LHS ];  // start timer measuring the communication of lhs

  // Pack matrix values for export
  std::map< std::size_t,
            std::map< std::size_t,
                      std::map< std::size_t, tk::real > > > expmat;
  for (const auto& c : exp) {
    auto& e = expmat[ c.first ];
    for (auto p : c.second) {
      const auto it = m_lhs.find( p );
      if (it != m_lhs.end())
        e.insert( *it );
      else
        Throw( "Performer chare " + std::to_string(thisIndex) +
               " can't find gid " + std::to_string(p) +
               " to be exported in its part of the lhs matrix" );
    }
  }

  if (m_toimport.empty()) trigger_lhs_complete();

  // Export matrix contributions
  for (const auto& c : expmat)
    thisProxy[ static_cast<int>(c.first) ].addLhs( thisIndex, c.second );
}

void
Performer::addLhs(
  int id,
  const std::map< std::size_t, std::map< std::size_t, tk::real > >& rows )
//******************************************************************************
// Receive matrix rows contribution from fellow Performer chare
//! \author J. Bakosi
//******************************************************************************
{
//   std::cout << "import on " << thisIndex << ": ";
//   for (const auto& r : rows) {
//     std::cout << "(" << r.first << ") ";
//     //for (const auto& c : r.second)
//     //  std::cout << c.first << " ";
//   }
//   std::cout << "\n";

  // Add contributions received
  for (const auto& r : rows) {
    auto& localrow = m_lhs[ r.first ];
    m_lhsimport[ static_cast<std::size_t>(id) ].push_back( r.first );
    for (const auto& c : r.second) localrow[ c.first ] += c.second;
  }

//   for (const auto& m : m_lhsimport) {
//     std::cout << m.first << " % ";
//     for (auto p : m.second)
//       std::cout << p << " ";
//     std::cout << '\n';
//   }

  if (lhscomplete()) trigger_lhs_complete();
}

void
Performer::contributeLhs()
//******************************************************************************
// Contribute our portion of the left-hand side matrix
//! \author J. Bakosi
//******************************************************************************
{
  // Keep only rows owned before send
  for (auto it=begin(m_lhs); it!=end(m_lhs); )
    if (!own( it->first ))
      it = m_lhs.erase( it );
    else
      ++it;

//   std::cout << thisIndex << ": ";
//   for (const auto& r : m_lhs) {
//     std::cout << r.first << ", ";
//     for (const auto& c : r.second)
//       std::cout << c.first << " ";
//   }
//   std::cout << '\n';

  m_lsmproxy.ckLocalBranch()->charelhs( m_lhs );
}

void
Performer::rhs(
  const std::vector< std::size_t >& gnode,
  const std::vector< std::size_t >& inpoel,
  const std::array< std::vector< tk::real >, 3 >& coord )
//******************************************************************************
// Compute right-hand side of PDE
//! \param[in] gnode Global node ids of owned elements
//! \param[in] inpoel Tetrahedron element connectivity of owned elements
//! \param[in] coord Mesh point coordinates
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const auto a = inpoel[e*4+0];
    const auto b = inpoel[e*4+1];
    const auto c = inpoel[e*4+2];
    const auto d = inpoel[e*4+3];

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

    // prescribed rotational velocity: [ 0.5-y, x-0.5, 0 ]
    std::array< std::array< tk::real, 4 >, 3 > vel;  // ndim*nnode [3][4]
    vel[0][0] = 0.5-y[a];  vel[1][0] = x[a]-0.5;  vel[2][0] = 0.0;
    vel[0][1] = 0.5-y[b];  vel[1][1] = x[b]-0.5;  vel[2][1] = 0.0;
    vel[0][2] = 0.5-y[c];  vel[1][2] = x[c]-0.5;  vel[2][2] = 0.0;
    vel[0][3] = 0.5-y[d];  vel[1][3] = x[d]-0.5;  vel[2][3] = 0.0;

    // shape function derivatives
    std::array< std::array< tk::real, 3 >, 4 > grad;  // nnode*ndim [4][3]
    grad[1] = tk::crossdiv( ca, da, J );
    grad[2] = tk::crossdiv( ba, da, J );
    grad[3] = tk::crossdiv( ba, ca, J );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

    // add mass contribution to rhs
    for (std::size_t i=0; i<4; ++i) {
      const auto g = gnode[ inpoel[e*4+i] ];    // global node id
      const auto& u = m_x[g];                   // ref to unknown at g
      auto& r = m_rhs[g];                       // ref to rhs at g
      for (std::size_t j=0; j<4; ++j)           // add contribution to rhs
        r += mass[i][j] * u;
    }

    // add advection contribution to rhs
    for (std::size_t i=0; i<4; ++i) {
      const auto g = gnode[ inpoel[e*4+i] ];    // global node id
      auto& r = m_rhs[g];                       // ref to rhs at g
      for (std::size_t j=0; j<4; ++j)           // add contribution to rhs
        for (std::size_t k=0; k<3; ++k)
          for (std::size_t l=0; l<4; ++l)
            r += mass[i][j] * vel[k][j] * grad[l][k] *
                   m_x[ gnode[ inpoel[e*4+l] ] ];
    }

    // add diffusion contribution to rhs
    for (std::size_t i=0; i<4; ++i) {
      const auto g = gnode[ inpoel[e*4+i] ];    // global node id
      auto& r = m_rhs[g];                       // ref to rhs at g
      for (std::size_t j=0; j<4; ++j) {         // add contribution to rhs
        const auto& u = m_x[ gnode[ inpoel[e*4+j] ] ];  // ref to unknown
        for (std::size_t k=0; k<3; ++k)
          r += grad[i][k] * grad[j][k] * u;
      }
    }
  }

  m_timestamp.emplace_back( "Compute right-hand side vector", t.dsec() );

  // Update right-hand side vector of PDE by communicating with fellow chares,
  // and send off our contribution to linear system merger
  commRhs( g_ecomm.size() > m_id ? g_ecomm[ m_id ] :
           std::map< std::size_t, std::vector< std::size_t > >() );
}

void
Performer::commRhs(
  const std::map< std::size_t, std::vector< std::size_t > >& exp )
//******************************************************************************
//  Perform the necessary communication among fellow Performers to update the
//  chare-boundaries for right-hand side vector of PDE
//! \param[in] exp Chare export map, see tk::elemCommMaps()
//! \author J. Bakosi
//******************************************************************************
{
  m_timer[ TimerTag::RHS ];  // start timer measuring the communication of rhs

  // Pack vector values for export
  std::map< std::size_t, std::map< std::size_t, tk::real > > expvec;
  for (const auto& c : exp) {
    auto& e = expvec[ c.first ];
    for (auto p : c.second) {
      const auto it = m_rhs.find( p );
      if (it != m_rhs.end())
        e.insert( *it );
      else
        Throw( "Performer chare " + std::to_string(thisIndex) +
               " can't find gid " + std::to_string(p) +
               " to be exported in its part of the rhs vector" );
    }
  }

  if (m_toimport.empty()) trigger_rhs_complete();

  // Export vector contributions
  for (const auto& c : expvec)
    thisProxy[ static_cast<int>(c.first) ].addRhs( thisIndex, c.second );
}

void
Performer::addRhs( int id, const std::map< std::size_t, tk::real >& rows )
//******************************************************************************
// Receive rhs vector rows contribution from fellow Performer chare
//! \author J. Bakosi
//******************************************************************************
{
  // Add contributions received
  for (const auto& r : rows) {
    m_rhs[ r.first ] += r.second;
    m_rhsimport[ static_cast<std::size_t>(id) ].push_back( r.first );
  }

  if (rhscomplete()) trigger_rhs_complete();
}

void
Performer::contributeRhs()
//******************************************************************************
// Contribute our portion of the right-hand side vector
//! \author J. Bakosi
//******************************************************************************
{
  // Keep only rows owned before send
  for (auto it=begin(m_rhs); it!=end(m_rhs); )
    if (!own( it->first ))
      it = m_rhs.erase( it );
    else
      ++it;

  m_lsmproxy.ckLocalBranch()->charerhs( m_rhs );
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
Performer::initIds( const std::vector< std::size_t >& gelem )
//******************************************************************************
//! Initialize local->global, global->local node ids, element connectivity
//! \param[in] gelem Set of unique owned global element ids
//! \return Unique global node ids of owned elements and tetrahedron element
//!   connectivity
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  // Build unique global node ids of owned elements
  std::vector< std::size_t > gnode;
  for (auto e : gelem) {
    gnode.push_back( g_tetinpoel[e*4] );
    gnode.push_back( g_tetinpoel[e*4+1] );
    gnode.push_back( g_tetinpoel[e*4+2] );
    gnode.push_back( g_tetinpoel[e*4+3] );
  }
  tk::unique( gnode );

  // Assign local node ids to global node ids
  const auto lnode = assignLid( gnode );

  // Generate element connectivity for owned elements using local point ids
  std::vector< std::size_t > inpoel;
  for (auto e : gelem) {
    inpoel.push_back( lid( lnode, g_tetinpoel[e*4] ) );
    inpoel.push_back( lid( lnode, g_tetinpoel[e*4+1] ) );
    inpoel.push_back( lid( lnode, g_tetinpoel[e*4+2] ) );
    inpoel.push_back( lid( lnode, g_tetinpoel[e*4+3] ) );
  }

  // Store number of owned elements
  m_nelem = inpoel.size() / 4;

  m_timestamp.emplace_back( "Initialize mesh point ids, element connectivity",
                            t.dsec() );
  return { gnode, inpoel };
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

std::size_t
Performer::lid( const std::map< std::size_t, std::size_t >& lnode,
                std::size_t gid ) const
//******************************************************************************
//! Find local for global node id
//! \param[in] lnode Global->local id map
//! \param[in] gid
//! \return Local id
//! \author J. Bakosi
//******************************************************************************
{
  const auto it = lnode.find( gid );

  if (it != lnode.end())
    return it->second;
  else
    Throw( "Can't find global node id " + std::to_string(gid) );
}

std::array< std::vector< tk::real >, 3 >
Performer::initCoords( const std::vector< std::size_t >& gnode )
//******************************************************************************
//  Read coordinates of mesh nodes from file
//! \param[in] gnode Global node ids of owned elements
//! \return Mesh point coordinates
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  tk::ExodusIIMeshReader
    er( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );

  std::vector< tk::real > x, y, z;
  if (!g_meshfilemap.empty())
    for (auto p : gnode) er.readNode( g_meshfilemap[p], x, y, z );
  else
    for (auto p : gnode) er.readNode( p, x, y, z );

  m_timestamp.emplace_back( "Read mesh point coordinates from file", t.dsec() );

  return { { x, y, z } };
}

void
Performer::writeMesh( const std::vector< std::size_t >& inpoel,
                      const std::array< std::vector< tk::real >, 3 >& coord )
//******************************************************************************
// Output chare mesh to file
//! \param[in] inpoel Tetrahedron element connectivity of owned elements
//! \param[in] coord Mesh point coordinates
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  // Create mesh object initializing element connectivity and point coords
  tk::UnsMesh mesh( inpoel, coord );

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
  std::vector< tk::real > chid( m_nelem, static_cast<tk::real>(m_id) );
  ew.writeElemScalar( m_it+1, 1, chid );
}

void
Performer::writeSolution( const tk::ExodusIIMeshWriter& ew ) const
//******************************************************************************
// Output solution to file
//! \param[in] ew ExodusII mesh-based writer object
//! \author J. Bakosi
//******************************************************************************
{
  std::vector< tk::real > x;
  for (const auto& p : m_x) x.push_back( p.second );
  ew.writeNodeScalar( m_it+1, 1, x );
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
Performer::writeFields()
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
  writeSolution( ew );

  m_timestamp.emplace_back( "Write mesh-based fields to file", t.dsec() );
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "performer.def.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
