//******************************************************************************
/*!
  \file      src/Inciter/Performer.C
  \author    J. Bakosi
  \date      Sat 30 May 2015 11:41:15 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Performer advances the Euler equations
  \details   Performer advances the Euler equations. There are a potentially
    large number of Performer Charm++ chares created by Conductor. Each
    performer gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances the Euler equations in time.
*/
//******************************************************************************

#include <string>
#include <tuple>
#include <iterator>

#include "Exception.h"
#include "Performer.h"
#include "Vector.h"
#include "ContainerUtil.h"
#include "UnsMesh.h"
#include "ExodusIIMeshReader.h"
#include "ExodusIIMeshWriter.h"
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
  // Initialize import maps
  initImports();
  // Take over global mesh point ids of owned nodes
  std::vector< std::size_t > gelem( g_element[ m_id ] );
  // Initialize local->global, global->local node ids, element connectivity
  std::vector< std::size_t > gnode, inpoel;
  std::tie( gnode, inpoel ) = initIds( gelem );
  // Read coordinates of owned and received mesh nodes
  auto coord = initCoords( gnode );
  // Compute left-hand side of PDE
  lhs( gnode, inpoel, coord );
  // Compute righ-hand side of PDE
  rhs( gnode, inpoel, coord );
  // Output chare mesh and nodal chare id field to file
  writeChareId( inpoel, coord );
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
Performer::lhs(
  const std::vector< std::size_t >& gnode,
  const std::vector< std::size_t >& inpoel,
  const std::array< std::vector< tk::real >, 3 >& coord )
//******************************************************************************
// Compute left-hand side of PDE
//! \param[in] gnode Global node ids whose coordinates to read from file
//! \param[in] inpoel Tetrahedron element connectivity
//! \param[in] coord Mesh point coordinates
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];
    const auto a = inpoel[e*4+0];
    const auto b = inpoel[e*4+1];
    const auto c = inpoel[e*4+2];
    const auto d = inpoel[e*4+3];
    std::array< tk::real, 3 > v1{{ x[b]-x[a], y[b]-y[a], z[b]-z[a] }},
                              v2{{ x[c]-x[a], y[c]-y[a], z[c]-z[a] }},
                              v3{{ x[d]-x[a], y[d]-y[a], z[d]-z[a] }};
    const auto vol = tk::triple( v1, v2, v3 ) / 6.0 / 20.0;
    const auto A = gnode[a];
    const auto B = gnode[b];
    const auto C = gnode[c];
    const auto D = gnode[d];
    auto& mA = m_lhs[A];
    mA[A] += 2.0*vol;
    mA[B] += vol;
    mA[C] += vol;
    mA[D] += vol;
    auto& mB = m_lhs[B];
    mB[A] += vol;
    mB[B] += 2.0*vol;
    mB[C] += vol;
    mB[D] += vol;
    auto& mC = m_lhs[C];
    mC[A] += vol;
    mC[B] += vol;
    mC[C] += 2.0*vol;
    mC[D] += vol;
    auto& mD = m_lhs[D];
    mD[A] += vol;
    mD[B] += vol;
    mD[C] += vol;
    mD[D] += 2.0*vol;
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
  m_timer.emplace_back();

  // Pack matrix values for export
  std::map< std::size_t,
            std::map< std::size_t,
                      std::map< std::size_t, tk::real > > > expmat;
  std::size_t ncommpts = 0;
  for (const auto& c : exp) {
    auto& e = expmat[ c.first ];
    for (auto p : c.second) {
      const auto it = m_lhs.find( p );
      if (it != m_lhs.end()) {
        e.insert( *it );
        ncommpts += it->second.size();
      } else
        Throw( "Performer chare " + std::to_string(thisIndex) +
               " can't find gid " + std::to_string(p) +
               " to be exported in its part of the matrix" );
    }
  }

  if (m_toimport.empty()) mass_complete();

  // Export matrix contributions
  for (const auto& c : expmat)
    thisProxy[ static_cast<int>(c.first) ].add( thisIndex, c.second );
}

void
Performer::add(
  int id,
  const std::map< std::size_t, std::map< std::size_t, tk::real > >& rows )
//******************************************************************************
// Receive matrix row contribution from fellow Performer chare
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
    m_import[ static_cast<std::size_t>(id) ].push_back( r.first );
    for (const auto& c : r.second) localrow[ c.first ] += c.second;
  }

//   for (const auto& m : m_import) {
//     std::cout << m.first << " % ";
//     for (auto p : m.second)
//       std::cout << p << " ";
//     std::cout << '\n';
//   }

  if (lhscomplete()) mass_complete();
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
//! \param[in] gnode Global node ids whose coordinates to read from file
//! \param[in] inpoel Tetrahedron element connectivity
//! \param[in] coord Mesh point coordinates
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];
    const auto a = inpoel[e*4+0];
    const auto b = inpoel[e*4+1];
    const auto c = inpoel[e*4+2];
    const auto d = inpoel[e*4+3];

    // compute element volume
    std::array< tk::real, 3 > v1{{ x[b]-x[a], y[b]-y[a], z[b]-z[a] }},
                              v2{{ x[c]-x[a], y[c]-y[a], z[c]-z[a] }},
                              v3{{ x[d]-x[a], y[d]-y[a], z[d]-z[a] }};
    const auto vol = tk::triple( v1, v2, v3 ) / 6.0 / 20.0;
    const auto A = gnode[a];
    const auto B = gnode[b];
    const auto C = gnode[c];
    const auto D = gnode[d];

    // construct tetrahedron element-level matrices

    // consistent mass
    std::array< tk::real, 16 > mass;                  // size: nnode*nnode (4*4)
    // diagonal
    mass[0*4+0] = mass[1*4+1] = mass[2*4+2] = mass[3*4+3] = 2.0*vol;
    // off-diagonal
    mass[0*4+1] = mass[0*4+2] = mass[0*4+3] =
    mass[1*4+0] = mass[1*4+2] = mass[1*4+3] =
    mass[2*4+0] = mass[2*4+1] = mass[2*4+3] =
    mass[3*4+0] = mass[3*4+1] = mass[3*4+2] = vol;

    // prescribed rotational velocity: [ 0.5-y, x-0.5, 0 ]
    std::array< tk::real, 12 > vel;                   // size: ndim*nnode (3*4)
    vel[0*4+0] = 0.5-y[a];  vel[1*4+0] = x[a]-0.5;  vel[2*4+0] = 0.0;
    vel[0*4+1] = 0.5-y[b];  vel[1*4+1] = x[b]-0.5;  vel[2*4+1] = 0.0;
    vel[0*4+2] = 0.5-y[c];  vel[1*4+2] = x[c]-0.5;  vel[2*4+2] = 0.0;
    vel[0*4+3] = 0.5-y[d];  vel[1*4+3] = x[d]-0.5;  vel[2*4+3] = 0.0;

    // shape function derivatives
    std::array< tk::real, 12 > der;                   // size: nnode*ndim (4*3)

    //...
    //std::array< tk::real, 4 > dh, dl, fh, fl, mh, ml; // size: nnode (4)

//     auto& mA = m_lhs[A];
//     mA[A] += 2.0*vol;
//     mA[B] += vol;
//     mA[C] += vol;
//     mA[D] += vol;
//     auto& mB = m_lhs[B];
//     mB[A] += vol;
//     mB[B] += 2.0*vol;
//     mB[C] += vol;
//     mB[D] += vol;
//     auto& mC = m_lhs[C];
//     mC[A] += vol;
//     mC[B] += vol;
//     mC[C] += 2.0*vol;
//     mC[D] += vol;
//     auto& mD = m_lhs[D];
//     mD[A] += vol;
//     mD[B] += vol;
//     mD[C] += vol;
//     mD[D] += 2.0*vol;

  }

  m_timestamp.emplace_back( "Compute right-hand side vector", t.dsec() );

//   // Update left-hand side matrix of PDE by communicating with fellow chares,
//   // and send off our contribution to linear system merger
//   commLhs( g_ecomm.size() > m_id ? g_ecomm[ m_id ] :
//            std::map< std::size_t, std::vector< std::size_t > >() );
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
//! \param[in] gnode Global node ids whose coordinates to read from file
//! \return Mesh point coordinates
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  tk::UnsMesh mesh;
  tk::ExodusIIMeshReader
    er( g_inputdeck.get< tag::cmd, tag::io, tag::input >(), mesh );

  std::vector< tk::real > x, y, z;
  if (!g_meshfilemap.empty())
    for (auto p : gnode) er.readNode( g_meshfilemap[p], x, y, z );
  else
    for (auto p : gnode) er.readNode( p, x, y, z );

  m_timestamp.emplace_back( "Read mesh point coordinates from file", t.dsec() );

  return { { x, y, z } };
}

void
Performer::writeChareId( const std::vector< std::size_t >& inpoel,
                         const std::array< std::vector< tk::real >, 3 >& coord )
//******************************************************************************
// Output chare mesh and nodal chare id field to file
//! \param[in] inpoel Tetrahedron element connectivity
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
        mesh );

  // Write chare mesh
  ew.write();
  ew.writeTimeStamp( 1, 1.0 );

  // Write elem chare id field to mesh
  std::vector< tk::real > chid( mesh.nelem(), static_cast<tk::real>(m_id) );
  ew.writeElemVarNames( { "Chare Id" } );
  ew.writeElemScalar( 1, 1, chid );

  m_timestamp.emplace_back( "Write chare id field to file", t.dsec() );
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "performer.def.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
