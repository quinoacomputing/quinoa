//******************************************************************************
/*!
  \file      src/Inciter/Performer.C
  \author    J. Bakosi
  \date      Mon 11 May 2015 02:48:00 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Performer advances the Euler equations
  \details   Performer advances the Euler equations. There are a potentially
    large number of Performer Charm++ chares created by Conductor. Each
    performer gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances the Euler equations in time.
*/
//******************************************************************************

#include <Performer.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <performer.def.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::pair< std::vector<std::size_t>, std::vector<std::size_t> > g_esup;
extern std::vector< std::size_t > g_tetinpoel;
extern std::vector< std::vector< std::size_t > > g_point;
extern std::vector< std::vector< std::size_t > > g_element;
extern std::vector< std::map< std::size_t, std::vector<std::size_t> > > g_comm;

} // inciter::

using inciter::Performer;

Performer::Performer( CProxy_Conductor& hostproxy,
                      LinSysMergerProxy& lsmproxy ) :
  m_id( static_cast< std::size_t >( thisIndex ) ),
  m_hostproxy( hostproxy ),
  m_lsmproxy( lsmproxy ),
  m_point( g_point[ m_id ] ),
  m_export( g_comm.size() > m_id ? g_comm[ m_id ] :
            std::map< std::size_t, std::vector< std::size_t > >() )
//******************************************************************************
// Constructor
//! \param[in] hostproxy Host proxy
//! \param[in] lsmproxy Linear system merger (LinSysMerger) proxy
//! \author J. Bakosi
//******************************************************************************
{
  // Contribute our matrix nonzero structure to the linear system merger
  m_lsmproxy.ckLocalBranch()->charenz( m_point, psup() );

  // Take over global mesh point ids of owned nodes
  std::vector< std::size_t > gelem( g_element[ m_id ] );
  // Initialize local->global, global->local node ids, element connectivity
  std::vector< std::size_t > gnode, inpoel;
  std::tie( gnode, inpoel ) = initIds( gelem );
  // Read coordinates of owned and received mesh nodes
  auto coord = initCoords( gnode );
  // Output chare mesh and nodal chare id field to file
  writeChareId( inpoel, coord );
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

//   // Add received node ids to local->global and global->local node ids
//   for (auto& m : m_import)
//     for (auto p : m.second) {
//       //m_point.push_back( p );
//       lnode[p] = l++;
//     }
//   Assert( m_import.empty() ? lnode.size() == m_point.size() : true,
//           "The size of global->local and local->global node id map must "
//           "be equal if there is only a single chare" );

  return { gnode, inpoel };
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
Performer::psup()
//******************************************************************************
//! Compute points surrounding points owned
//! \return Linked vectors storing points surrounding points owned
//! \details This function computes a special version of points surrounding
//!   points, derived from the mesh connectivity. The data structure returned
//!   contains only those points surrounding points that we own. Note that we
//!   only store points that we own, but of course, their surrounding points
//!   might list point ids that we do not own. This data structure describes our
//!   portion of the nonzero structure of a distributed matrix.
//! \author J. Bakosi
//******************************************************************************
{
  // Build unique global point ids of elements with at least one owned point
  std::vector< std::size_t > gnode;
  for (auto p : m_point)
    for (auto i=g_esup.second[p]+1; i<=g_esup.second[p+1]; ++i) {
      auto e = g_esup.first[i];
      if (g_tetinpoel[e*4+0] == p || g_tetinpoel[e*4+1] == p ||
          g_tetinpoel[e*4+2] == p || g_tetinpoel[e*4+3] == p) {
        gnode.push_back( g_tetinpoel[e*4+0] );
        gnode.push_back( g_tetinpoel[e*4+1] );
        gnode.push_back( g_tetinpoel[e*4+2] );
        gnode.push_back( g_tetinpoel[e*4+3] );
      }
    }
  tk::unique( gnode );

  // Assign local point ids to global point ids
  const auto lnode = assignLid( gnode );

  // Get element connectivity of those containing at least one owned point
  std::vector< std::size_t > inpoel;
  for (auto p : m_point)
    for (auto i=g_esup.second[p]+1; i<=g_esup.second[p+1]; ++i) {
      auto e = g_esup.first[i];
      if (g_tetinpoel[e*4+0] == p || g_tetinpoel[e*4+1] == p ||
          g_tetinpoel[e*4+2] == p || g_tetinpoel[e*4+3] == p) {
        inpoel.push_back( lid( lnode, g_tetinpoel[e*4+0] ) );
        inpoel.push_back( lid( lnode, g_tetinpoel[e*4+1] ) );
        inpoel.push_back( lid( lnode, g_tetinpoel[e*4+2] ) );
        inpoel.push_back( lid( lnode, g_tetinpoel[e*4+3] ) );
      }
    }

  // Generate points surrounding points based on inpoel with local ids
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );
  auto& psup1 = psup.first;
  auto& psup2 = psup.second;

  // Lambda to find out if a point is owned
  auto own = [&]( std::size_t gid ) -> bool {
    for (auto p : m_point) if (p == gid) return true;
    return false;
  };

  // Create new points surrounding points of owned points only
  std::vector< std::size_t > p1( 1, 0 ), p2( 1, 0 );
  std::size_t k = 0;
  for (std::size_t p=0; p<psup2.size()-1; ++p)
    if ( own( gnode[p] ) ) {
      for (auto i=psup2[p]+1; i<=psup2[p+1]; ++i)
        p1.push_back( psup1[i] );
      p2.push_back( p2.back() + psup2[p+1] - psup2[p] );
    }
  psup1 = std::move( p1 );
  psup2 = std::move( p2 );

  // Convert local to global point ids in points surrounding points
  for (auto& p : psup1) p = gnode[ p ];

  return psup;
}

std::map< std::size_t, std::size_t >
Performer::assignLid( const std::vector< std::size_t >& gid )
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
                std::size_t gid )
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
  tk::UnsMesh mesh;
  tk::ExodusIIMeshReader
    er( g_inputdeck.get< tag::cmd, tag::io, tag::input >(), mesh );

  std::vector< tk::real > x, y, z;
  for (auto p : gnode) er.readNode( p, x, y, z );

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
 }
