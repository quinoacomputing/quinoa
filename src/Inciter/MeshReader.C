// *****************************************************************************
/*!
  \file      src/Inciter/MeshReader.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Mesh reader class for inciter connecting to various readers
  \details   Mesh reader ckass for inciter connecting to various lower level
             mesh readers.
*/
// *****************************************************************************

#include <numeric>

#include "MeshReader.h"

using inciter::MeshReader;

void
MeshReader::readGraph( std::vector< std::size_t >& ginpoel, int n, int m )
// *****************************************************************************
//  Read our chunk of the mesh graph (connectivity) from file
//! \param[in] n Total number of PEs
//! \param[in] m This PE
//! \param[in,out] ginpoel Vector to read tetrtahedron element connectivity of
//!    our chunk of the mesh into
// *****************************************************************************
{
  // Get number of mesh points and number of tetrahedron elements in file
  m_er.readElemBlockIDs();
  auto nel = m_er.nelem( tk::ExoElemType::TET );

  // Read our contiguously-numbered chunk of tetrahedron element
  // connectivity from file and also generate and store the list of global
  // element indices for our chunk of the mesh

  // Compute extents of element IDs of our mesh chunk to read
  auto npes = static_cast< std::size_t >( n );
  auto mype = static_cast< std::size_t >( m );
  auto chunk = nel / npes;
  auto from = mype * chunk;
  auto till = from + chunk;
  if (mype == npes-1) till += nel % npes;

  // Read tetrahedron connectivity between from and till
  m_er.readElements( {{from, till-1}}, tk::ExoElemType::TET, ginpoel );
}

std::array< std::vector< tk::real >, 3 >
MeshReader::readCoords( const std::vector< std::size_t > gid )
// *****************************************************************************
//  Read coordinates of a number of mesh nodes from ExodusII file
//! \param[in] gid Global node IDs whose coordinates to read
//! \return Vector of node coordinates read from file
// *****************************************************************************
{
  Assert( !gid.empty(), "Global ID vector empty" );

  // Read node coordinates from file with global node IDs given in gid
  return m_er.readNodes( gid );
}
