// *****************************************************************************
/*!
  \file      src/IO/Omega_h_MeshReader.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Omega_h mesh reader
  \details   Omega_h mesh reader class definition.
*/
// *****************************************************************************

//#include "Omega_h_library.hpp"

#include "Omega_h_MeshReader.h"

#include "Macro.h"

using tk::Omega_h_MeshReader;

void
Omega_h_MeshReader::readGraph( std::vector< std::size_t >& ginpoel,
                               int n, int m )
// *****************************************************************************
//  Read a chunk of the mesh graph (connectivity) from Omega h file
//! \param[in] n Total number of PEs
//! \param[in] m This PE
//! \param[in,out] ginpoel Vector to read tetrtahedron element connectivity of
//!    our chunk of the mesh into
// *****************************************************************************
{
IGNORE(ginpoel);
IGNORE(n);
IGNORE(m);
}

std::size_t
Omega_h_MeshReader::readHeader()
// *****************************************************************************
//  Read header from Omega_h
//! \return Number of nodes in mesh
// *****************************************************************************
{
  return 0;
}

std::array< std::vector< tk::real >, 3 >
Omega_h_MeshReader::readCoords( const std::vector< std::size_t >& gid ) const
// *****************************************************************************
//  Read coordinates of a number of mesh nodes from Omega h file
//! \param[in] gid Global node IDs whose coordinates to read
//! \return Vector of node coordinates read from file
// *****************************************************************************
{
IGNORE(gid);
  std::array< std::vector< tk::real >, 3 > coords;
  return coords;
}

std::size_t
Omega_h_MeshReader::readSidesetFaces(
  std::map< int, std::vector< std::size_t > >& bface,
  std::map< int, std::vector< int > >& faceid )
// *****************************************************************************
//  Read face list of all side sets from Omega_h file
//! \param[in,out] bface Face-Element lists mapped to side set ids
//! \param[in,out] faceid Side set side lists associated to side set ids
//! \return Total number of boundary faces
// *****************************************************************************
{
IGNORE(bface);
IGNORE(faceid);
  return 0;
}

void
Omega_h_MeshReader::readFaces( std::size_t nbfac,
                               std::vector< std::size_t >& conn ) const
// *****************************************************************************
//  Read face connectivity of a number of boundary faces from Omega_h file
//! \param[in] nbfac Number of boundary faces
//! \param[inout] conn Connectivity vector to push to
//! \details This function reads in the total number of boundary faces,
//!   also called triangle-elements, and their connectivity.
//! \note It is okay to call this function with zero nbfac: it will be no-op.
// *****************************************************************************
{
  // Return if no boundary faces in file
  if (nbfac == 0) return;
IGNORE(conn);
}

std::map< int, std::vector< std::size_t > >
Omega_h_MeshReader::readSidesets()
// *****************************************************************************
//  Read node list of all side sets from Omega_h file
//! \return Node lists mapped to side set ids
// *****************************************************************************
{
  // Node lists mapped to side set ids
  std::map< int, std::vector< std::size_t > > side;

  return side;
}
