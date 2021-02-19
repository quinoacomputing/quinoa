// *****************************************************************************
/*!
  \file      src/IO/RDGFLOMeshReader.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     RDGFLO mesh reader class declaration
  \details   RDGFLO mesh reader class declaration. Mesh reader facilitating
             reading a mesh from a simple text file used by Prof. Hong Luo at
             North Carolina State University.
*/
// *****************************************************************************

#include <array>
#include <istream>
#include <string>
#include <vector>
#include <cstddef>

#include "Types.hpp"
#include "Exception.hpp"
#include "UnsMesh.hpp"
#include "Reorder.hpp"
#include "RDGFLOMeshReader.hpp"

using tk::RDGFLOMeshReader;

void
RDGFLOMeshReader::readHeader()
// *****************************************************************************
//  Read RDGFLO mesh header
// *****************************************************************************
{
  // read in header: "npoin, ntetr, npyra, npris, nhexa, ntria, nquad, time"
  std::string s;
  for (int i=0; i<8; ++i) m_inFile >> s;

  // number of points
  m_inFile >> m_nnode;

  // number of tetrahedra
  m_inFile >> m_ntet;
  if (m_ntet == 0) Throw( "No tetrahedra in input mesh" );

  // number of pyramids
  std::size_t npyra;
  m_inFile >> npyra;
  if (npyra > 0) Throw( "Pyramids not supported. Need a mesh with only tetrahedra." );

  // number of prisms
  std::size_t npris;
  m_inFile >> npris;
  if (npris > 0) Throw( "Prisms not supported. Need a mesh with only tetrahedra." );

  // number of hexahedra
  std::size_t nhexa;
  m_inFile >> nhexa;
  if (nhexa > 0) Throw( "Hexahedra not supported. Need a mesh with only tetrahedra." );

  // number of triangles
  m_inFile >> m_ntri;

  // number of quads
  std::size_t nquad;
  m_inFile >> nquad;
  if (nquad > 0) Throw( "Quads not supported. Need a mesh with only tetrahedra." );

  // time
  tk::real time;
  m_inFile >> time;
}

void
RDGFLOMeshReader::readMesh( UnsMesh& mesh )
// *****************************************************************************
//  Read RDGFLO mesh
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  // Read header
  readHeader();
  // Read elements
  readElements( mesh );
  // Read nodes
  readNodes( mesh );
}

void
RDGFLOMeshReader::readNodes( UnsMesh& mesh )
// *****************************************************************************
//  Read nodes
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  std::string s;

  std::getline( m_inFile, s );
  std::getline( m_inFile, s );
  Assert( s == " grid point coordinates",
          "Invalid keyword, expected: ' grid point coordinates'" );

  // Read in node coordinates: x-coord y-coord z-coord
  auto& xcoord = mesh.x();
  auto& ycoord = mesh.y();
  auto& zcoord = mesh.z();
  xcoord.resize( m_nnode );
  ycoord.resize( m_nnode );
  zcoord.resize( m_nnode );
  for (std::size_t i=0; i<m_nnode; ++i) {
    std::size_t id;
    tk::real x, y, z;
    m_inFile >> id >> x >> y >> z;
    --id;       // convert to zero-based node ids
    xcoord[ id ] = x;
    ycoord[ id ] = y;
    zcoord[ id ] = z;
  }
}

void
RDGFLOMeshReader::readElements( UnsMesh& mesh )
// *****************************************************************************
//  Read element connectivity
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  std::string s;

  std::getline( m_inFile, s );
  std::getline( m_inFile, s );
  Assert( s == " element connectivity",
          "Invalid keyword, expected: ' element connectivity'" );

  // Read in tetrahedra element connectivity
  auto& tetinpoel = mesh.tetinpoel();
  tetinpoel.resize( m_ntet * 4 );
  for (std::size_t i=0; i<m_ntet; ++i) {
    std::size_t id;
    std::array< std::size_t, 4 > n;
    // ignore cell id, a, b
    m_inFile >> id >> n[0] >> n[1] >> n[2] >> n[3];
    --id;       // convert to zero-based element ids
    tetinpoel[ id*4+0 ] = n[0]-1;       // store zero-based node ids
    tetinpoel[ id*4+1 ] = n[1]-1;
    tetinpoel[ id*4+2 ] = n[2]-1;
    tetinpoel[ id*4+3 ] = n[3]-1;
  }

  std::getline( m_inFile, s );
  std::getline( m_inFile, s );
  Assert( s == " boundary conditions & connectivity for boundary faces",
          "Invalid keyword, expected: ' boundary conditions & connectivity "
          "for boundary faces'" );

  // Read in triangle element connectivity and side set ids
  auto& triinpoel = mesh.triinpoel();
  auto& bface = mesh.bface();
  auto& faceid = mesh.faceid();
  triinpoel.resize( m_ntri * 3 );
  for (std::size_t i=0; i<m_ntri; ++i) {
    std::array< std::size_t, 15 > n;
    for (std::size_t j=0; j<n.size(); ++j)  m_inFile >> n[j];
    --n[0];     // convert to zero-based node ids
    triinpoel[ n[0]*3+0 ] = n[6]-1;       // store zero-based node ids
    triinpoel[ n[0]*3+1 ] = n[7]-1;
    triinpoel[ n[0]*3+2 ] = n[8]-1;
    auto setid = static_cast< int >( n[1] );
    bface[ setid ].push_back( m_ntet + n[0] );  // tets written first
    faceid[ setid ].push_back( 0 );
  }
}
