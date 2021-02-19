// *****************************************************************************
/*!
  \file      src/IO/UGRIDMeshReader.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     UGRID mesh reader class declaration
  \details   UGRID mesh reader class declaration. Mesh reader facilitating
             reading a mesh from a simple text file used by NASA.
  \see       http://www.simcenter.msstate.edu/software/downloads/doc/ug_io/3d_grid_file_type_ugrid.html, http://www.simcenter.msstate.edu/software/downloads/doc/aflr3/aflr3_io_summary.html
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
#include "UGRIDMeshReader.hpp"

using tk::UGRIDMeshReader;

void
UGRIDMeshReader::readHeader()
// *****************************************************************************
//  Read UGRID mesh header
// *****************************************************************************
{
  std::string s;

  // Number_of_Nodes
  m_inFile >> m_nnode;

  // Number_of_Surf_Trias
  m_inFile >> m_ntri;

  // Number_of_Surf_Quads
  int nquad;
  m_inFile >> nquad;

  // Number_of_Vol_Tets
  m_inFile >> m_ntet;

  // Number_of_Vol_Pents_5
  int pent5;
  m_inFile >> pent5;

  // Number_of_Vol_Pents_6
  int pent6;
  m_inFile >> pent6;

  // Number_of_Vol_Hexs
  int nhex;
  m_inFile >> nhex;
}

void
UGRIDMeshReader::readMesh( UnsMesh& mesh )
// *****************************************************************************
//  Read UGRID mesh
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  // Read header
  readHeader();
  // Read nodes
  readNodes( mesh );
  // Read elements
  readElements( mesh );
}

void
UGRIDMeshReader::readNodes( UnsMesh& mesh )
// *****************************************************************************
//  Read nodes
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  // Read in node coordinates: x-coord y-coord z-coord
  for (std::size_t i=0; i<m_nnode; ++i) {
    tk::real x, y, z;
    m_inFile >> x >> y >> z;
    mesh.x().push_back( x );
    mesh.y().push_back( y );
    mesh.z().push_back( z );
  }
}

void
UGRIDMeshReader::readElements( UnsMesh& mesh )
// *****************************************************************************
//  Read element connectivity
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  // Read in triangle element connectivity
  for (std::size_t i=0; i<m_ntri; ++i) {
    std::array< std::size_t, 3 > n;
    m_inFile >> n[0] >> n[1] >> n[2];
    mesh.triinpoel().push_back( n[0] );
    mesh.triinpoel().push_back( n[1] );
    mesh.triinpoel().push_back( n[2] );
  }

  // Read side sets of triangle elements
  for (std::size_t i=0; i<m_ntri; ++i) {
    int setid;
    m_inFile >> setid;
    mesh.bface()[ setid ].push_back( m_ntet + i );
    mesh.faceid()[ setid ].push_back( 0 );
  }

  // Read in tetrahedra element connectivity
  for (std::size_t i=0; i<m_ntet; ++i) {
    std::array< std::size_t, 4 > n;
    m_inFile >> n[0] >> n[1] >> n[2] >> n[3];
    mesh.tetinpoel().push_back( n[0] );
    mesh.tetinpoel().push_back( n[1] );
    mesh.tetinpoel().push_back( n[2] );
    mesh.tetinpoel().push_back( n[3] );
  }

  // Shift node IDs to start from zero
  shiftToZero( mesh.triinpoel() );
  shiftToZero( mesh.tetinpoel() );
}
