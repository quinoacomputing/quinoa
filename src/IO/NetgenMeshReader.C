//******************************************************************************
/*!
  \file      src/IO/NetgenMeshReader.C
  \author    J. Bakosi
  \date      Tue 17 Mar 2015 02:51:59 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Netgen mesh reader class definition
  \details   Netgen mesh reader class definition. Only supports tetrahedra.
*/
//******************************************************************************

#include <limits>
#include <array>
#include <cmath>

#include <NetgenMeshReader.h>

using tk::NetgenMeshReader;

void
NetgenMeshReader::read()
//******************************************************************************
//  Public interface for reading Netgen mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Read nodes
  readNodes();
  // Read elements
  readElements();

  // Clear failbit triggered by eof, so close() won't throw a false FAILED_CLOSE
  m_inFile.clear();
}

void
NetgenMeshReader::readNodes()
//******************************************************************************
//  Read nodes
//! \author J. Bakosi
//******************************************************************************
{
  int nnode;
  m_inFile >> nnode;
  ErrChk( nnode > 0,
          "Number of nodes must be greater than zero in file " + m_filename  );

  // Read in node coordinates: x-coord y-coord z-coord
  for (int i=0; i<nnode; ++i) {
    tk::real x, y, z;
    m_inFile >> x >> y >> z;
    m_mesh.x().push_back( x );
    m_mesh.y().push_back( y );
    m_mesh.z().push_back( z );
  }

  std::string s;
  getline( m_inFile, s );  // finish reading the last line
}

void
NetgenMeshReader::readElements()
//******************************************************************************
//  Read element connectivity
//! \author J. Bakosi
//******************************************************************************
{
  std::string s;
  int nel;

  // Read in number of tetrahedra
  m_inFile >> nel;
  ErrChk( nel > 0,
          "Number of tetrahedra (volume elements) must be greater than zero "
          "in file " + m_filename );
  getline( m_inFile, s );  // finish reading the last line

  // Read in tetrahedra element tags and connectivity
  for (int i=0; i<nel; ++i) {
    int tag;
    std::array< int, 4 > n;
    // tag n[1-4]
    m_inFile >> tag >> n[3] >> n[0] >> n[1] >> n[2];
    m_mesh.tettag().push_back( { tag } );
    m_mesh.tetinpoel().push_back( n[0] );
    m_mesh.tetinpoel().push_back( n[1] );
    m_mesh.tetinpoel().push_back( n[2] );
    m_mesh.tetinpoel().push_back( n[3] );
  }

  // Read in number of triangles
  m_inFile >> nel;
  ErrChk( nel > 0,
          "Number of triangles (surface elements) must be greater than zero in "
          "file " + m_filename );
  getline( m_inFile, s );  // finish reading the last line

  // Read in triangle element tags and connectivity
  for (int i=0; i<nel; ++i) {
    int tag;
    std::array< int, 3 > n;
    // tag n[1-3]
    m_inFile >> tag >> n[0] >> n[1] >> n[2];
    m_mesh.tritag().push_back( { tag } );
    m_mesh.triinpoel().push_back( n[0] );
    m_mesh.triinpoel().push_back( n[1] );
    m_mesh.triinpoel().push_back( n[2] );
  }
}
