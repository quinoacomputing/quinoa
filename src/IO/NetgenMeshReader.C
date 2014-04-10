//******************************************************************************
/*!
  \file      src/IO/NetgenMeshReader.C
  \author    J. Bakosi
  \date      Wed 09 Apr 2014 07:20:12 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Netgen mesh reader class definition
  \details   Netgen mesh reader class definition
*/
//******************************************************************************

#include <limits>
#include <cmath>

#include <GmshMesh.h>
#include <NetgenMeshReader.h>

using quinoa::NetgenMeshReader;

void
NetgenMeshReader::read()
//******************************************************************************
//  Public interface for read Netgen mesh
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
//  Read "$Nodes--$EndNodes" section
//! \author J. Bakosi
//******************************************************************************
{
  std::size_t nnode;
  m_inFile >> nnode;
  ErrChk( nnode > 0, tk::ExceptType::FATAL,
          "Number of nodes must be greater than zero in file " + m_filename  );

  // Read in node coordinates: x-coord y-coord z-coord
  for ( std::size_t i=0; i<nnode; ++i ) {
    tk::point coord;

    m_inFile >> coord[0] >> coord[1] >> coord[2];

    m_mesh.nodeId().push_back( i );
    m_mesh.coord().push_back( coord );
  }

  std::string s;
  getline( m_inFile, s );  // finish reading the last line
}

void
NetgenMeshReader::readElements()
//******************************************************************************
//  Read "$Elements--$EndElements" section
//! \author J. Bakosi
//******************************************************************************
{
  std::string s;

  // Read in number of elements
  int nel;
  m_inFile >> nel;
  ErrChk( nel > 0, tk::ExceptType::FATAL,
          "Number of elements must be greater than zero in file " +
          m_filename );

  getline( m_inFile, s );  // finish reading the last line

  for (int i=0; i<nel; ++i) {
    std::vector< int > n( 4, 0 );
    // mat n[1-4], throw away mat
    m_inFile >> n[3] >> n[3] >> n[0] >> n[1] >> n[2];
    n[0] -= 1;
    n[1] -= 1;
    n[2] -= 1;
    n[3] -= 1;
    m_mesh.tetrahedronId().push_back( i );
    m_mesh.tettag().push_back( {1} );
    m_mesh.tetinpoel().push_back( n );
  }
}
