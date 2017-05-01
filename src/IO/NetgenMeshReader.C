// *****************************************************************************
/*!
  \file      src/IO/NetgenMeshReader.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Netgen mesh reader class definition
  \details   Netgen mesh reader class definition. Only supports tetrahedra.
*/
// *****************************************************************************

#include <array>
#include <istream>
#include <string>
#include <vector>
#include <cstddef>

#include "Types.h"
#include "Exception.h"
#include "UnsMesh.h"
#include "Reorder.h"
#include "NetgenMeshReader.h"

using tk::NetgenMeshReader;

void
NetgenMeshReader::readMesh( UnsMesh& mesh )
// *****************************************************************************
//  Read Netgen mesh
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  // Read nodes
  readNodes( mesh );
  // Read elements
  readElements( mesh );
}

void
NetgenMeshReader::readNodes( UnsMesh& mesh )
// *****************************************************************************
//  Read nodes
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  int nnode;
  m_inFile >> nnode;
  ErrChk( nnode > 0,
          "Number of nodes must be greater than zero in file " + m_filename  );

  // Read in node coordinates: x-coord y-coord z-coord
  for (int i=0; i<nnode; ++i) {
    tk::real x, y, z;
    m_inFile >> x >> y >> z;
    mesh.x().push_back( x );
    mesh.y().push_back( y );
    mesh.z().push_back( z );
  }

  std::string s;
  getline( m_inFile, s );  // finish reading the last line
}

void
NetgenMeshReader::readElements( UnsMesh& mesh )
// *****************************************************************************
//  Read element connectivity
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  int nel;

  // Read in number of tetrahedra
  m_inFile >> nel;
  if (!m_inFile.eof()) {
    ErrChk( nel > 0, "Number of tetrahedra (volume elements) must be greater "
                     "than zero in file " + m_filename );
    std::string s;
    getline( m_inFile, s );  // finish reading the last line

    // Read in tetrahedra element tags and connectivity
    for (int i=0; i<nel; ++i) {
      int tag;
      std::array< std::size_t, 4 > n;
      // tag n[1-4]
      m_inFile >> tag >> n[3] >> n[0] >> n[1] >> n[2];
      mesh.tetinpoel().push_back( n[0] );
      mesh.tetinpoel().push_back( n[1] );
      mesh.tetinpoel().push_back( n[2] );
      mesh.tetinpoel().push_back( n[3] );
    }

    // Shift node IDs to start from zero
    shiftToZero( mesh.tetinpoel() );
  }

  // Read in number of triangles
  m_inFile >> nel;
  if (!m_inFile.eof()) {
    ErrChk( nel > 0, "Number of triangles (surface elements) must be greater "
                     "than zero in file " + m_filename );
    std::string s;
    getline( m_inFile, s );  // finish reading the last line

    // Read in triangle element tags and connectivity
    for (int i=0; i<nel; ++i) {
      int tag;
      std::array< std::size_t, 3 > n;
      // tag n[1-3]
      m_inFile >> tag >> n[0] >> n[1] >> n[2];
      mesh.triinpoel().push_back( n[0] );
      mesh.triinpoel().push_back( n[1] );
      mesh.triinpoel().push_back( n[2] );
    }

    // Shift node IDs to start from zero
    shiftToZero( mesh.triinpoel() );
  }
}
