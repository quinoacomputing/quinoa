// *****************************************************************************
/*!
  \file      src/IO/NetgenMeshWriter.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Netgen mesh writer class definition
  \details   Netgen mesh writer class definition. Only supports tetrahedra.
*/
// *****************************************************************************

#include <iomanip>
#include <ostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>
#include <cstring>
#include <cstddef>

#include "Types.h"
#include "UnsMesh.h"
#include "Exception.h"
#include "NetgenMeshWriter.h"

using tk::NetgenMeshWriter;

void
NetgenMeshWriter::writeMesh( const UnsMesh& mesh )
// *****************************************************************************
//  Public interface for writing Netgen mesh
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  writeNodes( mesh );
  writeElements( mesh );
}

void
NetgenMeshWriter::writeNodes( const UnsMesh& mesh )
// *****************************************************************************
//  Write nodes
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  // Write out number of nodes
  m_outFile << mesh.nnode() << std::endl;

  // Write node coordinates: x-coord y-coord z-coord
  m_outFile << std::setprecision(6) << std::fixed;
  for ( std::size_t i=0; i<mesh.nnode(); ++i ) {
    m_outFile << '\t' << mesh.x()[i]
              << '\t' << mesh.y()[i]
              << '\t' << mesh.z()[i] << std::endl;
  }
}

void
NetgenMeshWriter::writeElements( const UnsMesh& mesh )
// *****************************************************************************
//  Write elements, i.e., connectivity
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  if (mesh.tetinpoel().empty()) return;

  // Make sure tetrahedron element connectivity starts with zero
  Assert( *std::minmax_element( begin(mesh.tetinpoel()),
                                end(mesh.tetinpoel()) ).first == 0,
          "tetrahedron node ids should start from zero" );

  // Get number of tetrahedra in mesh
  auto n = mesh.tetinpoel().size()/4;  

  // Write out number of tetrahedra
  m_outFile << n << std::endl;

  // Create empty tag vector
  std::vector< std::vector< int > > tg;
  tg.resize( n );
  for (auto& t : tg) t.push_back( 0 );

  // Write out tetrehadra element tags and connectivity
  for (std::size_t i=0; i<n; ++i) {
    // tag n[1-4]
    m_outFile << '\t' << tg[i][0]
              << '\t' << mesh.tetinpoel()[i*4+3]+1
              << '\t' << mesh.tetinpoel()[i*4+0]+1
              << '\t' << mesh.tetinpoel()[i*4+1]+1
              << '\t' << mesh.tetinpoel()[i*4+2]+1 << std::endl;
  }

  if (mesh.triinpoel().empty()) return;

  // Make sure triangle element connectivity starts with zero
  Assert( *std::minmax_element( begin(mesh.triinpoel()),
                                end(mesh.triinpoel()) ).first == 0,
          "triangle node ids should start from zero" );

  // Get number of triangles in mesh
  n = mesh.triinpoel().size()/3;

  // Write out number of triangles
  m_outFile << n << std::endl;

  // Create empty tag vector if there is no tag
  tg.clear();
  tg.resize( n );
  for (auto& t : tg) t.push_back( 0 );

  // Write out triangle element tags and connectivity
  for (std::size_t i=0; i<n; ++i) {
    // tag n[1-4]
    m_outFile << '\t' << tg[i][0]
              << '\t' << mesh.triinpoel()[i*3+0]+1
              << '\t' << mesh.triinpoel()[i*3+1]+1
              << '\t' << mesh.triinpoel()[i*3+2]+1 << std::endl;
  }
}
