// *****************************************************************************
/*!
  \file      src/IO/HyperMeshReader.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Hyper mesh reader class definition
  \details   Hyper mesh reader class definition. Only supports tetrahedra.
*/
// *****************************************************************************

#include <array>
#include <string>
#include <sstream>

#include "NoWarning/pugixml.h"

#include "Types.h"
#include "Exception.h"
#include "UnsMesh.h"
#include "HyperMeshReader.h"

using tk::HyperMeshReader;

void
HyperMeshReader::readMesh( UnsMesh& mesh )
// *****************************************************************************
//  Read Hyper mesh
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  auto filenames = getFileNames();

  // Read nodes
  readNodes( filenames.first, mesh );
  // Read elements
  readElements( filenames.second, mesh );
}


std::pair< std::string, std::string >
HyperMeshReader::getFileNames() const
// *****************************************************************************
//  Read Hyper mesh metadata and extract filenames we need to read
//! \return Vector of strings containing the filenames
//! \author J. Bakosi
// *****************************************************************************
{
  // Read XML file for metadata
  pugi::xml_document meta;
  auto res = meta.load_file( m_filename.c_str() );
  ErrChk( res == true, "Could not parse XML attempting to read HyperMesh "
          "metadata file '" + m_filename + "', at character position " +
          std::to_string(res.offset) + " : " + res.description() );

  // Extract path from metadata filename
  std::string fn = m_filename;
  std::string path( fn.substr(0,fn.find_last_of("/")+1) );
  // Prepend path to filenames to be extracted
  std::pair< std::string, std::string > filenames( path, path );
 
  for (const auto& child : meta.children("mesh"))
    for (const auto& group : child.children()) {

      if (group.name() == std::string("coordinates")) {

        filenames.first += group.attribute("file").value();

      } else if (group.name() == std::string("element_set")) {

        filenames.second += group.attribute("file").value();
        Assert( group.attribute("topology").value() ==
                  std::string("four_node_tet"), "Only pure tetrahedron-element "
                  "meshes are supported" );
      }
    }

  return filenames;
}
        
void
HyperMeshReader::readNodes( const std::string& filename, UnsMesh& mesh ) const
// *****************************************************************************
//  Read nodes
//! \param[in] filename Filename to read nodes from
//! \param[in] mesh Unstructured mesh object to put nodes coordinates
//! \note We throw away the node ID, which means the nodes must be in order.
//! \author J. Bakosi
// *****************************************************************************
{
  // Read in node coordinates: x-coord y-coord z-coord
  for (auto& l : tk::Reader(filename).lines()) {
    std::stringstream ss(l);
    int id;
    tk::real x, y, z;
    ss >> id >> x >> y >> z;
    mesh.x().push_back( x );
    mesh.y().push_back( y );
    mesh.z().push_back( z );
  }
}

void
HyperMeshReader::readElements( const std::string& filename, UnsMesh& mesh )
const
// *****************************************************************************
//  Read element connectivity
//! \param[in] filename Filename to read nodes from
//! \param[in] mesh Unstructured mesh object to put element connectivity
//! \note We throw away the element ID.
//! \author J. Bakosi
// *****************************************************************************
{
  for (auto& l : tk::Reader(filename).lines()) {
    std::stringstream ss(l);
    int id;
    std::array< std::size_t, 4 > n;
    ss >> id >> n[0] >> n[1] >> n[2] >> n[3];
    mesh.tetinpoel().push_back( n[0] );
    mesh.tetinpoel().push_back( n[1] );
    mesh.tetinpoel().push_back( n[2] );
    mesh.tetinpoel().push_back( n[3] );
  }
}
