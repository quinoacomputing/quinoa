// *****************************************************************************
/*!
  \file      src/IO/MeshFactory.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unstructured mesh reader and writer factory
  \details   Unstructured mesh reader and writer factory.
*/
// *****************************************************************************

#include <string>

#include "MeshFactory.hpp"
#include "MeshDetect.hpp"
#include "Timer.hpp"
#include "Reader.hpp"
#include "GmshMeshReader.hpp"
#include "NetgenMeshReader.hpp"
#include "ExodusIIMeshReader.hpp"
#include "HyperMeshReader.hpp"
#include "UGRIDMeshReader.hpp"
#include "ASCMeshReader.hpp"
#include "NetgenMeshWriter.hpp"
#include "GmshMeshWriter.hpp"
#include "ExodusIIMeshWriter.hpp"
#include "DerivedData.hpp"
#include "Reorder.hpp"
#include "QuinoaConfig.hpp"

#ifdef HAS_OMEGA_H
  #include "Omega_h_MeshReader.hpp"
#endif

namespace tk {

UnsMesh
readUnsMesh( const tk::Print& print,
             const std::string& filename,
             std::pair< std::string, tk::real >& timestamp )
// *****************************************************************************
//  Read unstructured mesh from file
//! \param[in] print Pretty printer
//! \param[in] filename Filename to read mesh from
//! \param[out] timestamp A time stamp consisting of a timer label (a string),
//!   and a time state (a tk::real in seconds) measuring the mesh read time
//! \return Unstructured mesh object
// *****************************************************************************
{
  print.diagstart( "Reading mesh from file ..." );

  // Read in mesh
  tk::Timer t;
 
  //! Create unstructured mesh to store mesh
  UnsMesh mesh;

  const auto meshtype = detectInput( filename );

  if (meshtype == MeshReaderType::GMSH)
    GmshMeshReader( filename ).readMesh( mesh );
  else if (meshtype == MeshReaderType::NETGEN)
    NetgenMeshReader( filename ).readMesh( mesh );
  else if (meshtype == MeshReaderType::EXODUSII)
    ExodusIIMeshReader( filename ).readMesh( mesh );
  else if (meshtype == MeshReaderType::ASC)
    ASCMeshReader( filename ).readMesh( mesh );
  else if (meshtype == MeshReaderType::UGRID)
    UGRIDMeshReader( filename ).readMesh( mesh );
  else if (meshtype == MeshReaderType::HYPER)
    HyperMeshReader( filename ).readMesh( mesh );

  timestamp = std::make_pair( "Read mesh from file", t.dsec() );

  print.diagend( "done" );

  // Return (move out) mesh object
  return mesh;
}

std::vector< std::pair< std::string, tk::real > >
writeUnsMesh( const tk::Print& print,
              const std::string& filename,
              UnsMesh& mesh,
              bool reorder )
// *****************************************************************************
//  Write unstructured mesh to file
//! \param[in] print Pretty printer
//! \param[in] filename Filename to write mesh to
//! \param[in] mesh Unstructured mesh object to write from
//! \param[in] reorder Whether to also reorder mesh nodes
//! \return Vector of time stamps consisting of a timer label (a string), and a
//!   time state (a tk::real in seconds) measuring the renumber and the mesh
//!   write time
// *****************************************************************************
{
  std::vector< std::pair< std::string, tk::real > > times;

  tk::Timer t;

  // If mesh has tetrahedra but no triangles, generate triangle connectivity
  if (!mesh.tetinpoel().empty() && mesh.triinpoel().empty()) {
    print.diagstart( "Generating missing surface mesh ..." );

    const auto& inpoel = mesh.tetinpoel();        // get tet connectivity
    auto esup = tk::genEsup( inpoel, 4 );         // elements surrounding points
    auto esuel = tk::genEsuelTet( inpoel, esup ); // elems surrounding elements
    auto& triinpoel = mesh.triinpoel();
    // collect boundary faces
    for (std::size_t e=0; e<esuel.size()/4; ++e) {
      auto mark = e*4;
      for (std::size_t f=0; f<4; ++f) {
        if (esuel[mark+f] == -1) {
          // extract triangle element connectivity from tetrahedron
          triinpoel.push_back( inpoel[ mark+tk::lpofa[f][0] ] );
          triinpoel.push_back( inpoel[ mark+tk::lpofa[f][1] ] );
          triinpoel.push_back( inpoel[ mark+tk::lpofa[f][2] ] );
        }
      }
    }

    print.diagend( "done" );
    times.emplace_back( "Generate surface mesh", t.dsec() );
    t.zero();
  }

  if (reorder) {
    print.diagstart( "Reordering mesh nodes ..." );

    // If mesh has tetrahedra elements, reorder based on those
    if (!mesh.tetinpoel().empty()) {

      auto& inpoel = mesh.tetinpoel();
      const auto psup = tk::genPsup( inpoel, 4, tk::genEsup( inpoel, 4 ) );
      auto map = tk::renumber( psup );
      tk::remap( inpoel, map );
      tk::remap( mesh.triinpoel(), map );
      tk::remap( mesh.x(), map );
      tk::remap( mesh.y(), map );
      tk::remap( mesh.z(), map );

    // If mesh has no tetrahedra elements, reorder based on triangle mesh if any
    } else if (!mesh.triinpoel().empty()) {

      auto& inpoel = mesh.triinpoel();
      const auto psup = tk::genPsup( inpoel, 3, tk::genEsup( inpoel, 3 ) );
      auto map = tk::renumber( psup );
      tk::remap( inpoel, map );
      tk::remap( mesh.x(), map );
      tk::remap( mesh.y(), map );
      tk::remap( mesh.z(), map );
    }

    print.diagend( "done" );
    times.emplace_back( "Reorder mesh", t.dsec() );
    t.zero();
  }

  print.diagstart( "Writing mesh to file ..." );

  const auto meshtype = pickOutput( filename );

  if (meshtype == MeshWriterType::GMSH)
    GmshMeshWriter( filename ).writeMesh( mesh );
  else if (meshtype == MeshWriterType::NETGEN)
    NetgenMeshWriter( filename ).writeMesh( mesh );
  else if (meshtype== MeshWriterType::EXODUSII)
    ExodusIIMeshWriter( filename, ExoWriter::CREATE ).writeMesh( mesh );

  print.diagend( "done" );
  times.emplace_back( "Write mesh to file", t.dsec() );

  return times;
}

} // tk::
