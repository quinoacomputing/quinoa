// *****************************************************************************
/*!
  \file      src/IO/MeshFactory.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Unstructured mesh reader and writer factory
  \details   Unstructured mesh reader and writer factory.
*/
// *****************************************************************************

#include <string>
#include <stdexcept>
#include <fstream>

#include "MeshFactory.h"
#include "Exception.h"
#include "Timer.h"
#include "Reader.h"
#include "GmshMeshReader.h"
#include "NetgenMeshReader.h"
#include "ExodusIIMeshReader.h"
#include "HyperMeshReader.h"
#include "ASCMeshReader.h"
#include "Omega_h_MeshReader.h"
#include "NetgenMeshWriter.h"
#include "GmshMeshWriter.h"
#include "ExodusIIMeshWriter.h"
#include "DerivedData.h"
#include "Reorder.h"

namespace tk {

MeshReaderType
detectInput( const std::string& filename )
// *****************************************************************************
//  Detect input mesh file type
//! \param[in] filename File to open and detect its type
//! \return enum specifying the mesh reader type
// *****************************************************************************
{
  std::ifstream inFile;

  // Check if file exists, throw exception if it does not
  inFile.open( filename, std::ifstream::in );
  ErrChk( inFile.good(), "Failed to open file: " + filename );

  // Attempt to read a character, throw if it fails
  // It is curious that on some systems opening a directory instead of a file
  // with the above ifstream::open() call does not set the failbit. Thus we get
  // here fine, so we try to read a character from it. If the read fails, we
  // assume it is a directory and assume that it is in Omega_h osh file format.
  // Read more at: http://stackoverflow.com/questions/9591036/
  // ifstream-open-doesnt-set-error-bits-when-argument-is-a-directory.
  inFile.get();
  if (!inFile.good()) return MeshReaderType::OMEGA_H;

  // Close it
  inFile.close();
  ErrChk( !inFile.fail(), "Failed to close file: " + filename );

  // Get first three letters from input file
  std::string s( Reader( filename ).firstline().substr(0,4) );

  if ( s.find("$Me") != std::string::npos ) {
    return MeshReaderType::GMSH;
  } else if ( s.find("CDF") != std::string::npos ||
              s.find("HDF") != std::string::npos ) {
    return MeshReaderType::EXODUSII;
  } else if ( s.find("<?x") != std::string::npos ) {
    return MeshReaderType::HYPER;
  } else if ( s.find("*nd") != std::string::npos ) {
    return MeshReaderType::ASC;
  } else {
    try {
      std::stoi(s);    // try to convert to an integer
    } catch ( std::invalid_argument ) {
      Throw( "Input mesh file type could not be determined from header: " +
             filename );
    }
    // could also catch std::out_of_range, the other exception potentially
    // thrown by std::stoi(), but a three-digit integer will always fit into int

    // if we got here, the above string-to-integer conversion succeeded
    return MeshReaderType::NETGEN;
  }
}

MeshWriterType
pickOutput( const std::string& filename )
// *****************************************************************************
//  Determine output mesh file type
//! \param[in] filename Filename to pick its type based on extension given
//! \return enum specifying the mesh writer type
// *****************************************************************************
{
  // Get extension of input file name
  std::string fn = filename;
  std::string ext( fn.substr(fn.find_last_of(".") + 1) );

  if ( ext == "msh" ) {
    return MeshWriterType::GMSH;
  } else if ( ext == "exo" || ext == "h5" ) {
    return MeshWriterType::EXODUSII;
  } else if ( ext == "mesh" ) {
    return MeshWriterType::NETGEN;
  } else {
    Throw( "Output mesh file type could not be determined from extension of "
           "filename '" + filename + "'; valid extensions are: "
           "'msh' for Gmsh, 'exo' or 'h5' for ExodusII, 'mesh' for Netgen's "
           "neutral" );
  }
}

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
