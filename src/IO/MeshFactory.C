// *****************************************************************************
/*!
  \file      src/IO/MeshFactory.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unstructured mesh reader and writer factory
  \details   Unstructured mesh reader and writer factory.
*/
// *****************************************************************************

#include <string>
#include <stdexcept>

#include "MeshFactory.h"
#include "Exception.h"
#include "Timer.h"
#include "Reader.h"
#include "GmshMeshReader.h"
#include "NetgenMeshReader.h"
#include "ExodusIIMeshReader.h"
#include "HyperMeshReader.h"
#include "NetgenMeshWriter.h"
#include "GmshMeshWriter.h"
#include "ExodusIIMeshWriter.h"
#include "DerivedData.h"
#include "Reorder.h"

namespace tk {

MeshReader
detectInput( const std::string& filename )
// *****************************************************************************
//  Detect input mesh file type
//! \param[in] filename File to open and detect its type
//! \return enum specifying the mesh reader type
//! \author J. Bakosi
// *****************************************************************************
{
  // Get first three letters from input file
  std::string s( Reader( filename ).firstline().substr(0,4) );

  if ( s.find("$Me") != std::string::npos ) {
    return MeshReader::GMSH;
  } else if ( s.find("CDF") != std::string::npos ||
              s.find("HDF") != std::string::npos ) {
    return MeshReader::EXODUSII;
  } else if ( s.find("<?x") != std::string::npos ) {
    return MeshReader::HYPERMESH;
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
    return MeshReader::NETGEN;
  }
}

MeshWriter
pickOutput( const std::string& filename )
// *****************************************************************************
//  Determine output mesh file type
//! \param[in] filename Filename to pick its type based on extension given
//! \return enum specifying the mesh writer type
//! \author J. Bakosi
// *****************************************************************************
{
  // Get extension of input file name
  std::string fn = filename;
  std::string ext( fn.substr(fn.find_last_of(".") + 1) );

  if ( ext == "msh" ) {
    return MeshWriter::GMSH;
  } else if ( ext == "exo" || ext == "h5" ) {
    return MeshWriter::EXODUSII;
  } else if ( ext == "mesh" ) {
    return MeshWriter::NETGEN;
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
//! \author J. Bakosi
// *****************************************************************************
{
  print.diagstart( "Reading mesh from file ..." );

  // Read in mesh
  tk::Timer t;
 
  //! Create unstructured mesh to store mesh
  UnsMesh mesh;

  const auto meshtype = detectInput( filename );

  if (meshtype == MeshReader::GMSH)
    GmshMeshReader( filename ).readMesh( mesh );
  else if (meshtype == MeshReader::NETGEN)
    NetgenMeshReader( filename ).readMesh( mesh );
  else if (meshtype == MeshReader::EXODUSII)
    ExodusIIMeshReader( filename ).readMesh( mesh );
  else if (meshtype == MeshReader::HYPERMESH)
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
//! \author J. Bakosi
// *****************************************************************************
{
  std::vector< std::pair< std::string, tk::real > > times;

  tk::Timer t;

  if (reorder) {
    print.diagstart( "Reordering mesh nodes ..." );

    auto& inpoel = mesh.tetinpoel();
    const auto psup = tk::genPsup( inpoel, 4, tk::genEsup( inpoel, 4 ) );
    std::vector< std::size_t > map, invmap;
    std::tie( map, invmap ) = tk::renumber( psup );
    tk::remap( inpoel, map );

    print.diagend( "done" );
    times.emplace_back( "Renumber mesh", t.dsec() );
    t.zero();
  }

  print.diagstart( "Writing mesh to file ..." );

  const auto meshtype = pickOutput( filename );

  if (meshtype == MeshWriter::GMSH)
    GmshMeshWriter( filename ).writeMesh( mesh );
  else if (meshtype == MeshWriter::NETGEN)
    NetgenMeshWriter( filename ).writeMesh( mesh );
  else if (meshtype== MeshWriter::EXODUSII)
    ExodusIIMeshWriter( filename, ExoWriter::CREATE ).writeMesh( mesh );

  print.diagend( "done" );
  times.emplace_back( "Write mesh to file", t.dsec() );

  return times;
}

} // tk::
