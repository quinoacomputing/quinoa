//******************************************************************************
/*!
  \file      src/IO/MeshFactory.C
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:25:08 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unstructured mesh reader and writer factory
  \details   Unstructured mesh reader and writer factory.
*/
//******************************************************************************

#include <string>
#include <stdexcept>

#include "MeshFactory.h"
#include "Exception.h"
#include "Timer.h"
#include "Reader.h"
#include "GmshMeshReader.h"
#include "NetgenMeshReader.h"
#include "ExodusIIMeshReader.h"
#include "NetgenMeshWriter.h"
#include "GmshMeshWriter.h"
#include "ExodusIIMeshWriter.h"

namespace tk {

MeshReader
detectInput( const std::string& filename )
//******************************************************************************
//  Detect input mesh file type
//! \param[in] filename File to open and detect its type
//! \return enum specifying the mesh reader type
//! \author J. Bakosi
//******************************************************************************
{
  // Get first three letters from input file
  std::string s( Reader( filename ).firstline().substr(0,3) );

  if ( s == "$Me" ) {
    return MeshReader::GMSH;
  } else if ( s == "CDF" || s == "HDF" ) {
    return MeshReader::EXODUSII;
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
//******************************************************************************
//  Determine output mesh file type
//! \param[in] filename Filename to pick its type based on extension given
//! \return enum specifying the mesh writer type
//! \author J. Bakosi
//******************************************************************************
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
readUnsMesh( const std::string& filename,
             std::pair< std::string, tk::real >& timestamp )
//******************************************************************************
//  Read unstructured mesh from file
//! \param[in] filename Filename to read mesh from
//! \param[out] timestamp A time stamp consisting of a timer label (a string),
//!   and a time state (a tk::real in seconds) measuring the mesh read time
//! \return Unstructured mesh object
//! \author J. Bakosi
//******************************************************************************
{
  // Read in mesh
  tk::Timer t;
 
  //! Create unstructured mesh to store mesh
  UnsMesh mesh;

  const auto meshtype = detectInput( filename );

  if (meshtype == MeshReader::GMSH)
    GmshMeshReader( filename ).readMesh( mesh );
  else if (meshtype == MeshReader::NETGEN)
    NetgenMeshReader( filename ).readMesh( mesh );
  else if (meshtype== MeshReader::EXODUSII)
    ExodusIIMeshReader( filename ).readMesh( mesh );

  timestamp = std::make_pair( "Read mesh from file", t.dsec() );

  // Return (move out) mesh object
  return mesh;
}

void
writeUnsMesh( const std::string& filename,
              const UnsMesh& mesh,
              std::pair< std::string, tk::real >& timestamp )
//******************************************************************************
//  Write unstructured mesh to file
//! \param[in] filename Filename to write mesh to
//! \param[in] mesh Unstructured mesh object to write from
//! \param[out] timestamp A time stamp consisting of a timer label (a string),
//!   and a time state (a tk::real in seconds) measuring the mesh write time
//! \author J. Bakosi
//******************************************************************************
{
  tk::Timer t;

  const auto meshtype = pickOutput( filename );

  if (meshtype == MeshWriter::GMSH)
    GmshMeshWriter( filename ).writeMesh( mesh );
  else if (meshtype == MeshWriter::NETGEN)
    NetgenMeshWriter( filename ).writeMesh( mesh );
  else if (meshtype== MeshWriter::EXODUSII)
    ExodusIIMeshWriter( filename, ExoWriter::CREATE ).writeMesh( mesh );

  timestamp = std::make_pair( "Write mesh to file", t.dsec() );
}

} // tk::
