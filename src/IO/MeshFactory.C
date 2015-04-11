//******************************************************************************
/*!
  \file      src/IO/MeshFactory.C
  \author    J. Bakosi
  \date      Sat 11 Apr 2015 06:34:20 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unstructured mesh reader and writer factory
  \details   Unstructured mesh reader and writer factory.
*/
//******************************************************************************

#include <stdexcept>

#include <MeshFactory.h>
#include <Factory.h>
#include <Reader.h>
#include <GmshMeshReader.h>
#include <NetgenMeshReader.h>
#include <ExodusIIMeshReader.h>
#include <NetgenMeshWriter.h>
#include <GmshMeshWriter.h>
#include <ExodusIIMeshWriter.h>

namespace tk {

MeshReaderType
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
    return MeshReaderType::GMSH;
  } else if ( s == "CDF" || s == "HDF" ) {
    return MeshReaderType::EXODUSII;
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
  //! Mesh readers factory
  std::map< MeshReaderType, std::function<Reader*()> > readers;

  //! Create unstructured mesh to store mesh
  UnsMesh mesh;

  // Register mesh readers
  record< GmshMeshReader >
        ( readers, MeshReaderType::GMSH, filename, std::ref(mesh) );
  record< NetgenMeshReader >
        ( readers, MeshReaderType::NETGEN, filename, std::ref(mesh) );
  record< ExodusIIMeshReader >
        ( readers, MeshReaderType::EXODUSII, filename, std::ref(mesh) );

  // Read in mesh
  tk::Timer t;
  instantiate( readers, detectInput( filename ) )->read();
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
  //! Mesh writers factory
  std::map< MeshWriterType, std::function<Writer*()> > writers;

  // Register mesh writers
  record< GmshMeshWriter >
        ( writers, MeshWriterType::GMSH, filename, std::ref(mesh) );
  record< NetgenMeshWriter >
        ( writers, MeshWriterType::NETGEN, filename, std::ref(mesh) );
  record< ExodusIIMeshWriter >
        ( writers, MeshWriterType::EXODUSII, filename, std::ref(mesh) );

  // Write out mesh
  tk::Timer t;
  instantiate( writers, pickOutput( filename ) )->write();
  timestamp = std::make_pair( "Write mesh to file", t.dsec() );
}

} // tk::
