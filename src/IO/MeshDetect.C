//******************************************************************************
/*!
  \file      src/IO/MeshDetect.C
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 03:37:15 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Mesh converter driver
  \details   Mesh converter driver.
*/
//******************************************************************************

#include <stdexcept>

#include <MeshDetect.h>
#include <Reader.h>

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

} // tk::
