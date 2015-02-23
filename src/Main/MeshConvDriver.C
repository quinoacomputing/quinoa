//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.C
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 08:41:01 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Mesh converter driver
  \details   Mesh converter driver.
*/
//******************************************************************************

#include <Print.h>
#include <Exception.h>
#include <Factory.h>
#include <MeshConvDriver.h>
#include <MeshConv/CmdLine/Parser.h>
#include <GmshMeshReader.h>
#include <NetgenMeshReader.h>
#include <ExodusIIMeshReader.h>
#include <NetgenMeshWriter.h>
#include <GmshMeshWriter.h>
#include <ExodusIIMeshWriter.h>
#include <meshconv.decl.h>

using meshconv::MeshConvDriver;

MeshConvDriver::MeshConvDriver( const tk::Print& print,
                                const ctr::CmdLine& cmdline )
//******************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
//! \author J. Bakosi
//******************************************************************************
{
  // Save input file name
  m_input = cmdline.get< tag::io, tag::input >();
  // Save output file name
  m_output = cmdline.get< tag::io, tag::output >();
}

void
MeshConvDriver::execute() const
//******************************************************************************
//  Execute: Convert Gmsh mesh to Exodus II mesh
//! \author J. Bakosi
//******************************************************************************
{
  //! Mesh readers factory
  std::map< MeshReaderType, std::function<tk::Reader*()> > readers;

  //! Mesh writers factory
  std::map< MeshWriterType, std::function<tk::Writer*()> > writers;

  //! Create unstructured mesh to store mesh
  tk::UnsMesh mesh;

  // Register mesh readers
  tk::record< tk::GmshMeshReader >( readers, MeshReaderType::GMSH,
                                    m_input, std::ref(mesh) );
  tk::record< tk::NetgenMeshReader >( readers, MeshReaderType::NETGEN,
                                      m_input, std::ref(mesh) );
  tk::record< tk::ExodusIIMeshReader >( readers, MeshReaderType::EXODUSII,
                                        m_input, std::ref(mesh) );

  // Register mesh writers
  tk::record< tk::GmshMeshWriter >( writers, MeshWriterType::GMSH,
                                    m_output, std::ref(mesh) );
  tk::record< tk::NetgenMeshWriter >( writers, MeshWriterType::NETGEN,
                                      m_output, std::ref(mesh) );
  tk::record< tk::ExodusIIMeshWriter >( writers, MeshWriterType::EXODUSII,
                                        m_output, std::ref(mesh) );

  // Read in mesh
  tk::instantiate( readers, detectInput() )->read();

  // Write out mesh
  tk::instantiate( writers, pickOutput() )->write();
}

MeshConvDriver::MeshReaderType
MeshConvDriver::detectInput() const
//******************************************************************************
//  Detect input mesh file type
//! \return enum specifying the mesh reader type
//! \author J. Bakosi
//******************************************************************************
{
  // Get first three letters from input file
  std::string s( tk::Reader( m_input ).firstline().substr(0,3) );

  if ( s == "$Me" ) {

    return MeshReaderType::GMSH;

  } else if ( s == "CDF" || s == "HDF" ) {

    return MeshReaderType::EXODUSII;

  } else {

    try {

      std::stoi(s);    // try to convert to an integer

    } catch ( std::invalid_argument ) {

      Throw( "Input mesh file type could not be determined from header: " +
             m_input );

    }

    // could also catch std::out_of_range, the other exception potentially
    // thrown by std::stoi(), but a three-digit integer will always fit into int

    // if we got here, the above string-to-integer conversion succeeded
    return MeshReaderType::NETGEN;

  }
}

MeshConvDriver::MeshWriterType
MeshConvDriver::pickOutput() const
//******************************************************************************
//  Determine output mesh file type
//! \return enum specifying the mesh writer type
//! \author J. Bakosi
//******************************************************************************
{
  // Get extension of input file name
  std::string fn = m_output;
  std::string ext( fn.substr(fn.find_last_of(".") + 1) );

  if ( ext == "msh" ) {

    return MeshWriterType::GMSH;

  } else if ( ext == "exo" || ext == "h5" ) {

    return MeshWriterType::EXODUSII;

  } else if ( ext == "mesh" ) {

    return MeshWriterType::NETGEN;

  } else {

    Throw( "Output mesh file type could not be determined from extension of "
           "filename '" + m_output + "'; valid extensions are: "
           "'msh' for Gmsh, 'exo' or 'h5' for ExodusII, 'mesh' for Netgen's "
           "neutral" );

  }
}
