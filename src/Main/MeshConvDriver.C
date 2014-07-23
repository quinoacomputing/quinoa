//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.C
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 03:43:45 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     MeshConvDriver that drives MeshConv
  \details   MeshConvDriver that drives MeshConv
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
#include <Handler.h>

using meshconv::MeshConvDriver;

MeshConvDriver::MeshConvDriver( const tk::Print& print,
                                const ctr::CmdLine& cmdline )
//******************************************************************************
//  Constructor
//! \author J. Bakosi
//******************************************************************************
{
  try {

    // Save input file name
    m_input = cmdline.get< tag::io, tag::input >();
    // Save output file name
    m_output = cmdline.get< tag::io, tag::output >();

  } catch (...) { tk::processException(); }
}

void
MeshConvDriver::execute() const
//******************************************************************************
//  Execute: Convert Gmsh mesh to Exodus II mesh
//! \author J. Bakosi
//******************************************************************************
{
  try {

    //! Mesh readers factory
    std::map< MeshReaderType, std::function<tk::Reader*()> > readers;

    //! Mesh writers factory
    std::map< MeshWriterType, std::function<tk::Writer*()> > writers;

    //! Create unstructured mesh to store mesh
    quinoa::UnsMesh mesh;

    // Register mesh readers
    tk::record< quinoa::GmshMeshReader >( readers, MeshReaderType::GMSH,
                                          m_input, std::ref(mesh) );
    tk::record< quinoa::NetgenMeshReader >( readers, MeshReaderType::NETGEN,
                                            m_input, std::ref(mesh) );
    tk::record< quinoa::ExodusIIMeshReader >( readers, MeshReaderType::EXODUSII,
                                              m_input, std::ref(mesh) );

    // Register mesh writers
    tk::record< quinoa::GmshMeshWriter >( writers, MeshWriterType::GMSH,
                                          m_output, std::ref(mesh) );
    tk::record< quinoa::NetgenMeshWriter >( writers, MeshWriterType::NETGEN,
                                            m_output, std::ref(mesh) );
    tk::record< quinoa::ExodusIIMeshWriter >( writers, MeshWriterType::EXODUSII,
                                              m_output, std::ref(mesh) );

    // Read in mesh
    tk::instantiate( readers, detectInput() )->read();

    // Write out mesh
    tk::instantiate( writers, pickOutput() )->write();

  } catch (...) { tk::processException(); }
}

MeshConvDriver::MeshReaderType
MeshConvDriver::detectInput() const
//******************************************************************************
//  Detect input mesh file type
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
