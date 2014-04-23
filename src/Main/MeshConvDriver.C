//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.C
  \author    J. Bakosi
  \date      Wed Apr 23 11:41:13 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshConvDriver that drives MeshConv
  \details   MeshConvDriver that drives MeshConv
*/
//******************************************************************************

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

using meshconv::MeshConvDriver;

MeshConvDriver::MeshConvDriver(int argc, char** argv, const tk::Print& print)
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \param[in] print     Simple pretty printer
//! \author J. Bakosi
//******************************************************************************
{
  // Parse command line into cmdline
  CmdLineParser cmdParser(argc, argv, print, m_cmdline);

  // Register mesh readers
  tk::record< quinoa::GmshMeshReader >( m_readers, MeshReaderType::GMSH,
    m_cmdline->get<tag::io, tag::input>(), std::ref(m_mesh) );
  tk::record< quinoa::NetgenMeshReader >( m_readers, MeshReaderType::NETGEN,
    m_cmdline->get<tag::io, tag::input>(), std::ref(m_mesh) );
  tk::record< quinoa::ExodusIIMeshReader >( m_readers, MeshReaderType::EXODUSII,
    m_cmdline->get<tag::io, tag::input>(), std::ref(m_mesh) );

  // Register mesh writers
  tk::record< quinoa::GmshMeshWriter >( m_writers, MeshWriterType::GMSH,
    m_cmdline->get<tag::io, tag::output>(), std::ref(m_mesh) );
  tk::record< quinoa::NetgenMeshWriter >( m_writers, MeshWriterType::NETGEN,
    m_cmdline->get<tag::io, tag::output>(), std::ref(m_mesh) );
  tk::record< quinoa::ExodusIIMeshWriter >( m_writers, MeshWriterType::EXODUSII,
    m_cmdline->get<tag::io, tag::output>(), std::ref(m_mesh) );
}

void
MeshConvDriver::execute() const
//******************************************************************************
//  Execute: Convert Gmsh mesh to Exodus II mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Create mesh reader
  std::unique_ptr< tk::Reader > m_reader =
    tk::instantiate( m_readers, detectInput() );
  // Create mesh writer
  std::unique_ptr< tk::Writer > m_writer =
    tk::instantiate( m_writers, pickOutput() );

  // Read in mesh
  m_reader->read();
  // Write out mesh
  m_writer->write();
}

meshconv::MeshReaderType
MeshConvDriver::detectInput() const
//******************************************************************************
//  Detect input mesh file type
//! \author J. Bakosi
//******************************************************************************
{
  // Get first three letters from input file
  std::string s(
    tk::Reader(m_cmdline->get<tag::io,tag::input>()).firstline().substr(0,3) );

  if ( s == "$Me" ) {

    return MeshReaderType::GMSH;

  } else if ( s == "CDF" || s == "HDF" ) {

    return MeshReaderType::EXODUSII;

  } else {

    try {

      std::stoi(s);    // try to convert to an integer

    } catch ( std::invalid_argument ) {

      Throw( tk::ExceptType::RUNTIME,
             "Input mesh file type could not be determined from header: " +
             m_cmdline->get<tag::io,tag::input>() );

    }

    // could also catch std::out_of_range, the other exception potentially
    // thrown by std::stoi(), but a three-digit integer will always fit into int

    // if we got here, the above string-to-integer conversion succeeded
    return MeshReaderType::NETGEN;

  }
}

meshconv::MeshWriterType
MeshConvDriver::pickOutput() const
//******************************************************************************
//  Determine output mesh file type
//! \author J. Bakosi
//******************************************************************************
{
  // Get extension of input file name
  std::string fn = m_cmdline->get< tag::io, tag::output >();
  std::string ext( fn.substr(fn.find_last_of(".") + 1) );

  if ( ext == "msh" ) {

    return MeshWriterType::GMSH;

  } else if ( ext == "exo" || ext == "h5" ) {

    return MeshWriterType::EXODUSII;

  } else if ( ext == "mesh" ) {

    return MeshWriterType::NETGEN;

  } else {

    Throw( tk::ExceptType::RUNTIME,
           "Output mesh file type could not be determined from extension of "
           "filename '" +
           m_cmdline->get<tag::io,tag::output>() + "'; valid extensions are: "
           "'msh' for Gmsh, 'exo' or 'h5' for ExodusII, 'mesh' for Netgen's "
           "neutral" );

  }
}
