//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.C
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 02:59:23 PM MST
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
#include <MeshDetect.h>
#include <GmshMeshReader.h>
#include <NetgenMeshReader.h>
#include <ExodusIIMeshReader.h>
#include <NetgenMeshWriter.h>
#include <GmshMeshWriter.h>
#include <ExodusIIMeshWriter.h>

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
  std::map< tk::MeshReaderType, std::function<tk::Reader*()> > readers;

  //! Mesh writers factory
  std::map< tk::MeshWriterType, std::function<tk::Writer*()> > writers;

  //! Create unstructured mesh to store mesh
  tk::UnsMesh mesh;

  // Register mesh readers
  tk::record< tk::GmshMeshReader >( readers, tk::MeshReaderType::GMSH,
                                    m_input, std::ref(mesh) );
  tk::record< tk::NetgenMeshReader >( readers, tk::MeshReaderType::NETGEN,
                                      m_input, std::ref(mesh) );
  tk::record< tk::ExodusIIMeshReader >( readers, tk::MeshReaderType::EXODUSII,
                                        m_input, std::ref(mesh) );

  // Register mesh writers
  tk::record< tk::GmshMeshWriter >( writers, tk::MeshWriterType::GMSH,
                                    m_output, std::ref(mesh) );
  tk::record< tk::NetgenMeshWriter >( writers, tk::MeshWriterType::NETGEN,
                                      m_output, std::ref(mesh) );
  tk::record< tk::ExodusIIMeshWriter >( writers, tk::MeshWriterType::EXODUSII,
                                        m_output, std::ref(mesh) );

  // Read in mesh
  tk::instantiate( readers, tk::detectInput( m_input ) )->read();

  // Write out mesh
  tk::instantiate( writers, tk::pickOutput( m_output ) )->write();
}
