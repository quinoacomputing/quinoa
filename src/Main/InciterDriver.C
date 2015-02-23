//******************************************************************************
/*!
  \file      src/Main/InciterDriver.C
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 03:40:10 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter driver
  \details   Inciter driver.
*/
//******************************************************************************

#include <InciterDriver.h>
#include <Inciter/InputDeck/Parser.h>
#include <MeshDetect.h>
#include <GmshMeshReader.h>
#include <NetgenMeshReader.h>
#include <ExodusIIMeshReader.h>
#include <NetgenMeshWriter.h>
#include <GmshMeshWriter.h>
#include <ExodusIIMeshWriter.h>
#include <inciter.decl.h>

extern CProxy_Main mainProxy;

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::InciterDriver;

InciterDriver::InciterDriver( const InciterPrint& print,
                              const ctr::CmdLine& cmdline ) :
  m_print( print ),
  m_input( cmdline.get< tag::io, tag::input >() )
//******************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
//! \author J. Bakosi
//******************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here (if any)

  // Parse input deck into g_inputdeck
  InputDeckParser inputdeckParser( m_print, cmdline, g_inputdeck );

  m_print.endpart();
  m_print.part( "Factory" );
}

void
InciterDriver::execute() const
//******************************************************************************
//  Execute: Read mesh from file, partition to multiple chunks
//! \author J. Bakosi
//******************************************************************************
{
  //! Mesh readers factory
  std::map< tk::MeshReaderType, std::function<tk::Reader*()> > readers;

  //! Create unstructured mesh to store mesh
  tk::UnsMesh mesh;

  // Register mesh readers
  tk::record< tk::GmshMeshReader >( readers, tk::MeshReaderType::GMSH,
                                    m_input, std::ref(mesh) );
  tk::record< tk::NetgenMeshReader >( readers, tk::MeshReaderType::NETGEN,
                                      m_input, std::ref(mesh) );
  tk::record< tk::ExodusIIMeshReader >( readers, tk::MeshReaderType::EXODUSII,
                                        m_input, std::ref(mesh) );

  // Read in mesh
  tk::instantiate( readers, tk::detectInput( m_input ) )->read();

  // Echo mesh statistics
  meshStats( mesh );

  mainProxy.finalize();
}

void
InciterDriver::meshStats( const tk::UnsMesh& mesh ) const
//******************************************************************************
//  Execute: Echo mesh statistics
//! \param[in] mesh Unstructured mesh object reference to echo stats of
//! \author J. Bakosi
//******************************************************************************
{
  // Echo mesh stats in verbose mode
  m_print.section( "Input mesh statistics" );
  m_print.item( "Number of element blocks", mesh.neblk() );
  m_print.item( "Number of elements", mesh.nelem() );
  m_print.item( "Number of nodes", mesh.nnode() );

  if (!mesh.lininpoel().empty())
    m_print.item( "Number of lines", mesh.lininpoel().size() );
  if (!mesh.triinpoel().empty())
    m_print.item( "Number of triangles", mesh.triinpoel().size() );
  if (!mesh.tetinpoel().empty())
    m_print.item( "Number of tetrahedra", mesh.tetinpoel().size() );
}
