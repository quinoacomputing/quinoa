//******************************************************************************
/*!
  \file      src/Main/InciterDriver.C
  \author    J. Bakosi
  \date      Fri 27 Feb 2015 03:27:40 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter driver
  \details   Inciter driver.
*/
//******************************************************************************

#include <InciterDriver.h>
#include <Inciter/InputDeck/Parser.h>
#include <LoadDistributor.h>
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
extern tk::UnsMesh g_mesh;

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

  // Register mesh readers
  tk::record< tk::GmshMeshReader >( readers, tk::MeshReaderType::GMSH,
                                    m_input, std::ref(g_mesh) );
  tk::record< tk::NetgenMeshReader >( readers, tk::MeshReaderType::NETGEN,
                                      m_input, std::ref(g_mesh) );
  tk::record< tk::ExodusIIMeshReader >( readers, tk::MeshReaderType::EXODUSII,
                                        m_input, std::ref(g_mesh) );

  // Start timer measuring mesh read time
  tk::Timer mesh_read;

  // Read in mesh
  tk::instantiate( readers, tk::detectInput( m_input ) )->read();

  // Report mesh read time to main proxy
  mainProxy.timestamp( "Mesh read", mesh_read.dsec() );

  // Compute load distribution given total work (= number of mesh cells) and
  // user-specified virtualization
  uint64_t chunksize, remainder;
  const auto nchare =
    tk::linearLoadDistributor(
       g_inputdeck.get< tag::cmd, tag::virtualization >(),
       g_mesh.tetinpoel().size(),
       chunksize,
       remainder );

  // Print out info on what will be done and how
  info( chunksize, remainder, nchare );

  mainProxy.finalize();
}

void
InciterDriver::info( uint64_t chunksize, uint64_t remainder, uint64_t nchare )
const
//******************************************************************************
//  Print information at startup
//! \param[in] chunksize Chunk size, see Base/LoadDistribution.h
//! \param[in] remainder Remainder, see Base/LoadDistribution.h
//! \param[in] nchare Number of work units (Charm++ chares)
//! \author J. Bakosi
//******************************************************************************
{
  // Print out mesh stats
  m_print.section( "Input mesh statistics" );
  m_print.item( "Number of element blocks", g_mesh.neblk() );
  m_print.item( "Number of elements", g_mesh.nelem() );
  m_print.item( "Number of nodes", g_mesh.nnode() );

  if (!g_mesh.lininpoel().empty())
    m_print.item( "Number of lines", g_mesh.lininpoel().size() );
  if (!g_mesh.triinpoel().empty())
    m_print.item( "Number of triangles", g_mesh.triinpoel().size() );
  if (!g_mesh.tetinpoel().empty())
    m_print.item( "Number of tetrahedra", g_mesh.tetinpoel().size() );

  // Print out info on load distribution
  m_print.section( "Load distribution" );
  m_print.item( "Virtualization [0.0...1.0]",
                g_inputdeck.get< tag::cmd, tag::virtualization >() );
  m_print.item( "Load (number of mesh cells)", g_mesh.tetinpoel().size() );
  m_print.item( "Number of processing elements", CkNumPes() );
  m_print.item( "Number of work units",
                std::to_string( nchare ) + " (" +
                std::to_string( nchare-1 ) + "*" +
                std::to_string( chunksize ) + "+" +
                std::to_string( chunksize+remainder ) + ")" );
}
