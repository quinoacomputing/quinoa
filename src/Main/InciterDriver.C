//******************************************************************************
/*!
  \file      src/Main/InciterDriver.C
  \author    J. Bakosi
  \date      Wed 25 Feb 2015 07:47:22 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter driver
  \details   Inciter driver.
*/
//******************************************************************************

#include <zoltan.h>

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
       mesh.tetinpoel().size(),
       chunksize,
       remainder );

  // Print out info on what will be done and how
  info( mesh, chunksize, remainder, nchare );

  // Partition mesh using Zoltan
  partition( mesh );

  mainProxy.finalize();
}

void
InciterDriver::partition( const tk::UnsMesh& mesh ) const
//******************************************************************************
//  Partition mesh using Zoltan
//! \param[in] mesh Unstructured mesh object reference to echo stats of
//! \author J. Bakosi
//******************************************************************************
{
//   // Initialize the Zoltan library
//   float ver = 0.0;
//   ErrChk( Zoltan_Initialize( 0, nullptr, &ver ) == ZOLTAN_OK,
//           "Zoltan could not be initialized" );
//
//   // Create Zoltan data structure
//   struct Zoltan_Struct *z;
//   z = Zoltan_Create( MPI_COMM_WORLD );
//   Assert( z != nullptr, "Zoltan_Create failed" );
//
//   // Destroy Zoltan data structure
//   Zoltan_Destroy( &z );
}

void
InciterDriver::info( const tk::UnsMesh& mesh,
                     uint64_t chunksize,
                     uint64_t remainder,
                     uint64_t nchare ) const
//******************************************************************************
//  Print information at startup
//! \param[in] mesh Unstructured mesh object reference to echo stats of
//! \param[in] chunksize Chunk size, see Base/LoadDistribution.h
//! \param[in] remainder Remainder, see Base/LoadDistribution.h
//! \param[in] nchare Number of work units (Charm++ chares)
//! \author J. Bakosi
//******************************************************************************
{
  // Print out mesh stats
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

  // Print out info on load distribution
  m_print.section( "Load distribution" );
  m_print.item( "Virtualization [0.0...1.0]",
                g_inputdeck.get< tag::cmd, tag::virtualization >() );
  m_print.item( "Load (number of mesh cells)", mesh.tetinpoel().size() );
  m_print.item( "Number of processing elements", CkNumPes() );
  m_print.item( "Number of work units",
                std::to_string( nchare ) + " (" +
                std::to_string( nchare-1 ) + "*" +
                std::to_string( chunksize ) + "+" +
                std::to_string( chunksize+remainder ) + ")" );
}
