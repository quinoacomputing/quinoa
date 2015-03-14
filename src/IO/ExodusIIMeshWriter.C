//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshWriter.C
  \author    J. Bakosi
  \date      Fri 13 Mar 2015 09:04:31 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     ExodusII mesh writer
  \details   ExodusII mesh writer class definition. Currently, this is a bare
     minimum functionality to interface with the ExodusII writer. It only writes
     3D meshes and only triangle and tetrahedron elements.
*/
//******************************************************************************

#include <iostream>

#include <exodusII.h>
#include <ne_nemesisI.h>

#include <Config.h>
#include <ExodusIIMeshWriter.h>
#include <Exception.h>

using tk::ExodusIIMeshWriter;

ExodusIIMeshWriter::ExodusIIMeshWriter( const std::string& filename,
                                        const UnsMesh& mesh,
                                        int cpuwordsize,
                                        int iowordsize ) :
  Writer( filename ), m_filename( filename ), m_mesh( mesh ), m_outFile( 0 )
//******************************************************************************
//  Constructor: create Exodus II file
//! \param[in] filename File to open as ExodusII file
//! \param[in] mesh Unstructured mesh object to write data from
//! \param[in] cpuwordsize Set CPU word size, see ExodusII documentation
//! \param[in] iowordsize Set I/O word size, see ExodusII documentation
//! \author J. Bakosi
//******************************************************************************
{
  m_outFile = ex_create( filename.c_str(),
                         EX_CLOBBER | EX_LARGE_MODEL,
                         &cpuwordsize,
                         &iowordsize );

  ErrChk( m_outFile > 0, "Failed to create file: " + filename );
}

ExodusIIMeshWriter::~ExodusIIMeshWriter() noexcept
//******************************************************************************
//  Destructor
//! \author J. Bakosi
//******************************************************************************
{
  if ( ex_close(m_outFile) < 0 )
    printf( ">>> WARNING: Failed to close file: %s\n", m_filename.c_str() );
}

void
ExodusIIMeshWriter::write()
//******************************************************************************
//  Write ExodusII mesh file
//! \author J. Bakosi
//******************************************************************************
{
  writeHeader();
  writeNodes();
  writeElements();
}

void
ExodusIIMeshWriter::writeHeader()
//******************************************************************************
//  Write ExodusII header
//! \author J. Bakosi
//******************************************************************************
{
  ErrChk(
    ex_put_init( m_outFile,
                 "Written by Quinoa",
                 3,     // number of dimensions
                 static_cast< int64_t >( m_mesh.nnode() ),
                 m_mesh.triinpoel().size()/3 + m_mesh.tetinpoel().size()/4,
                 static_cast< int64_t >( m_mesh.neblk() ),
                 0,     // number of node sets
                 0 ) == 0,
    "Failed to write header to file: " + m_filename );
}

void
ExodusIIMeshWriter::writeNodes()
//******************************************************************************
//  Write node coordinates to ExodusII file
//! \author J. Bakosi
//******************************************************************************
{
  ErrChk( ex_put_coord( m_outFile, m_mesh.x().data(), m_mesh.y().data(),
                        m_mesh.z().data() ) == 0,
          "Failed to write coordinates to file: " + m_filename );
}

void
ExodusIIMeshWriter::writeElements()
//******************************************************************************
//  Write element connectivity to ExodusII file
//! \author J. Bakosi
//******************************************************************************
{
  writeElemBlock( 1, 3, "TRIANGLES", m_mesh.triinpoel() );
  writeElemBlock( 2, 4, "TETRAHEDRA", m_mesh.tetinpoel() );
}

void
ExodusIIMeshWriter::writeElemBlock( int elclass,
                                    int nnpe,
                                    const std::string& eltype,
                                    const std::vector< int >& inpoel )
//******************************************************************************
//  Write element block to ExodusII file
//! \param[in] elclass Element class id
//! \param[in] nnpe Number of nodes per element for block
//! \param[in] eltype String describing element type
//! \param[in] inpoel Element connectivity
//! \author J. Bakosi
//******************************************************************************
{
  if (inpoel.empty()) return;

  // Write element block information
  ErrChk(
    ex_put_elem_block( m_outFile,
                       elclass,
                       eltype.c_str(),
                       static_cast< int64_t >( inpoel.size() ) / nnpe,
                       nnpe,
                       0 ) == 0,
    "Failed to write " + eltype + " element block to file: " + m_filename );

  // Write element connectivity with 1-based element ids
  for (int e=0; e<static_cast<int>(inpoel.size()); e+=nnpe) {
    ErrChk(
      ne_put_n_elem_conn( m_outFile,
                          elclass,
                          e/nnpe+1,
                          1,
                          &inpoel[ static_cast< std::size_t >( e ) ] ) == 0,
      "Failed to write " + eltype + " element connectivity to file: " +
        m_filename );
  }
}
