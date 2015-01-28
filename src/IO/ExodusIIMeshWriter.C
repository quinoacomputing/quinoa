//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshWriter.C
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 08:32:01 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
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

using quinoa::ExodusIIMeshWriter;

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
                 (std::string("Written by Quinoa::") +
                   MESHCONV_EXECUTABLE).c_str(),
                 3,     // number of dimensions
                 m_mesh.nnode(),
                 m_mesh.triinpoel().size() + m_mesh.tetinpoel().size(),
                 m_mesh.neblk(),
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
ExodusIIMeshWriter::writeElemBlock(
  int elclass,
  int nnpe,
  const std::string& eltype,
  const std::vector< std::vector< int > >& inpoel )
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
  ErrChk( ex_put_elem_block( m_outFile, elclass, eltype.c_str(), inpoel.size(),
                             nnpe, 0 ) == 0,
          "Failed to write " + eltype + " element block to file: " +
          m_filename );

  // Write element connectivity
  for (std::size_t e=0; e<inpoel.size(); ++e) {
    ErrChk( ne_put_n_elem_conn( m_outFile, elclass, e+1, 1,
                                inpoel[e].data() ) == 0,  // 1-based element ids
            "Failed to write " + eltype + " element connectivity to file: " +
            m_filename );
  }
}
