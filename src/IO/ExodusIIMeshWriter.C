//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshWriter.C
  \author    J. Bakosi
  \date      Mon 27 Apr 2015 03:40:11 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     ExodusII mesh writer
  \details   ExodusII mesh writer class definition. Currently, this is a bare
     minimum functionality to interface with the ExodusII writer. It only writes
     3D meshes and only triangle and tetrahedron elements.
*/
//******************************************************************************

#include <iostream>

#include <exodusII.h>

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

  ErrChk( m_outFile > 0, "Failed to create ExodusII file: " + filename );
}

ExodusIIMeshWriter::~ExodusIIMeshWriter() noexcept
//******************************************************************************
//  Destructor
//! \author J. Bakosi
//******************************************************************************
{
  if ( ex_close(m_outFile) < 0 )
    printf( ">>> WARNING: Failed to close ExodusII file: %s\n",
            m_filename.c_str() );
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
          "Failed to write coordinates to ExodusII file: " + m_filename );
}

void
ExodusIIMeshWriter::writeElements()
//******************************************************************************
//  Write element connectivity to ExodusII file
//! \author J. Bakosi
//******************************************************************************
{
  int elclass = 0;

  writeElemBlock( elclass, 3, "TRIANGLES", m_mesh.triinpoel() );
  writeElemBlock( elclass, 4, "TETRAHEDRA", m_mesh.tetinpoel() );
}

void
ExodusIIMeshWriter::writeElemBlock( int& elclass,
                                    int nnpe,
                                    const std::string& eltype,
                                    const std::vector< std::size_t >& inpoel )
//******************************************************************************
//  Write element block to ExodusII file
//! \param[inout] elclass Count element class ids in file
//! \param[in] nnpe Number of nodes per element for block
//! \param[in] eltype String describing element type
//! \param[in] inpoel Element connectivity.
//! \author J. Bakosi
//******************************************************************************
{
  if (inpoel.empty()) return;

  // increase number of element classes in file
  ++elclass;

  // Make sure element connectivity starts with zero
  Assert( *std::minmax_element( begin(inpoel), end(inpoel) ).first == 0,
          "node ids should start from zero" );

  // Write element block information
  ErrChk(
    ex_put_elem_block( m_outFile,
                       elclass,
                       eltype.c_str(),
                       static_cast< int64_t >( inpoel.size() ) / nnpe,
                       nnpe,
                       0 ) == 0,
    "Failed to write " + eltype + " element block to ExodusII file: " +
    m_filename );

  // Write element connectivity with 1-based element ids
  std::vector< int > inp;
  for (auto p : inpoel) inp.push_back( static_cast< int >( p+1 ) );
  ErrChk( ex_put_elem_conn( m_outFile, elclass, inp.data() ) == 0,
          "Failed to write " + eltype + " element connectivity to ExodusII "
          "file: " + m_filename );
}

void
ExodusIIMeshWriter::writeTimeStamp( int it, tk::real time )
//******************************************************************************
//  Write time stamp to ExodusII file
//! \param[in] it Iteration number
//! \param[in] time Time
//! \author J. Bakosi
//******************************************************************************
{
  ErrChk( ex_put_time( m_outFile, it, &time ) == 0,
          "Failed to time stamp to ExodusII file: " + m_filename );
}

void
ExodusIIMeshWriter::writeNodeVarNames( const std::vector< std::string >& nv )
//******************************************************************************
//  Write the names of nodal output variables to ExodusII file
//! \param[in] nv Nodal variable names
//! \author J. Bakosi
//******************************************************************************
{
  ErrChk(
    ex_put_var_param( m_outFile, "n", static_cast<int>(nv.size()) ) == 0,
    "Failed to write nodal output variable parameters to ExodusII file: " +
    m_filename );

  std::vector< const char* > names;
  std::transform( std::begin(nv), std::end(nv), std::back_inserter(names),
                  std::mem_fn(&std::string::c_str) );

  ErrChk( ex_put_var_names( m_outFile,
                            "n",
                            static_cast<int>(nv.size()),
                            const_cast<char**>(names.data()) ) == 0,
          "Failed to write nodal output variable names to ExodusII file: " +
          m_filename );
}

void
ExodusIIMeshWriter::writeElemVarNames( const std::vector< std::string >& ev )
//******************************************************************************
//  Write the names of element output variables to ExodusII file
//! \param[in] ev Elem variable names
//! \author J. Bakosi
//******************************************************************************
{
  ErrChk(
    ex_put_var_param( m_outFile, "e", static_cast<int>(ev.size()) ) == 0,
    "Failed to write element output variable parameters to ExodusII file: " +
    m_filename );

  std::vector< const char* > names;
  std::transform( std::begin(ev), std::end(ev), std::back_inserter(names),
                  std::mem_fn(&std::string::c_str) );

  ErrChk( ex_put_var_names( m_outFile,
                            "e",
                            static_cast<int>(ev.size()),
                            const_cast<char**>(names.data()) ) == 0,
          "Failed to write element output variable names to ExodusII file: " +
          m_filename );
}

void
ExodusIIMeshWriter::writeNodeScalar( int it,
                                     int varid,
                                     const std::vector< tk::real >& var )
//******************************************************************************
//  Write node scalar field to ExodusII file
//! \param[in] it Iteration number
//! \param[in] varid Variable id
//! \param[in] var Vector of variable to output
//! \author J. Bakosi
//******************************************************************************
{
  ErrChk( ex_put_nodal_var( m_outFile,
                            it,
                            varid,
                            static_cast< int64_t >( var.size() ),
                            var.data() ) == 0,
          "Failed to write node scalar to ExodusII file: " + m_filename );
}

void
ExodusIIMeshWriter::writeElemScalar( int it,
                                     int varid,
                                     const std::vector< tk::real >& var )
//******************************************************************************
//  Write elem scalar field to ExodusII file
//! \param[in] it Iteration number
//! \param[in] varid Variable id
//! \param[in] var Vector of variable to output
//! \author J. Bakosi
//******************************************************************************
{
  ErrChk( ex_put_elem_var( m_outFile,
                           it,
                           varid,
                           1,
                           static_cast< int64_t >( var.size() ),
                           var.data() ) == 0,
          "Failed to write elem scalar to ExodusII file: " + m_filename );
}
