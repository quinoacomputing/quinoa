//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshWriter.C
  \author    J. Bakosi
  \date      Sun 12 Apr 2015 08:39:25 AM MDT
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

  // Note that we cast away the constness of the element connectivities
  // explicitly below. This could be avoided if the mesh object reference,
  // m_mesh, was simply held as a non-const reference. However, that would allow
  // all other member functions in this class to modify the mesh object via
  // their non-const member functions. Instead, that is not allowed (after all a
  // writer should not modify the internal state of the object being written),
  // and the constness of the mesh connectivity references, obtained by the
  // element connectivity mesh object member functions, triinpoel() and
  // tetinpoel(), are specifically and locally casted away so that the element
  // block writer can temporarily shift the zero-based point ids to one-based
  // ones, write the connectivities out to file, and then shift the point ids
  // back. If everything goes well, this indeed does not modify the
  // connectivities (as is visible to the outside. The rationale is that this
  // avoids copies just for outputing to file and is thus significantly faster.
  // Also the full connectivity is output at once, using the single ExodusII I/O
  // call, ex_put_elem_conn, instead of having to use the NemesisI I/O call,
  // ne_put_n_elem_conn, outputing a single element's connectivity at a time.
  // See also member function writeElemBlock().

  writeElemBlock( elclass, 3, "TRIANGLES",
                  const_cast< std::vector< int >& >( m_mesh.triinpoel() ) );
  writeElemBlock( elclass, 4, "TETRAHEDRA",
                  const_cast< std::vector< int >& >( m_mesh.tetinpoel() ) );
}

void
ExodusIIMeshWriter::writeElemBlock( int& elclass,
                                    int nnpe,
                                    const std::string& eltype,
                                    std::vector< int >& inpoel )
//******************************************************************************
//  Write element block to ExodusII file
//! \param[inout] elclass Count element class ids in file
//! \param[in] nnpe Number of nodes per element for block
//! \param[in] eltype String describing element type
//! \param[in] inpoel Element connectivity. Note that inpoel is modified
//!   in-place to avoid copying memory: the node coordinates are shifted to be
//!   one-based, written out, then shifted back to zero-based.
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
  for (auto& p : inpoel) ++p;   // shift node ids to one-based
  ErrChk( ex_put_elem_conn( m_outFile, elclass, inpoel.data() ) == 0,
          "Failed to write " + eltype + " element connectivity to ExodusII "
          "file: " + m_filename );
  for (auto& p : inpoel) --p;   // shift back node ids to zero-based
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
ExodusIIMeshWriter::writeVarNames( const std::vector< std::string >& nvar )
//******************************************************************************
//  Write the number and names of output variables to ExodusII file
//! \param[in] nvar Variable names
//! \author J. Bakosi
//******************************************************************************
{
  ErrChk(
    ex_put_var_param( m_outFile, "n", static_cast<int>(nvar.size()) ) == 0,
    "Failed to write the number of output variables to ExodusII file: " +
    m_filename );

  std::vector< const char* > names;
  std::transform( std::begin(nvar), std::end(nvar), std::back_inserter(names),
                  std::mem_fn(&std::string::c_str) );

  ErrChk( ex_put_var_names( m_outFile,
                            "n",
                            static_cast<int>(nvar.size()),
                            const_cast<char**>(names.data()) ) == 0,
          "Failed to write the number of output variables to ExodusII file: " +
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
