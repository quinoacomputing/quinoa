// *****************************************************************************
/*!
  \file      src/IO/ExodusIIMeshWriter.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     ExodusII mesh-based data writer
  \details   ExodusII mesh-based data writer class definition.
*/
// *****************************************************************************

#include <algorithm>
#include <functional>
#include <iterator>
#include <string>
#include <utility>
#include <cstdint>
#include <cstdio>

#include "NoWarning/exodusII.h"

#include "ExodusIIMeshWriter.h"
#include "Exception.h"
#include "UnsMesh.h"

using tk::ExodusIIMeshWriter;

ExodusIIMeshWriter::ExodusIIMeshWriter( const std::string& filename,
                                        ExoWriter mode,
                                        int cpuwordsize,
                                        int iowordsize ) :
  m_filename( filename ), m_outFile( 0 )
// *****************************************************************************
//  Constructor: create/open Exodus II file
//! \param[in] filename File to open as ExodusII file
//! \param[in] mode ExodusII writer constructor mode: ExoWriter::CREATE for
//!   creating a new file, ExoWriter::OPEN for opening an existing file for
//!   appending
//! \param[in] cpuwordsize Set CPU word size, see ExodusII documentation
//! \param[in] iowordsize Set I/O word size, see ExodusII documentation
//! \author J. Bakosi
// *****************************************************************************
{
  // Increase verbosity from ExodusII library in debug mode
  #ifdef NDEBUG
  ex_opts( EX_DEBUG | EX_VERBOSE );
  #endif

  if (mode == ExoWriter::CREATE) {

    m_outFile = ex_create( filename.c_str(),
                           EX_CLOBBER | EX_LARGE_MODEL,
                           &cpuwordsize,
                           &iowordsize );

  } else if (mode == ExoWriter::OPEN) {

    float version;
    m_outFile = ex_open( filename.c_str(),
                         EX_WRITE, 
                         &cpuwordsize,
                         &iowordsize,
                         &version );

  } else Throw( "Unknown ExodusII writer constructor mode" );

  ErrChk( m_outFile > 0, "Failed to create/open ExodusII file: " + filename );
}

ExodusIIMeshWriter::~ExodusIIMeshWriter() noexcept
// *****************************************************************************
//  Destructor
//! \author J. Bakosi
// *****************************************************************************
{
  if ( ex_close(m_outFile) < 0 )
    printf( ">>> WARNING: Failed to close ExodusII file: %s\n",
            m_filename.c_str() );
}

void
ExodusIIMeshWriter::writeMesh( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write ExodusII mesh file
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  writeHeader( mesh );
  writeNodes( mesh );
  writeElements( mesh );
}

void
ExodusIIMeshWriter::writeHeader( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write ExodusII header
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  ErrChk(
    ex_put_init( m_outFile,
                 "Written by Quinoa",
                 3,     // number of dimensions
                 static_cast< int64_t >( mesh.nnode() ),
                 static_cast< int64_t >( mesh.triinpoel().size()/3 +
                                         mesh.tetinpoel().size()/4 ),
                 static_cast< int64_t >( mesh.neblk() ),
                 0,     // number of node sets
                 0 ) == 0,
    "Failed to write header to file: " + m_filename );
}

void
ExodusIIMeshWriter::writeNodes( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write node coordinates to ExodusII file
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  ErrChk( ex_put_coord( m_outFile, mesh.x().data(), mesh.y().data(),
                        mesh.z().data() ) == 0,
          "Failed to write coordinates to ExodusII file: " + m_filename );
}

void
ExodusIIMeshWriter::writeElements( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write element connectivity to ExodusII file
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  int elclass = 0;

  writeElemBlock( elclass, 3, "TRIANGLES", mesh.triinpoel() );
  writeElemBlock( elclass, 4, "TETRAHEDRA", mesh.tetinpoel() );
}

void
ExodusIIMeshWriter::writeElemBlock( int& elclass,
                                    int64_t nnpe,
                                    const std::string& eltype,
                                    const std::vector< std::size_t >& inpoel )
const
// *****************************************************************************
//  Write element block to ExodusII file
//! \param[inout] elclass Count element class ids in file
//! \param[in] nnpe Number of nodes per element for block
//! \param[in] eltype String describing element type
//! \param[in] inpoel Element connectivity.
//! \author J. Bakosi
// *****************************************************************************
{
  if (inpoel.empty()) return;

  // increase number of element classes in file
  ++elclass;

  // Make sure element connectivity starts with zero
  Assert( *std::minmax_element( begin(inpoel), end(inpoel) ).first == 0,
          "node ids should start from zero" );

  // Write element block information
  ErrChk(
    ex_put_block( m_outFile,            // exo file handle
                  EX_ELEM_BLOCK,        // block type: elem block
                  elclass,              // element class id
                  eltype.c_str(),       // element block description
                  static_cast< int64_t >( inpoel.size() ) / nnpe, // nof cells
                  nnpe, // number of nodes per element
                  6,    // number of edges per element
                  4,    // number of faces per element
                  0     // number of attributes per element
                ) == 0,
    "Failed to write " + eltype + " element block to ExodusII file: " +
    m_filename );

  // Write element connectivity with 1-based node ids
  std::vector< int > inp;
  for (auto p : inpoel) inp.push_back( static_cast< int >( p+1 ) );
  ErrChk( ex_put_conn( m_outFile, EX_ELEM_BLOCK, elclass, inp.data(),
                       nullptr, nullptr ) == 0,
          "Failed to write " + eltype + " element connectivity to ExodusII "
          "file: " + m_filename );
}

void
ExodusIIMeshWriter::writeTimeStamp( uint64_t it, tk::real time ) const
// *****************************************************************************
//  Write time stamp to ExodusII file
//! \param[in] it Iteration number
//! \param[in] time Time
//! \author J. Bakosi
// *****************************************************************************
{
  ErrChk( ex_put_time( m_outFile, static_cast<int>(it), &time ) == 0,
          "Failed to time stamp to ExodusII file: " + m_filename );
}

void
ExodusIIMeshWriter::writeNodeVarNames( const std::vector< std::string >& nv )
const
// *****************************************************************************
//  Write the names of nodal output variables to ExodusII file
//! \param[in] nv Nodal variable names
//! \author J. Bakosi
// *****************************************************************************
{
  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  ErrChk(
    ex_put_variable_param( m_outFile, EX_NODE_BLOCK,
                           static_cast< int >( nv.size() ) ) == 0,
    "Failed to write nodal output variable parameters to ExodusII file: " +
    m_filename );

  char* names[ nv.size() ];
  std::size_t i = 0;
  for (const auto& n : nv) names[ i++ ] = const_cast< char* >( n.c_str() );

  ErrChk( ex_put_variable_names( m_outFile,
                                 EX_NODE_BLOCK,
                                 static_cast< int >( nv.size() ),
                                 names ) == 0,
          "Failed to write nodal output variable names to ExodusII file: " +
          m_filename );

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

void
ExodusIIMeshWriter::writeElemVarNames( const std::vector< std::string >& ev )
const
// *****************************************************************************
//  Write the names of element output variables to ExodusII file
//! \param[in] ev Elem variable names
//! \author J. Bakosi
// *****************************************************************************
{
  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  ErrChk(
    ex_put_variable_param( m_outFile, EX_ELEM_BLOCK,
                           static_cast< int >( ev.size() ) ) == 0,
    "Failed to write element output variable parameters to ExodusII file: " +
    m_filename );

  char* names[ ev.size() ];
  std::size_t i = 0;
  for (const auto& n : ev) names[ i++ ] = const_cast< char* >( n.c_str() );

  ErrChk( ex_put_variable_names( m_outFile,
                                 EX_ELEM_BLOCK,
                                 static_cast< int >( ev.size() ),
                                 names ) == 0,
          "Failed to write element output variable names to ExodusII file: " +
          m_filename );

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

void
ExodusIIMeshWriter::writeNodeScalar( uint64_t it,
                                     int varid,
                                     const std::vector< tk::real >& var ) const
// *****************************************************************************
//  Write node scalar field to ExodusII file
//! \param[in] it Iteration number
//! \param[in] varid Variable id
//! \param[in] var Vector of variable to output
//! \author J. Bakosi
// *****************************************************************************
{
  ErrChk( ex_put_var( m_outFile,
                      static_cast< int >( it ),
                      EX_NODE_BLOCK,
                      varid,
                      1,
                      static_cast< int64_t >( var.size() ),
                      var.data() ) == 0,
          "Failed to write node scalar to ExodusII file: " + m_filename );
}

void
ExodusIIMeshWriter::writeElemScalar( uint64_t it,
                                     int varid,
                                     const std::vector< tk::real >& var ) const
// *****************************************************************************
//  Write elem scalar field to ExodusII file
//! \param[in] it Iteration number
//! \param[in] varid Variable id
//! \param[in] var Vector of variable to output
//! \author J. Bakosi
// *****************************************************************************
{
  ErrChk( ex_put_var( m_outFile,
                      static_cast< int >( it ),
                      EX_ELEM_BLOCK,
                      varid,
                      1,
                      static_cast< int64_t >( var.size() ),
                      var.data() ) == 0,
          "Failed to write elem scalar to ExodusII file: " + m_filename );
}
