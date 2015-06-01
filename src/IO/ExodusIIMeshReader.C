//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.C
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:24:34 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     ExodusII mesh reader
  \details   ExodusII mesh reader class definition. Currently, this is a bare
     minimum functionality to interface with the ExodusII reader. It only reads
     3D meshes and only triangle and tetrahedron elements.
*/
//******************************************************************************

#include <iostream>
#include <array>

#include <exodusII.h>
#include <ne_nemesisI.h>

#include "ExodusIIMeshReader.h"
#include "Exception.h"
#include "Reorder.h"

using tk::ExodusIIMeshReader;

ExodusIIMeshReader::ExodusIIMeshReader( const std::string& filename,
                                        UnsMesh& mesh,
                                        int cpuwordsize,
                                        int iowordsize ) :
  Reader( filename ), m_mesh( mesh ), m_inFile( 0 )
//******************************************************************************
//  Constructor: open Exodus II file
//! \param[in] filename File to open as ExodusII file
//! \param[inout] mesh Unstructured mesh object to load the data to
//! \param[in] cpuwordsize Set CPU word size, see ExodusII documentation
//! \param[in] iowordsize Set I/O word size, see ExodusII documentation
//! \author J. Bakosi
//******************************************************************************
{
  float version;

  m_inFile = ex_open( filename.c_str(), EX_READ, &cpuwordsize, &iowordsize,
                      &version );

  ErrChk( m_inFile > 0, "Failed to open ExodusII file: " + filename );
}

ExodusIIMeshReader::~ExodusIIMeshReader() noexcept
//******************************************************************************
//  Destructor
//! \author J. Bakosi
//******************************************************************************
{
  if ( ex_close(m_inFile) < 0 )
    printf( ">>> WARNING: Failed to close ExodusII file: %s\n",
            m_filename.c_str() );
}

void
ExodusIIMeshReader::read()
//******************************************************************************
//  Read ExodusII mesh file
//! \author J. Bakosi
//******************************************************************************
{
  readHeader();
  readElements();
  readNodes();
}

void
ExodusIIMeshReader::readGraph()
//******************************************************************************
//  Read only connectivity graph from file
//! \author J. Bakosi
//******************************************************************************
{
  readHeader();
  readElements();
}

void
ExodusIIMeshReader::readHeader()
//******************************************************************************
//  Read ExodusII header
//! \author J. Bakosi
//******************************************************************************
{
  char title[MAX_LINE_LENGTH+1];
  int ndim, nel, nnodeset, nelemset, nnode, neblk;

  ErrChk(
    ex_get_init( m_inFile, title, &ndim, &nnode, &nel, &neblk, &nnodeset,
                 &nelemset ) == 0,
    "Failed to read header from ExodusII file: " + m_filename );

  ErrChk( nnode > 0,
          "Number of nodes read from ExodusII file must be larger than zero" );
  ErrChk( neblk > 0,
          "Number of element blocks read from ExodusII file must be larger "
          "than zero" );
  ErrChk( ndim == 3, "Need a 3D mesh from ExodusII file " + m_filename);

  m_neblk = static_cast< std::size_t >( neblk );

  // set mesh graph size
  m_mesh.size() = m_nnode = static_cast< std::size_t >( nnode );
}

void
ExodusIIMeshReader::readNodes()
//******************************************************************************
//  Read all node coordinates from ExodusII file
//! \author J. Bakosi
//******************************************************************************
{
  m_mesh.x().resize( m_nnode );
  m_mesh.y().resize( m_nnode );
  m_mesh.z().resize( m_nnode );

  ErrChk( ex_get_coord( m_inFile, m_mesh.x().data(), m_mesh.y().data(),
                        m_mesh.z().data() ) == 0,
          "Failed to read coordinates from ExodusII file: " + m_filename );
}

void
ExodusIIMeshReader::readNode( std::size_t id,
                              std::vector< tk::real >& x,
                              std::vector< tk::real >& y,
                              std::vector< tk::real >& z )
//******************************************************************************
//  Read coordinates of a single mesh node from ExodusII file
//! \param[in] id Node id whose coordinates to read
//! \param[inout] x Vector of x coordinates to push to
//! \param[inout] y Vector of y coordinates to push to
//! \param[inout] z Vector of z coordinates to push to
//! \author J. Bakosi
//******************************************************************************
{
  tk::real px, py, pz;

  ErrChk(
    ne_get_n_coord(
      m_inFile, static_cast<int64_t>(id)+1, 1, &px, &py, &pz) == 0,
      "Failed to read coordinates of node " + std::to_string( id ) +
      " from ExodusII file: " + m_filename );

  x.push_back( px );
  y.push_back( py );
  z.push_back( pz );
}

void
ExodusIIMeshReader::readElements()
//******************************************************************************
//  Read element blocks and connectivity from ExodusII file
//! \author J. Bakosi
//******************************************************************************
{
  std::vector< int > id( m_neblk );

  // Read element block ids
  ErrChk( ex_get_elem_blk_ids( m_inFile, id.data()) == 0,
          "Failed to read element block ids from ExodusII file: " +
          m_filename );

  for (std::size_t i=0; i<m_neblk; ++i) {
    char eltype[MAX_STR_LENGTH+1];
    int nel, nnpe, nattr;

    // Read element block information
    ErrChk(
      ex_get_elem_block( m_inFile, id[i], eltype, &nel, &nnpe, &nattr ) == 0,
      "Failed to read element block information from ExodusII file: " +
      m_filename );

    // Read element connectivity
    auto connectsize = static_cast< std::size_t >( nel*nnpe );
    if (nnpe == 4) {    // tetrahedra

      m_mesh.tettag().resize( connectsize, { 1 } );
      std::vector< int > inpoel( connectsize );
      ErrChk(
        ex_get_elem_conn( m_inFile, id[i], inpoel.data() ) == 0,
        "Failed to read " + std::string(eltype) + " element connectivity from "
        "ExodusII file: " + m_filename );
      for (auto n : inpoel)
        m_mesh.tetinpoel().push_back( static_cast< std::size_t >( n ) );

    } else if (nnpe == 3) {    // triangles

      m_mesh.tritag().resize( connectsize, { 1 } );
      std::vector< int > inpoel( connectsize );
      ErrChk(
        ex_get_elem_conn( m_inFile, id[i], inpoel.data() ) == 0,
        "Failed to read " + std::string(eltype) + " element connectivity from "
        "ExodusII file: " + m_filename );
      for (auto n : inpoel)
        m_mesh.triinpoel().push_back( static_cast< std::size_t >( n ) );

    }
  }

  // Shift node IDs to start from zero
  shiftToZero( m_mesh.triinpoel() );
  shiftToZero( m_mesh.tetinpoel() );
}
