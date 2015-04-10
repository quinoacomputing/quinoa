//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.C
  \author    J. Bakosi
  \date      Thu 09 Apr 2015 12:38:44 PM MDT
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

#include <ExodusIIMeshReader.h>
#include <Exception.h>
#include <DerivedData.h>

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
    printf( ">>> WARNING: Failed to close file: %s\n", m_filename.c_str() );
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
  int ndim, nel, nnodeset, nelemset;

  ErrChk(
    ex_get_init( m_inFile, title, &ndim, &m_nnode, &nel, &m_neblk, &nnodeset,
                 &nelemset ) == 0,
    "Failed to read header from file: " + m_filename );

  ErrChk( m_nnode > 0,
          "Number of nodes read from file must be larger than zero" );
  ErrChk( m_neblk > 0,
          "Number of element blocks read from file must be larger than zero" );
  ErrChk( ndim == 3, "Need a 3D mesh from file " + m_filename);

  // set mesh graph size
  m_mesh.size() = static_cast< std::size_t >( m_nnode );
}

void
ExodusIIMeshReader::readNodes()
//******************************************************************************
//  Read node coordinates from ExodusII file
//! \author J. Bakosi
//******************************************************************************
{
  m_mesh.x().resize( static_cast< std::size_t >( m_nnode ) );
  m_mesh.y().resize( static_cast< std::size_t >( m_nnode ) );
  m_mesh.z().resize( static_cast< std::size_t >( m_nnode ) );

  ErrChk( ex_get_coord( m_inFile, m_mesh.x().data(), m_mesh.y().data(),
                        m_mesh.z().data() ) == 0,
          "Failed to read coordinates from file: " + m_filename );
}

void
ExodusIIMeshReader::readNode( std::size_t id,
                              std::vector< tk::real >& x,
                              std::vector< tk::real >& y,
                              std::vector< tk::real >& z )
//******************************************************************************
//  Read the coordinates of a single mesh node from ExodusII file
//! \param[in] id Node id whose coordinates to read
//! \param[inout] x Vector of x coordinates to push to
//! \param[inout] y Vector of y coordinates to push to
//! \param[inout] z Vector of z coordinates to push to
//! \author J. Bakosi
//******************************************************************************
{
  tk::real px, py, pz;

  ErrChk(
    ne_get_n_coord(m_inFile, static_cast<int64_t>(id)+1, 1, &px, &py, &pz) == 0,
    "Failed to read coordinates of node " + std::to_string( id ) +
      " from file: " + m_filename );

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
  std::vector< int > id( static_cast< std::size_t >( m_neblk ) );
  char eltype[MAX_STR_LENGTH+1];
  int nel, nnpe, nattr;

  // Read element block ids
  ErrChk( ex_get_elem_blk_ids( m_inFile, id.data()) == 0,
          "Failed to read element block ids from file: " + m_filename );

  for (std::size_t i=0; i<static_cast<std::size_t>(m_neblk); ++i) {
    // Read element block information
    ErrChk(
      ex_get_elem_block( m_inFile, id[i], eltype, &nel, &nnpe, &nattr ) == 0,
      "Failed to read element block information from file: " + m_filename );

    // Read element connectivity
    if (nnpe == 4) {    // tetrahedra
      for (int e=0; e<nel; ++e) {
        std::array< int, 4 > n;
        ErrChk( ne_get_n_elem_conn( m_inFile, id[i], e+1, 1, n.data() ) == 0,
                "Failed to read " + std::string(eltype) +
                " element connectivity from file: " + m_filename );
        m_mesh.tettag().push_back( { 1 } );
        m_mesh.tetinpoel().push_back( n[0] );
        m_mesh.tetinpoel().push_back( n[1] );
        m_mesh.tetinpoel().push_back( n[2] );
        m_mesh.tetinpoel().push_back( n[3] );
      }
    } else if (nnpe == 3) {    // triangles
      for (int e=0; e<nel; ++e) {
        std::array< int, 3 > n;
        ErrChk( ne_get_n_elem_conn( m_inFile, id[i], e+1, 1, n.data() ) == 0,
                "Failed to read " + std::string(eltype) +
                " element connectivity from file: " + m_filename );
        m_mesh.tritag().push_back( { 1 } );
        m_mesh.triinpoel().push_back( n[0] );
        m_mesh.triinpoel().push_back( n[1] );
        m_mesh.triinpoel().push_back( n[2] );
      }
    }
  }

  // Shift node IDs to start from zero
  shiftToZero( m_mesh.triinpoel() );
  shiftToZero( m_mesh.tetinpoel() );
}
