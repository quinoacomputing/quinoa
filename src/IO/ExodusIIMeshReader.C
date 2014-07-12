//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.C
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:30:52 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     ExodusII mesh reader
  \details   ExodusII mesh reader
*/
//******************************************************************************

#include <iostream>

#include <exodusII.h>
#include <ne_nemesisI.h>

#include <Config.h>
#include <ExodusIIMeshReader.h>
#include <Exception.h>

using quinoa::ExodusIIMeshReader;

ExodusIIMeshReader::ExodusIIMeshReader( const std::string& filename,
                                        UnsMesh& mesh,
                                        int cpuwordsize,
                                        int iowordsize ) :
  Reader(filename), m_filename(filename), m_mesh(mesh), m_inFile(0)
//******************************************************************************
//  Constructor: create Exodus II file
//! \author J. Bakosi
//******************************************************************************
{
  float version;

  m_inFile = ex_open( filename.c_str(), EX_READ, &cpuwordsize, &iowordsize,
                      &version );

  ErrChk( m_inFile > 0, "Failed to open file: " + filename );
}

ExodusIIMeshReader::~ExodusIIMeshReader() noexcept
//******************************************************************************
//  Destructor
//! \author J. Bakosi
//******************************************************************************
{
  if ( ex_close(m_inFile) < 0 )
    printf( ">>> WARNING: Failed to close file: %s\n", m_filename );
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

  ErrChk( ndim == 3, "Need a 3D mesh from file " + m_filename);
}

void
ExodusIIMeshReader::readNodes()
//******************************************************************************
//  Read node coordinates from ExodusII file
//! \author J. Bakosi
//******************************************************************************
{
  m_mesh.x().reserve( m_nnode );
  m_mesh.y().reserve( m_nnode );
  m_mesh.z().reserve( m_nnode );

  ErrChk( ex_get_coord( m_inFile, m_mesh.x().data(), m_mesh.y().data(),
                        m_mesh.z().data() ) == 0,
          "Failed to read coordinates from file: " + m_filename );

  for (int i=1; i<=m_nnode; ++i) {
    m_mesh.nodeId().push_back( i );
  }
}

void
ExodusIIMeshReader::readElements()
//******************************************************************************
//  Read element blocks and connectivity from ExodusII file
//! \author J. Bakosi
//******************************************************************************
{
  std::vector< int > id( m_neblk );
  char eltype[MAX_STR_LENGTH+1];
  int nel, nnpe, nattr;

  // Read element block ids
  ErrChk( ex_get_elem_blk_ids( m_inFile, id.data()) == 0,
          "Failed to read element block ids from file: " + m_filename );

  for (int i=0; i<m_neblk; ++i) {
    // Read element block information
    ErrChk(
      ex_get_elem_block( m_inFile, id[i], eltype, &nel, &nnpe, &nattr ) == 0,
      "Failed to read element block information from file: " + m_filename );
    
    // Read element connectivity
    if (nnpe == 4) {    // tetrahedra
      for (int e=0; e<nel; ++e) {
        std::vector< int > n( 4 );
        ErrChk( ne_get_n_elem_conn( m_inFile, id[i], e+1, 1, n.data() ) == 0,
                "Failed to read " + std::string(eltype) +
                " element connectivity from file: " + m_filename );
        m_mesh.tetId().push_back( e );
        m_mesh.tettag().push_back( { 1 } );
        m_mesh.tetinpoel().push_back( n );
      }
    } else if (nnpe == 3) {    // triangles
      for (int e=0; e<nel; ++e) {
        std::vector< int > n( 3 );
        ErrChk( ne_get_n_elem_conn( m_inFile, id[i], e+1, 1, n.data() ) == 0,
                "Failed to read " + std::string(eltype) +
                " element connectivity from file: " + m_filename );
        m_mesh.triId().push_back( e );
        m_mesh.tritag().push_back( { 1 } );
        m_mesh.triinpoel().push_back( n );
      }
    }
  }
}
