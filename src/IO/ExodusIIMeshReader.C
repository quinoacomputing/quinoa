//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.C
  \author    J. Bakosi
  \date      Thu 19 Nov 2015 03:18:28 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     ExodusII mesh reader
  \details   ExodusII mesh reader class definition. Currently, this is a bare
     minimum functionality to interface with the ExodusII reader. It only reads
     3D meshes and only triangle and tetrahedron elements.
*/
//******************************************************************************

#include <cstdint>
#include <cstdio>
#include <string>
#include <numeric>
#include <unordered_map>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <exodusII.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <ne_nemesisI.h>

#include "ExodusIIMeshReader.h"
#include "Exception.h"
#include "UnsMesh.h"
#include "Reorder.h"

using tk::ExodusIIMeshReader;

ExodusIIMeshReader::ExodusIIMeshReader( const std::string& filename,
                                        int cpuwordsize,
                                        int iowordsize ) :
  m_filename( filename ),
  m_inFile( 0 ),
  m_eidt( m_nnpe.size(), -1 ),
  m_nel( m_nnpe.size(), -1 )
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
ExodusIIMeshReader::readMesh( UnsMesh& mesh )
//******************************************************************************
//  Read ExodusII mesh file
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
//******************************************************************************
{
  readHeader( mesh );
  readAllElements( mesh );
  readAllNodes( mesh );
}

void
ExodusIIMeshReader::readGraph( UnsMesh& mesh )
//******************************************************************************
//  Read only connectivity graph from file
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
//******************************************************************************
{
  readHeader( mesh );
  readAllElements( mesh );
}

std::size_t
ExodusIIMeshReader::readHeader()
//******************************************************************************
//  Read ExodusII header without setting mesh size
//! \return Number of nodes in mesh
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

  return static_cast< std::size_t >( nnode );
}

void
ExodusIIMeshReader::readHeader( UnsMesh& mesh )
//******************************************************************************
//  Read ExodusII header with setting mesh size
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
//******************************************************************************
{
  // Read ExodusII file header and set mesh graph size
  mesh.size() = m_nnode = static_cast< std::size_t >( readHeader() );
}

void
ExodusIIMeshReader::readAllNodes( UnsMesh& mesh ) const
//******************************************************************************
//  Read all node coordinates from ExodusII file
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
//******************************************************************************
{
  mesh.x().resize( m_nnode );
  mesh.y().resize( m_nnode );
  mesh.z().resize( m_nnode );

  ErrChk( ex_get_coord( m_inFile, mesh.x().data(), mesh.y().data(),
                        mesh.z().data() ) == 0,
          "Failed to read coordinates from ExodusII file: " + m_filename );
}

std::unordered_map< std::size_t, std::array< tk::real, 3 > >
ExodusIIMeshReader::readNodes( const std::array< std::size_t, 2 >& ext ) const
//******************************************************************************
//  Read coordinates of a number of mesh nodes from ExodusII file
//! \param[in] ext Extents of node ids whose coordinates to read
//! \param[inout] coord Unordered map of node coordinates associated to ids
//! \author J. Bakosi
//******************************************************************************
{
  auto num = ext[1] - ext[0] + 1;
  std::vector< tk::real > px( num ), py( num ), pz( num );

  ErrChk(
    ne_get_n_coord(
      m_inFile, static_cast<int64_t>(ext[0])+1, static_cast<int64_t>(num),
      px.data(), py.data(), pz.data() ) == 0,
      "Failed to read coordinates of nodes [" + std::to_string(ext[0]) +
      "..." + std::to_string(ext[1]) + "] from ExodusII file: " +
      m_filename );

  std::unordered_map< std::size_t, std::array< tk::real, 3 > > coord;
  for (std::size_t p=0; p<num; ++p)
    coord.emplace( ext[0]+p, std::array<tk::real,3>{{px[p],py[p],pz[p]}} );
  return coord;
}

std::size_t
ExodusIIMeshReader::readElemBlockIDs()
//******************************************************************************
//  Read element block IDs from ExodusII file
//! \return Total number of nodes in mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Read ExodusII file header
  auto nnode = readHeader();

  std::vector< int > eid( m_neblk );

  // Read element block ids
  ErrChk( ex_get_elem_blk_ids( m_inFile, eid.data()) == 0,
          "Failed to read element block ids from ExodusII file: " +
          m_filename );

  // Fill element block ID vector
  for (auto id : eid) {
    char eltype[MAX_STR_LENGTH+1];
    int nel, nnpe, nattr;

    // Read element block information
    ErrChk(
      ex_get_elem_block( m_inFile, id, eltype, &nel, &nnpe, &nattr ) == 0,
      "Failed to read element block information from ExodusII file: " +
      m_filename );

    // Store ExodusII element block ID
    m_eid.push_back( id );

    // Store ExodusII element block ID mapped to tk::ExoElemType enum, and
    // number of elements per block mapped to tk::ExoElemType enum
    if (nnpe == 4) {        // tetrahedra
      m_eidt[ static_cast<std::size_t>(ExoElemType::TET) ] = id;
      m_nel[ static_cast<std::size_t>(ExoElemType::TET) ] = nel;
    } else if (nnpe == 3) { // triangles
      m_eidt[ static_cast<std::size_t>(ExoElemType::TRI) ] = id;
      m_nel[ static_cast<std::size_t>(ExoElemType::TRI) ] = nel;
    }
  }

  return nnode;
}


void
ExodusIIMeshReader::readAllElements( UnsMesh& mesh )
//******************************************************************************
//  Read all element blocks and mesh connectivity from ExodusII file
//! \param[inout] mesh Unstructured mesh object to store mesh in
//! \author J. Bakosi
//******************************************************************************
{
  // Read element block ids
  readElemBlockIDs();

  for (auto id : m_eid) {
    char eltype[MAX_STR_LENGTH+1];
    int nel, nnpe, nattr;

    // Read element block information
    ErrChk(
      ex_get_elem_block( m_inFile, id, eltype, &nel, &nnpe, &nattr ) == 0,
      "Failed to read element block information from ExodusII file: " +
      m_filename );

    // Read element connectivity
    auto connectsize = static_cast< std::size_t >( nel*nnpe );
    if (nnpe == 4) {    // tetrahedra

      mesh.tettag().resize( connectsize, { 1 } );
      std::vector< int > inpoel( connectsize );
      ErrChk(
        ex_get_elem_conn( m_inFile, id, inpoel.data() ) == 0,
        "Failed to read " + std::string(eltype) + " element connectivity from "
        "ExodusII file: " + m_filename );
      for (auto n : inpoel)
        mesh.tetinpoel().push_back( static_cast< std::size_t >( n ) );

    } else if (nnpe == 3) {    // triangles

      mesh.tritag().resize( connectsize, { 1 } );
      std::vector< int > inpoel( connectsize );
      ErrChk(
        ex_get_elem_conn( m_inFile, id, inpoel.data() ) == 0,
        "Failed to read " + std::string(eltype) + " element connectivity from "
        "ExodusII file: " + m_filename );
      for (auto n : inpoel)
        mesh.triinpoel().push_back( static_cast< std::size_t >( n ) );

    }
  }

  // Shift node IDs to start from zero
  shiftToZero( mesh.triinpoel() );
  shiftToZero( mesh.tetinpoel() );
}

void
ExodusIIMeshReader::readElement( std::size_t id,
                                 tk::ExoElemType elemtype,
                                 std::vector< std::size_t >& conn ) const
//******************************************************************************
//  Read element connectivity of a single mesh cell from ExodusII file
//! \param[in] id Element id whose connectivity to read
//! \param[in] elemtype Element type
//! \param[inout] conn Connectivity vector to push to
//! \note Must be preceded by a call to readElemBlockIDs()
//! \author J. Bakosi
//******************************************************************************
{
  Assert( std::accumulate(begin(m_eidt), end(m_eidt), 0) != -m_nnpe.size(),
          "A call to ExodusIIMeshReader::readElement() must be preceded by a "
          "call to ExodusIIMeshReader::readElemBlockIDs()" );

  auto bid = static_cast< std::size_t >( elemtype );

  std::vector< int > c( m_nnpe[bid] );

  // Read element connectivity from file
  ErrChk(
    ne_get_n_elem_conn(
      m_inFile, m_eidt[bid], static_cast<int64_t>(id)+1, 1, c.data() ) == 0,
      "Failed to read element connectivity of element " + std::to_string( id ) +
      " from block " + std::to_string(m_eidt[bid]) + " from ExodusII file: " +
      m_filename );

  // Put in element connectivity using zero-based node indexing
  for (auto i : c) conn.push_back( static_cast<std::size_t>(i)-1 );
}

void
ExodusIIMeshReader::readElements( const std::array< std::size_t, 2 >& extent,
                                  tk::ExoElemType elemtype,
                                  std::vector< std::size_t >& conn ) const
//******************************************************************************
//  Read element connectivity of a single mesh cell from ExodusII file
//! \param[in] extent Extents of element ids whose connectivity to read
//! \param[in] elemtype Element type
//! \param[inout] conn Connectivity vector to push to
//! \note Must be preceded by a call to readElemBlockIDs()
//! \author J. Bakosi
//******************************************************************************
{
  Assert( std::accumulate(begin(m_eidt), end(m_eidt), 0) != -m_nnpe.size(),
          "A call to ExodusIIMeshReader::readElement() must be preceded by a "
          "call to ExodusIIMeshReader::readElemBlockIDs()" );

  auto bid = static_cast< std::size_t >( elemtype );

  auto num = extent[1] - extent[0];

  std::vector< int > c( num * m_nnpe[bid] );

  // Read element connectivity from file
  ErrChk(
    ne_get_n_elem_conn(
      m_inFile, m_eidt[bid], static_cast<int64_t>(extent[0])+1,
      static_cast<int64_t>(num), c.data() ) == 0,
      "Failed to read element connectivity of elements [" +
      std::to_string(extent[0]) + "..." + std::to_string(extent[1]) +
      "] from block " + std::to_string(m_eidt[bid]) + " from ExodusII file: " +
      m_filename );

  // Put in element connectivity using zero-based node indexing
  for (auto i : c) conn.push_back( static_cast<std::size_t>(i)-1 );
}

int
ExodusIIMeshReader::nel( tk::ExoElemType elemtype ) const
//******************************************************************************
//  Return number of elements in a mesh block in the ExodusII file
//! \param[in] elemtype Element type
//! \return Number of elements in elemtype given
//! \note Must be preceded by a call to readElemBlockIDs()
//! \author J. Bakosi
//******************************************************************************
{
  Assert( std::accumulate(begin(m_eidt),end(m_eidt),0) != -m_nnpe.size(),
          "A call to ExodusIIMeshReader::readElement() must be preceded by a "
          "call to ExodusIIMeshReader::readElemBlockIDs()" );

  return m_nel[ static_cast< std::size_t >( elemtype ) ];
}
