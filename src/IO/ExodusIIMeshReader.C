// *****************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     ExodusII mesh reader
  \details   ExodusII mesh reader class definition. Currently, this is a bare
     minimum functionality to interface with the ExodusII reader. It only reads
     3D meshes and only triangle and tetrahedron elements.
*/
// *****************************************************************************

#include <cstdint>
#include <cstdio>
#include <string>
#include <numeric>
#include <unordered_map>

#include "NoWarning/exodusII.h"

#include "ExodusIIMeshReader.h"
#include "ContainerUtil.h"
#include "Exception.h"
#include "UnsMesh.h"
#include "Reorder.h"

using tk::ExodusIIMeshReader;

ExodusIIMeshReader::ExodusIIMeshReader( const std::string& filename,
                                        int cpuwordsize,
                                        int iowordsize ) :
  m_filename( filename ),
  m_inFile( 0 ),
  m_nnode( 0 ),
  m_neblk( 0 ),
  m_neset( 0 ),
  m_eid(),
  m_eidt( m_nnpe.size() ),
  m_nel( m_nnpe.size() )
// *****************************************************************************
//  Constructor: open Exodus II file
//! \param[in] filename File to open as ExodusII file
//! \param[inout] mesh Unstructured mesh object to load the data to
//! \param[in] cpuwordsize Set CPU word size, see ExodusII documentation
//! \param[in] iowordsize Set I/O word size, see ExodusII documentation
//! \author J. Bakosi
// *****************************************************************************
{
  float version;

  m_inFile = ex_open( filename.c_str(), EX_READ, &cpuwordsize, &iowordsize,
                      &version );

  ErrChk( m_inFile > 0, "Failed to open ExodusII file: " + filename );
}

ExodusIIMeshReader::~ExodusIIMeshReader() noexcept
// *****************************************************************************
//  Destructor
//! \author J. Bakosi
// *****************************************************************************
{
  if ( ex_close(m_inFile) < 0 )
    printf( ">>> WARNING: Failed to close ExodusII file: %s\n",
            m_filename.c_str() );
}

void
ExodusIIMeshReader::readMesh( UnsMesh& mesh )
// *****************************************************************************
//  Read ExodusII mesh file
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  readHeader( mesh );
  readAllElements( mesh );
  readAllNodes( mesh );
}

void
ExodusIIMeshReader::readGraph( UnsMesh& mesh )
// *****************************************************************************
//  Read only connectivity graph from file
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  readHeader( mesh );
  readAllElements( mesh );
}

std::size_t
ExodusIIMeshReader::readHeader()
// *****************************************************************************
//  Read ExodusII header without setting mesh size
//! \return Number of nodes in mesh
//! \author J. Bakosi
// *****************************************************************************
{
  char title[MAX_LINE_LENGTH+1];
  int ndim, n, nnodeset, nelemset, nnode, neblk;

  ErrChk(
    ex_get_init( m_inFile, title, &ndim, &nnode, &n, &neblk, &nnodeset,
                 &nelemset ) == 0,
    "Failed to read header from ExodusII file: " + m_filename );

  ErrChk( nnode > 0,
          "Number of nodes read from ExodusII file must be larger than zero" );
  ErrChk( neblk > 0,
          "Number of element blocks read from ExodusII file must be larger "
          "than zero" );
  ErrChk( ndim == 3, "Need a 3D mesh from ExodusII file " + m_filename);

  m_neblk = static_cast< std::size_t >( neblk );
  m_neset = static_cast< std::size_t >( nelemset );

  return static_cast< std::size_t >( nnode );
}

void
ExodusIIMeshReader::readHeader( UnsMesh& mesh )
// *****************************************************************************
//  Read ExodusII header with setting mesh size
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  // Read ExodusII file header and set mesh graph size
  mesh.size() = m_nnode = static_cast< std::size_t >( readHeader() );
}

void
ExodusIIMeshReader::readAllNodes( UnsMesh& mesh ) const
// *****************************************************************************
//  Read all node coordinates from ExodusII file
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  mesh.x().resize( m_nnode );
  mesh.y().resize( m_nnode );
  mesh.z().resize( m_nnode );

  ErrChk( ex_get_coord( m_inFile, mesh.x().data(), mesh.y().data(),
                        mesh.z().data() ) == 0,
          "Failed to read coordinates from ExodusII file: " + m_filename );
}

std::array< std::vector< tk::real >, 3 >
ExodusIIMeshReader::readNodes( const std::array< std::size_t, 2 >& ext ) const
// *****************************************************************************
//  Read coordinates of a number of mesh nodes from ExodusII file
//! \param[in] ext Extents of element ids whose connectivity to read, both
//!   inclusive
//! \return Mesh node coordinates
//! \author J. Bakosi
// *****************************************************************************
{
  auto num = ext[1] - ext[0] + 1;
  std::vector< tk::real > px( num ), py( num ), pz( num );

  ErrChk(
    ex_get_partial_coord(
      m_inFile, static_cast<int64_t>(ext[0])+1, static_cast<int64_t>(num),
      px.data(), py.data(), pz.data() ) == 0,
      "Failed to read coordinates of nodes [" + std::to_string(ext[0]) +
      "..." + std::to_string(ext[1]) + "] from ExodusII file: " +
      m_filename );

  return {{ std::move(px), std::move(py), std::move(pz) }};
}

std::size_t
ExodusIIMeshReader::readElemBlockIDs()
// *****************************************************************************
//  Read element block IDs from ExodusII file
//! \return Total number of nodes in mesh
//! \author J. Bakosi
// *****************************************************************************
{
  // Read ExodusII file header
  auto nnode = readHeader();

  std::vector< int > eid( m_neblk );

  // Read element block ids
  ErrChk( ex_get_ids( m_inFile, EX_ELEM_BLOCK, eid.data()) == 0,
          "Failed to read element block ids from ExodusII file: " +
          m_filename );

  // Fill element block ID vector
  for (auto id : eid) {
    char eltype[MAX_STR_LENGTH+1];
    int n, nnpe, nattr;

    // Read element block information
    ErrChk( ex_get_block( m_inFile, EX_ELEM_BLOCK, id, eltype, &n, &nnpe,
                          &nattr, nullptr, nullptr ) == 0,
      "Failed to read element block information from ExodusII file: " +
      m_filename );

    // Store ExodusII element block ID
    m_eid.push_back( id );

    // Store ExodusII element block IDs mapped to elem type, and add up the
    // number of elements per for each elem type
    if (nnpe == 4) {        // tetrahedra
      auto e = static_cast< std::size_t >( ExoElemType::TET );
      m_eidt[ e ].push_back( id );
      m_nel[ e ].push_back( static_cast< std::size_t >( n ) );
      Assert( m_eidt[e].size() == m_nel[e].size(), "Size mismatch" );
    } else if (nnpe == 3) { // triangles
      auto e = static_cast< std::size_t >( ExoElemType::TRI );
      m_eidt[ e ].push_back( id );
      m_nel[ e ].push_back( static_cast< std::size_t >( n ) );
      Assert( m_eidt[e].size() == m_nel[e].size(), "Size mismatch" );
    }
  }

  return nnode;
}


void
ExodusIIMeshReader::readAllElements( UnsMesh& mesh )
// *****************************************************************************
//  Read all element blocks and mesh connectivity from ExodusII file
//! \param[inout] mesh Unstructured mesh object to store mesh in
//! \author J. Bakosi
// *****************************************************************************
{
  // Read element block ids
  readElemBlockIDs();

  for (auto id : m_eid) {
    char eltype[MAX_STR_LENGTH+1];
    int nel, nnpe, nattr;

    // Read element block information
    ErrChk( ex_get_block( m_inFile, EX_ELEM_BLOCK, id, eltype, &nel, &nnpe,
                          &nattr, nullptr, nullptr ) == 0,
      "Failed to read element block information from ExodusII file: " +
      m_filename );

    // Read element connectivity
    auto connectsize = static_cast< std::size_t >( nel*nnpe );
    if (nnpe == 4) {    // tetrahedra

      std::vector< int > inpoel( connectsize );
      ErrChk( ex_get_conn( m_inFile, EX_ELEM_BLOCK, id, inpoel.data(),
                           nullptr, nullptr ) == 0,
        "Failed to read " + std::string(eltype) + " element connectivity from "
        "ExodusII file: " + m_filename );
      for (auto n : inpoel)
        mesh.tetinpoel().push_back( static_cast< std::size_t >( n ) );

    } else if (nnpe == 3) {    // triangles

      std::vector< int > inpoel( connectsize );
      ErrChk( ex_get_conn( m_inFile, EX_ELEM_BLOCK, id, inpoel.data(),
                           nullptr, nullptr ) == 0,
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
ExodusIIMeshReader::readElements( const std::array< std::size_t, 2 >& ext,
                                  tk::ExoElemType elemtype,
                                  std::vector< std::size_t >& conn ) const
// *****************************************************************************
//  Read element connectivity of a number of mesh cells from ExodusII file
//! \param[in] ext Extents of element IDs whose connectivity to read:
//!   [from...till), using zero-based element IDs, where 'from' >=0, inclusive
//!   and 'till < 'maxelements', where 'maxelements' is the total number of
//!   elements of all element blocks in the file of the requested cell type.
//!   Note that 'maxelements' can be queried by nelem().
//! \param[in] elemtype Element type
//! \param[inout] conn Connectivity vector to push to
//! \note Must be preceded by a call to readElemBlockIDs()
//! \details This function takes the extents of element IDs in a zero-based
//!   fashion. These input extents can be thought of "absolute" extents that
//!   denote lowest and the largest-1 element IDs to be read from file.
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( tk::sumsize(m_eidt) > 0,
          "A call to this function must be preceded by a call to "
          "ExodusIIMeshReader::readElemBlockIDs()" );
  Assert( ext[0] <= ext[1] &&
          ext[0] < nelem(elemtype) &&
          ext[1] < nelem(elemtype),
          "Invalid element ID extents. Of the requested extents [from...till), "
          "'from' must be lower than or equal to 'till', and they must be in "
          "the range [0...maxelements), where 'maxelements' is the total "
          "number of elements of all element blocks in the file of the "
          "requested cell type. Requested element ID extents: ["
          + std::to_string(ext[0]) + "..." + std::to_string(ext[1])
          + "), 'maxelements' of cell type with "
          + std::to_string( m_nnpe[ static_cast<std::size_t>(elemtype) ] )
          + " nodes per cell in file '" + m_filename + "': "
          + std::to_string( nelem( elemtype ) ) );

  auto e = static_cast< std::size_t >( elemtype );
  // List of number of elements of all blocks of element type requested
  const auto& nel = m_nel[e];
  // List of element block IDs for element type requested
  const auto& bid = m_eidt[e];

  // Compute lower and upper element block ids to read from based on extents
  std::size_t lo_bid = 0, hi_bid = 0, offset = 0;
  for (std::size_t b=0; b<nel.size(); ++b) {
    std::size_t lo = offset;                    // lo (min) elem ID in block
    std::size_t hi = offset + nel[b] - 1;       // hi (max) elem ID in block
    if (ext[0] >= lo && ext[0] <= hi) lo_bid = b;
    if (ext[1] >= lo && ext[1] <= hi) hi_bid = b;
    offset += nel[b];
  }

  Assert( lo_bid < nel.size() && lo_bid < bid.size(),
          "Invalid start block ID" );
  Assert( hi_bid < nel.size() && hi_bid < bid.size(),
          "Invalid end block ID" );

  // Compute relative extents based on absolute ones for each block to read from
  std::vector< std::array< std::size_t, 2 > > rext;
  offset = 0;
  for (std::size_t b=0; b<lo_bid; ++b) offset += nel[b];
  for (std::size_t b=lo_bid; b<=hi_bid; ++b) {
    std::size_t lo = offset;
    std::size_t hi = offset + nel[b] - 1;
    std::size_t le = 1, he = nel[b];
    if (ext[0] >= lo && ext[0] <= hi) le = ext[0] - lo + 1;
    if (ext[1] >= lo && ext[1] <= hi) he = ext[1] - lo + 1;
    Assert( le >= 1 && le <= nel[b] && he >= 1 && he <= nel[b],
            "Relative index out of block" );
    rext.push_back( {{ le, he }} );
    offset += nel[b];
  }

  Assert( std::accumulate(
            std::next(rext.cbegin()), rext.cend(), rext[0][1]-rext[0][0]+1,
            []( std::size_t n, const std::array< std::size_t, 2 >& r )
            { return n + r[1] - r[0] + 1; }
          ) == ext[1]-ext[0]+1,
          "Total number of elements to read incorrect, requested extents: " +
          std::to_string(ext[0]) + " ... " + std::to_string(ext[1]) );

  std::vector< int > inpoel;

  // Read element connectivity from file
  std::size_t B = 0;
  for (auto b=lo_bid; b<=hi_bid; ++b, ++B) {
    const auto& r = rext[B];
    std::vector< int > c( (r[1]-r[0]+1) * m_nnpe[e] );
    ErrChk( ex_get_partial_conn( m_inFile,
                                 EX_ELEM_BLOCK,
                                 bid[b],
                                 static_cast< int64_t >( r[0] ),
                                 static_cast< int64_t >( r[1]-r[0]+1 ),
                                 c.data(),
                                 nullptr,
                                 nullptr ) == 0,
            "Failed to read element connectivity of elements [" +
            std::to_string(r[0]) + "..." + std::to_string(r[1]) +
            "] from element block " + std::to_string(bid[b]) + " in ExodusII "
            "file: " + m_filename );
    inpoel.reserve( inpoel.size() + c.size() );
    std::move( begin(c), end(c), std::back_inserter(inpoel) );
  }

  Assert( inpoel.size() == (ext[1]-ext[0]+1)*4,
          "Failed to read element connectivity of elements [" +
          std::to_string(ext[0]) + "..." + std::to_string(ext[1]) + ") from "
          "ExodusII file: " + m_filename );

  // Put in element connectivity using zero-based node indexing
  for (auto& i : inpoel) --i;
  conn.reserve( conn.size() + inpoel.size() );
  std::move( begin(inpoel), end(inpoel), std::back_inserter(conn) );
}

std::map< int, std::vector< std::size_t > >
ExodusIIMeshReader::readSidesets()
// *****************************************************************************
//  Read node list of all side sets from ExodusII file
//! \return Elem and side lists mapped to side set ids
//! \author J. Bakosi
// *****************************************************************************
{
  // Read ExodusII file header (fills m_neset)
  readHeader();

  // Node lists mapped to side set ids
  std::map< int, std::vector< std::size_t > > side;

  if (m_neset > 0) {
    // Read all side set ids from file
    std::vector< int > ids( m_neset );
    ErrChk( ex_get_ids( m_inFile, EX_SIDE_SET, ids.data() ) == 0,
            "Failed to read side set ids from ExodusII file: " + m_filename );
    // Read in node list for all side sets
    for (auto i : ids) {
      int nface, nnode;
      // Read number of faces and number of distribution factors in side set i
      ErrChk( ex_get_set_param( m_inFile, EX_SIDE_SET, i, &nface, &nnode ) == 0,
              "Failed to read side set " + std::to_string(i) + " parameters "
              "from ExodusII file: " + m_filename );
      // Read number of nodes in side set i (overwrite nnode)
      ErrChk( ex_get_side_set_node_list_len( m_inFile, i, &nnode ) == 0,
              "Failed to read side set " + std::to_string(i) + " node list "
              "length from ExodusII file: " + m_filename );
      Assert(nnode > 0, "Number of nodes = 0 in side set" + std::to_string(i));
      std::vector< int > df( static_cast< std::size_t >( nface ) );
      std::vector< int > nodes( static_cast< std::size_t >( nnode ) );
      // Read in node list for side set i
      ErrChk( ex_get_side_set_node_list( m_inFile, i, df.data(), nodes.data() )
                == 0, "Failed to read node list of side set " +
                      std::to_string(i) + " from ExodusII file: " +
                      m_filename );
      // Make node list unique
      tk::unique( nodes );
      // Store 0-based node ID list as std::size_t vector instead of ints
      auto& list = side[ i ];
      for (auto&& n : nodes) list.push_back( static_cast<std::size_t>(n-1) );
    }
  }

  return side;
}

std::size_t
ExodusIIMeshReader::nelem( tk::ExoElemType elemtype ) const
// *****************************************************************************
//  Return number of elements in all mesh blocks for a given elem type in file
//! \param[in] elemtype Element type
//! \return Number of elements in all blocks for the elem type
//! \note Must be preceded by a call to readElemBlockIDs()
//! \author J. Bakosi
// *****************************************************************************
{
  auto e = static_cast< std::size_t >( elemtype );
  std::size_t sum = 0;
  for (auto n : m_nel[e]) sum += n;
  return sum;
}
