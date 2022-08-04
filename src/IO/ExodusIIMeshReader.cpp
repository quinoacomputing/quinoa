// *****************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     ExodusII mesh reader
  \details   ExodusII mesh reader class definition.
*/
// *****************************************************************************

#include <numeric>

#include "NoWarning/exodusII.hpp"

#include "ExodusIIMeshReader.hpp"
#include "ContainerUtil.hpp"
#include "Exception.hpp"
#include "UnsMesh.hpp"
#include "Reorder.hpp"

using tk::ExodusIIMeshReader;

ExodusIIMeshReader::ExodusIIMeshReader( const std::string& filename,
                                        int cpuwordsize,
                                        int iowordsize ) :
  m_filename( filename ),
  m_cpuwordsize( cpuwordsize ),
  m_iowordsize( iowordsize ),
  m_inFile( 0 ),
  m_nnode( 0 ),
  m_neblk( 0 ),
  m_neset( 0 ),
  m_from( 0 ),
  m_till( 0 ),
  m_blockid(),
  m_blockid_by_type( ExoNnpe.size() ),
  m_nel( ExoNnpe.size() ),
  m_elemblocks(),
  m_tri()
// *****************************************************************************
//  Constructor: open Exodus II file
//! \param[in] filename File to open as ExodusII file
//! \param[in] cpuwordsize Set CPU word size, see ExodusII documentation
//! \param[in] iowordsize Set I/O word size, see ExodusII documentation
// *****************************************************************************
{
  // Increase verbosity from ExodusII library in debug mode
  #ifndef NDEBUG
  ex_opts( EX_DEBUG | EX_VERBOSE );
  #endif

  float version;

  m_inFile = ex_open( filename.c_str(), EX_READ, &cpuwordsize, &iowordsize,
                      &version );

  ErrChk( m_inFile > 0, "Failed to open ExodusII file: " + filename );
}

ExodusIIMeshReader::~ExodusIIMeshReader() noexcept
// *****************************************************************************
//  Destructor
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
// *****************************************************************************
{
  readHeader( mesh );
  readAllElements( mesh );
  readAllNodes( mesh );
  readSidesetFaces( mesh.bface(), mesh.faceid() );
  readTimeValues( mesh.vartimes() );
  readNodeVarNames( mesh.nodevarnames() );
  readNodeScalars( mesh.vartimes().size(),
                   mesh.nodevarnames().size(),
                   mesh.nodevars() );
}

void
ExodusIIMeshReader::readGraph( UnsMesh& mesh )
// *****************************************************************************
//  Read only connectivity graph from file
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  readHeader( mesh );
  readAllElements( mesh );
}

void
ExodusIIMeshReader::readMeshPart(
  std::vector< std::size_t >& ginpoel,
  std::vector< std::size_t >& inpoel,
  std::vector< std::size_t >& triinp,
  std::unordered_map< std::size_t, std::size_t >& lid,
  tk::UnsMesh::Coords& coord,
  std::unordered_map< std::size_t, std::set< std::size_t > >& elemBlockId,
  int numpes, int mype )
// *****************************************************************************
//  Read a part of the mesh (graph and coordinates) from ExodusII file
//! \param[in,out] ginpoel Container to store element connectivity of this PE's
//!   chunk of the mesh (global ids)
//! \param[in,out] inpoel Container to store element connectivity with local
//!   node IDs of this PE's mesh chunk
//! \param[in,out] triinp Container to store triangle element connectivity
//!   (if exists in file) with global node indices
//! \param[in,out] lid Container to store global->local node IDs of elements of
//!   this PE's mesh chunk
//! \param[in,out] coord Container to store coordinates of mesh nodes of this
//!   PE's mesh chunk
//! \param[in,out] elemBlockId List of elements for each block-id.
//! \param[in] numpes Total number of PEs (default n = 1, for a single-CPU read)
//! \param[in] mype This PE (default m = 0, for a single-CPU read)
// *****************************************************************************
{
  Assert( mype < numpes, "Invalid input: PE id must be lower than NumPEs" );
  Assert( ginpoel.empty() && inpoel.empty() && lid.empty() &&
          coord[0].empty() && coord[1].empty() && coord[2].empty(),
          "Containers to store mesh must be empty" );

  // Read info on element blocks from ExodusII file
  readElemBlockIDs();
  // Get number of number of tetrahedron elements in file
  auto nel = nelem( tk::ExoElemType::TET );

  // Compute extents of element IDs of this PE's mesh chunk to read
  auto npes = static_cast< std::size_t >( numpes );
  auto pe = static_cast< std::size_t >( mype );
  auto chunk = nel / npes;
  m_from = pe * chunk;
  m_till = m_from + chunk;
  if (pe == npes-1) m_till += nel % npes;

  // Read tetrahedron connectivity between from and till
  readElements( {{m_from, m_till-1}}, tk::ExoElemType::TET, ginpoel );
  elemBlockId = m_elemInBlockId;

  // Compute local data from global mesh connectivity
  std::vector< std::size_t > gid;
  std::tie( inpoel, gid, lid ) = tk::global2local( ginpoel );

  // Read this PE's chunk of the mesh node coordinates from file
  coord = readCoords( gid );

  // Generate set of unique faces
  tk::UnsMesh::FaceSet faces;
  for (std::size_t e=0; e<ginpoel.size()/4; ++e)
    for (std::size_t f=0; f<4; ++f) {
      const auto& tri = tk::expofa[f];
      faces.insert( {{{ ginpoel[ e*4+tri[0] ],
                        ginpoel[ e*4+tri[1] ],
                        ginpoel[ e*4+tri[2] ] }}} );
    }

  // Read triangle element connectivity (all triangle blocks in file)
  auto ntri = nelem( tk::ExoElemType::TRI );
  if ( ntri !=0 ) readElements( {{0,ntri-1}}, tk::ExoElemType::TRI, triinp );

  // Keep triangles shared in (partially-read) tetrahedron mesh
  std::vector< std::size_t > triinp_own;
  std::size_t ltrid = 0;        // local triangle id
  for (std::size_t e=0; e<triinp.size()/3; ++e) {
    auto i = faces.find( {{ triinp[e*3+0], triinp[e*3+1], triinp[e*3+2] }} );
    if (i != end(faces)) {
      m_tri[e] = ltrid++;       // generate global->local triangle ids
      triinp_own.push_back( triinp[e*3+0] );
      triinp_own.push_back( triinp[e*3+1] );
      triinp_own.push_back( triinp[e*3+2] );
    }
  }
  triinp = std::move(triinp_own);
}

std::array< std::vector< tk::real >, 3 >
ExodusIIMeshReader::readCoords( const std::vector< std::size_t >& gid ) const
// *****************************************************************************
//  Read coordinates of a number of mesh nodes from ExodusII file
//! \param[in] gid Global node IDs whose coordinates to read
//! \return Vector of node coordinates read from file
// *****************************************************************************
{
  // Read node coordinates from file with global node IDs given in gid
  return readNodes( gid );
}

std::size_t
ExodusIIMeshReader::readHeader()
// *****************************************************************************
//  Read ExodusII header without setting mesh size
//! \return Number of nodes in mesh
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
  ErrChk( ndim == 3, "Need a 3D mesh from ExodusII file " + m_filename );

  m_neblk = static_cast< std::size_t >( neblk );
  m_neset = static_cast< std::size_t >( nelemset );

  return static_cast< std::size_t >( nnode );
}

void
ExodusIIMeshReader::readHeader( UnsMesh& mesh )
// *****************************************************************************
//  Read ExodusII header with setting mesh size
//! \param[in] mesh Unstructured mesh object
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
// *****************************************************************************
{
  mesh.x().resize( m_nnode );
  mesh.y().resize( m_nnode );
  mesh.z().resize( m_nnode );

  ErrChk( ex_get_coord( m_inFile, mesh.x().data(), mesh.y().data(),
                        mesh.z().data() ) == 0,
          "Failed to read coordinates from ExodusII file: " + m_filename );
}

void
ExodusIIMeshReader::readNode( std::size_t fid,
                              std::size_t mid,
                              std::vector< tk::real >& x,
                              std::vector< tk::real >& y,
                              std::vector< tk::real >& z ) const
// *****************************************************************************
//  Read coordinates of a single mesh node from ExodusII file
//! \param[in] fid Node id in file whose coordinates to read
//! \param[in] mid Node id in memory to which to put new cordinates
//! \param[in,out] x Vector of x coordinates to push to
//! \param[in,out] y Vector of y coordinates to push to
//! \param[in,out] z Vector of z coordinates to push to
// *****************************************************************************
{
  Assert( x.size() == y.size() && x.size() == z.size(), "Size mismatch" );
  Assert( mid < x.size() && mid < y.size() && mid < z.size(),
          "Indexing out of bounds" );

  readNode( fid, x[mid], y[mid], z[mid] );
}

void
ExodusIIMeshReader::readNode( std::size_t id,
                              std::array< tk::real, 3 >& coord ) const
// *****************************************************************************
//  Read coordinates of a single mesh node from ExodusII file
//! \param[in] id Node id whose coordinates to read
//! \param[in,out] coord Array of x, y, and z coordinates
// *****************************************************************************
{
  readNode( id, coord[0], coord[1], coord[2] );
}

void
ExodusIIMeshReader::readNode( std::size_t id,
                              tk::real& x,
                              tk::real& y,
                              tk::real& z ) const
// *****************************************************************************
// Read coordinates of a single mesh node from file
//! \param[in] id Node id whose coordinates to read
//! \param[in,out] x X coordinate to write to
//! \param[in,out] y Y coordinate to write to
//! \param[in,out] z Z coordinate to write to
// *****************************************************************************
{
  ErrChk(
    ex_get_partial_coord( m_inFile, static_cast<int64_t>(id)+1, 1,
                          &x, &y, &z ) == 0,
    "Failed to read coordinates of node " + std::to_string(id) +
    " from ExodusII file: " + m_filename );
}

std::array< std::vector< tk::real >, 3 >
ExodusIIMeshReader::readNodes( const std::vector< std::size_t >& gid ) const
// *****************************************************************************
//  Read coordinates of a number of mesh nodes from ExodusII file
//! \param[in] gid Node IDs whose coordinates to read
//! \return Mesh node coordinates
// *****************************************************************************
{
  std::vector< tk::real > px( gid.size() ), py( gid.size() ), pz( gid.size() );

  std::size_t i=0;
  for (auto g : gid) readNode( g, i++, px, py, pz );

  return {{ std::move(px), std::move(py), std::move(pz) }};
}

std::size_t
ExodusIIMeshReader::readElemBlockIDs()
// *****************************************************************************
//  Read element block IDs from ExodusII file
//! \return Total number of nodes in mesh
// *****************************************************************************
{
  // Read ExodusII file header
  auto nnode = readHeader();

  std::vector< int > bid( m_neblk );

  // Read element block ids
  ErrChk( ex_get_ids( m_inFile, EX_ELEM_BLOCK, bid.data()) == 0,
          "Failed to read element block ids from ExodusII file: " +
          m_filename );

  m_elemblocks.clear();
  m_nel.clear();
  m_nel.resize( ExoNnpe.size() );
  m_blockid_by_type.clear();
  m_blockid_by_type.resize( ExoNnpe.size() );

  // Fill element block ID vector
  for (auto id : bid) {
    char eltype[MAX_STR_LENGTH+1];
    int n, nnpe, nattr;

    // Read element block information
    ErrChk( ex_get_block( m_inFile, EX_ELEM_BLOCK, id, eltype, &n, &nnpe,
                          &nattr, nullptr, nullptr ) == 0,
      "Failed to read element block information from ExodusII file: " +
      m_filename );

    // Store ExodusII element block ID
    m_blockid.push_back( id );

    auto nel = static_cast< std::size_t >( n );

    // Store info on ExodusII element blocks
    if (nnpe == 4) {        // tetrahedra

      m_elemblocks.push_back( { ExoElemType::TET, nel } );
      auto e = static_cast< std::size_t >( ExoElemType::TET );
      m_blockid_by_type[ e ].push_back( id );
      m_nel[ e ].push_back( nel );
      Assert( m_blockid_by_type[e].size() == m_nel[e].size(), "Size mismatch" );

    } else if (nnpe == 3) { // triangles

      m_elemblocks.push_back( { ExoElemType::TRI, nel } );
      auto e = static_cast< std::size_t >( ExoElemType::TRI );
      m_blockid_by_type[ e ].push_back( id );
      m_nel[ e ].push_back( nel );
      Assert( m_blockid_by_type[e].size() == m_nel[e].size(), "Size mismatch" );

    }
  }

  return nnode;
}

void
ExodusIIMeshReader::readAllElements( UnsMesh& mesh )
// *****************************************************************************
//  Read all element blocks and mesh connectivity from ExodusII file
//! \param[inout] mesh Unstructured mesh object to store mesh in
// *****************************************************************************
{
  // Read element block ids
  readElemBlockIDs();

  for (auto id : m_blockid) {
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
                                  std::vector< std::size_t >& conn )
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
//!   The mesh block-wise element set is also updated.
// *****************************************************************************
{
  Assert( tk::sumsize(m_blockid_by_type) > 0,
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
          + std::to_string( ExoNnpe[ static_cast<std::size_t>(elemtype) ] )
          + " nodes per cell in file '" + m_filename + "': "
          + std::to_string( nelem( elemtype ) ) );

  auto e = static_cast< std::size_t >( elemtype );
  // List of number of elements of all blocks of element type requested
  const auto& nel = m_nel[e];
  // List of element block IDs for element type requested
  const auto& bid = m_blockid_by_type[e];

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
    std::vector< int > c( (r[1]-r[0]+1) * ExoNnpe[e] );
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

    // Store tet-elements under their respective mesh block ids
    if (elemtype == ExoElemType::TET) {
      for (std::size_t i=0; i<c.size()/ExoNnpe[e]; ++i) {
        auto& tetblk = m_elemInBlockId[bid[b]];
        tetblk.insert((inpoel.size()/ExoNnpe[e]) + i);
      }
    }

    inpoel.reserve( inpoel.size() + c.size() );
    std::move( begin(c), end(c), std::back_inserter(inpoel) );
  }

  Assert( inpoel.size() == (ext[1]-ext[0]+1)*ExoNnpe[e],
          "Failed to read element connectivity of elements [" +
          std::to_string(ext[0]) + "..." + std::to_string(ext[1]) + ") from "
          "ExodusII file: " + m_filename );

  // Put in element connectivity using zero-based node indexing
  for (auto& i : inpoel) --i;
  conn.reserve( conn.size() + inpoel.size() );
  std::move( begin(inpoel), end(inpoel), std::back_inserter(conn) );
}

void
ExodusIIMeshReader::readFaces( std::vector< std::size_t >& conn )
// *****************************************************************************
//  Read face connectivity of a number of boundary faces from ExodusII file
//! \param[inout] conn Connectivity vector to push to
//! \details This function reads in the total number of boundary faces,
//!   also called triangle-elements in the EXO2 file, and their connectivity.
// *****************************************************************************
{
  // Return quietly if no triangle elements in file
  if (nelem(tk::ExoElemType::TRI) == 0) return;

  // Read triangle boundary-face connectivity (all the triangle element block)
  readElements( {{0,nelem(tk::ExoElemType::TRI)-1}}, tk::ExoElemType::TRI,
                conn );
}

std::vector< std::size_t >
ExodusIIMeshReader::readNodemap()
// *****************************************************************************
//  Read local to global node-ID map from ExodusII file
//! \return node_map Vector mapping the local Exodus node-IDs to global node-IDs
//! \details The node-map is required to get the "Exodus-global" node-IDs from
//!   the "Exodus-internal" node-IDs, which are returned from the exodus APIs.
//!   The node-IDs in the exodus file are referred to as the "Exodus-global"
//!   node-IDs or "fileIDs" in Quinoa.
// *****************************************************************************
{
  // Read triangle boundary-face connectivity
  auto nnode = readElemBlockIDs();

  // Create array to store node-number map
  std::vector< int > node_map( nnode );

  // Read in the node number map to map the above nodes to the global node-IDs
  ErrChk( ex_get_id_map( m_inFile, EX_NODE_MAP, node_map.data() ) == 0,
          "Failed to read node map length from ExodusII file: " );

  std::vector< std::size_t > node_map1( nnode );

  for (std::size_t i=0; i<nnode; ++i)
  {
          node_map1[i] = static_cast< std::size_t >(node_map[i]-1);
  }

  return node_map1;
}

std::map< int, std::vector< std::size_t > >
ExodusIIMeshReader::readSidesetNodes()
// *****************************************************************************
//  Read node list of all side sets from ExodusII file
//! \return Node lists mapped to side set ids
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
      for (auto n : nodes) list.push_back( static_cast<std::size_t>(n-1) );
    }
  }

  return side;
}

void
ExodusIIMeshReader::readSidesetFaces(
  std::map< int, std::vector< std::size_t > >& bface,
  std::map< int, std::vector< std::size_t > >& faces )
// *****************************************************************************
//  Read side sets from ExodusII file
//! \param[in,out] bface Elem ids of side sets to read into
//! \param[in,out] faces Elem-relative face ids of tets of side sets
// *****************************************************************************
{
  // Read element block ids
  readElemBlockIDs();

  if (m_neset > 0) {
    // Read side set ids from file
    std::vector< int > ids( m_neset );
    ErrChk( ex_get_ids( m_inFile, EX_SIDE_SET, ids.data() ) == 0,
            "Failed to read side set ids from ExodusII file: " + m_filename );

    // Read all side sets from file
    for (auto i : ids) {
      int nface, nnode;

      // Read number of faces in side set
      ErrChk( ex_get_set_param( m_inFile, EX_SIDE_SET, i, &nface, &nnode ) == 0,
              "Failed to read side set " + std::to_string(i) + " parameters "
              "from ExodusII file: " + m_filename );

      Assert(nface > 0, "Number of faces = 0 in side set" + std::to_string(i));

      std::vector< int > exoelem( static_cast< std::size_t >( nface ) );
      std::vector< int > exoface( static_cast< std::size_t >( nface ) );

      // Read in file-internal element ids and relative face ids for side set
      ErrChk( ex_get_set( m_inFile, EX_SIDE_SET, i, exoelem.data(),
                          exoface.data() ) == 0,
              "Failed to read side set " + std::to_string(i) );

      // Store file-internal element ids of side set
      auto& elem = bface[i];
      elem.resize( exoelem.size() );
      std::size_t j = 0;
      for (auto e : exoelem) elem[j++] = static_cast< std::size_t >( e-1 );

      // Store zero-based relative face ids of side set
      auto& face = faces[i];
      face.resize( exoface.size() );
      j = 0;
      for (auto n : exoface) face[j++] = static_cast< std::size_t >( n-1 );

      Assert( std::all_of( begin(face), end(face),
                           [](std::size_t f){ return f<4; } ),
              "Relative face id of side set must be between 0 and 3" );
      Assert( elem.size() == face.size(), "Size mismatch" );
    }
  }
}

std::pair< tk::ExoElemType, std::size_t >
ExodusIIMeshReader::blkRelElemId( std::size_t id ) const
// *****************************************************************************
// Compute element-block-relative element id and element type
//! \param[in] id (ExodusII) file-internal element id
//! \return Element type the internal id points to and element id relative to
//!   cell-type
//! \details This function takes an internal element id, which in general can
//!   point to any element block in the ExodusII file and thus we do not know
//!   which element type a block contains. It then computes which cell type the
//!   id points to and computes the relative index for the given cell type. This
//!   is necessary because elements are read in from file by from potentially
//!   multiple blocks by cell type.
//! \note Must be preceded by a call to readElemBlockIDs()
// *****************************************************************************
{
  auto TRI = tk::ExoElemType::TRI;
  auto TET = tk::ExoElemType::TET;

  std::size_t e = 0;            // counts elements (independent of cell type)
  std::size_t ntri = 0;         // counts triangle elements
  std::size_t ntet = 0;         // counts tetrahedron elements

  for (const auto& b : m_elemblocks) {  // walk all element blocks in order
    e += b.second;                      // increment file-internal element id
    if (e > id) {                       // found element block for internal id
      if (b.first == TRI) {             // if triangle block
        return { TRI, id-ntet };        // return cell type and triangle id
      } else if (b.first == TET) {      // if tetrahedron block
        return { TET, id-ntri };        // return cell type and tetrahedron id
      }
    }
    // increment triangle and tetrahedron elements independently
    if (b.first == TRI)
      ntri += b.second;
    else if (b.first == TET)
      ntet += b.second;
  }

  Throw( " Exodus internal element id not found" );
}

std::vector< std::size_t >
ExodusIIMeshReader::triinpoel(
  std::map< int, std::vector< std::size_t > >& belem,
  const std::map< int, std::vector< std::size_t > >& faces,
  const std::vector< std::size_t >& ginpoel,
  const std::vector< std::size_t >& triinp ) const
// *****************************************************************************
//  Generate triangle face connectivity for side sets
//! \param[in,out] belem File-internal elem ids of side sets
//! \param[in] faces Elem-relative face ids of side sets
//! \param[in] ginpoel Tetrahedron element connectivity with global nodes
//! \param[in] triinp Triangle element connectivity with global nodes
//!   (if exists in file)
//! \return Triangle face connectivity with global node IDs of side sets
//! \details This function takes lists of file-internal element ids (in belem)
//!   for side sets and does two things: (1) generates face connectivity (with
//!   global node IDs) for side sets, and (2) converts the (ExodusII)
//!   file-internal element IDs to face ids so that they can be used to index
//!   into the face connectivity. The IDs in belem are modified and the face
//!   connectivity (for boundary faces only) is returned.
//! \note Must be preceded by a call to readElemBlockIDs()
// *****************************************************************************
{
  Assert( !(m_from == 0 && m_till == 0),
          "Lower and upper tetrahedron id bounds must not both be zero" );

  // This will contain one of our final results: face (triangle) connectivity
  // for the side sets only. The difference between bnd_triinpoel and triinpoel
  // is that triinpoel is a triangle element connectivity, independent of side
  // sets, while bnd_triinpoel is a triangle connectivity only for side sets.
  std::vector< std::size_t > bnd_triinpoel;

  // Storage for boundary face lists for each side set on this PE
  std::map< int, std::vector< std::size_t > > belem_own;

  std::size_t f = 0;            // counts all faces
  for (auto& ss : belem) {      // for all side sets

    // insert side set id into new map
    auto& b = belem_own[ ss.first ];
    // get element-relative face ids for side set
    const auto& face = tk::cref_find( faces, ss.first );
    std::size_t s = 0;          // counts side set faces
    for (auto& i : ss.second) { // for all faces on side set

      // compute element-block-relative element id and element type
      auto r = blkRelElemId( i );

      // extract boundary face connectivity based on element type
      bool localface = false;
      if (r.first == tk::ExoElemType::TRI) {

        auto t = m_tri.find(r.second);
        if (t != end(m_tri)) {  // only if triangle id exists on this PE
          Assert( t->second < triinp.size()/3,
                  "Indexing out of triangle connectivity" );
          // generate triangle (face) connectivity using global node ids
          bnd_triinpoel.push_back( triinp[ t->second*3 + 0 ] );
          bnd_triinpoel.push_back( triinp[ t->second*3 + 1 ] );
          bnd_triinpoel.push_back( triinp[ t->second*3 + 2 ] );
          localface = true;
        }

      } else if (r.first == tk::ExoElemType::TET) {

        if (r.second >= m_from && r.second < m_till) {  // if tet is on this PE
          auto t = r.second - m_from;
          Assert( t < ginpoel.size()/4,
                  "Indexing out of tetrahedron connectivity" );
          // get ExodusII face-node numbering for side sets, see ExodusII
          // manual figure on "Sideset side Numbering"
          const auto& tri = tk::expofa[ face[s] ];
          // generate triangle (face) connectivity using global node ids, note
          // the switched node order, 0,2,1, as lpofa is different from expofa
          bnd_triinpoel.push_back( ginpoel[ t*4 + tri[0] ] );
          bnd_triinpoel.push_back( ginpoel[ t*4 + tri[1] ] );
          bnd_triinpoel.push_back( ginpoel[ t*4 + tri[2] ] );
          localface = true;
        }

      }

      ++s;

      // generate PE-local face id for side set (this is to be used to index
      // into triinpoel)
      if (localface) b.push_back( f++ );
    }

    // if no faces on this side set (on this PE), remove side set id
    if (b.empty()) belem_own.erase( ss.first );
  }

  belem = std::move(belem_own);

  return bnd_triinpoel;
}

void
ExodusIIMeshReader::readNodeVarNames( std::vector< std::string >& nv ) const
// *****************************************************************************
//  Read the names of nodal output variables from ExodusII file
//! \param[in,out] nv Nodal variable names
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

  int numvars = 0;

  ErrChk(
    ex_get_variable_param( m_inFile, EX_NODE_BLOCK, &numvars ) == 0,
    "Failed to read nodal output variable parameters from ExodusII file: " +
    m_filename );

  if (numvars) {

    char* names[ static_cast< std::size_t >( numvars ) ];
    for (int i=0; i<numvars; ++i)
      names[i] = static_cast<char*>( calloc((MAX_STR_LENGTH+1), sizeof(char)) );

    ErrChk( ex_get_variable_names( m_inFile,
                                   EX_NODAL,
                                   numvars,
                                   names ) == 0,
            "Failed to read nodal variable names from ExodusII file: " +
            m_filename );

    nv.resize( static_cast< std::size_t >( numvars ) );
    std::size_t i = 0;
    for (auto& n : nv) n = names[ i++ ];

  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

void
ExodusIIMeshReader::readTimeValues( std::vector< tk::real >& tv ) const
// *****************************************************************************
//  Read time values from ExodusII file
//! \param[in] tv Vector of time values at which field data is saved
// *****************************************************************************
{
  auto num_time_steps =
    static_cast< std::size_t >( ex_inquire_int( m_inFile, EX_INQ_TIME ) );

  if (num_time_steps) {
    tv.resize( num_time_steps, 0.0 );
    ErrChk( ex_get_all_times( m_inFile, tv.data() ) == 0,
             "Failed to read time values from ExodusII file: " + m_filename );
  }
}

void
ExodusIIMeshReader::readNodeScalars(
  std::size_t ntime,
  std::size_t nvar,
  std::vector< std::vector< std::vector< tk::real > > >& var ) const
// *****************************************************************************
//  Read node scalar fields from ExodusII file
//! \param[in] ntime Number of time steps to read
//! \param[in] nvar Number of variables to read
//! \param[in] var Vector of nodal variables to read to: inner vector: nodes,
//!   middle vector: (physics) variable, outer vector: time step
// *****************************************************************************
{
  var.resize( ntime );
  for (auto& v : var) {
    v.resize( nvar );
    for (auto& n : v) n.resize( m_nnode );
  }

  for (std::size_t t=0; t<var.size(); ++t) {
    for (std::size_t id=0; id<var[t].size(); ++id) {
      ErrChk( ex_get_var( m_inFile,
                          static_cast< int >( t+1 ),
                          EX_NODAL,
                          static_cast< int >( id+1 ),
                          1,
                          static_cast< int64_t >( var[t][id].size() ),
                          var[t][id].data() ) == 0,
              "Failed to read node scalar from ExodusII file: " + m_filename );
    }
  }
}

std::size_t
ExodusIIMeshReader::nelem( tk::ExoElemType elemtype ) const
// *****************************************************************************
//  Return number of elements in all mesh blocks for a given elem type in file
//! \param[in] elemtype Element type
//! \return Number of elements in all blocks for the elem type
//! \note Must be preceded by a call to readElemBlockIDs()
// *****************************************************************************
{
  auto e = static_cast< std::size_t >( elemtype );
  return std::accumulate( m_nel[e].cbegin(), m_nel[e].cend(), 0u );
}
