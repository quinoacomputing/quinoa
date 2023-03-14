// *****************************************************************************
/*!
  \file      src/IO/ExodusIIMeshWriter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     ExodusII mesh-based data writer
  \details   ExodusII mesh-based data writer class definition.
*/
// *****************************************************************************

#include <numeric>

#include "NoWarning/exodusII.hpp"

#include "ExodusIIMeshWriter.hpp"
#include "Exception.hpp"
#include "UnsMesh.hpp"

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
// *****************************************************************************
{
  // Increase verbosity from ExodusII library in debug mode
  #ifndef NDEBUG
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
// *****************************************************************************
{
  if ( ex_close(m_outFile) < 0 )
    printf( ">>> WARNING: Failed to close ExodusII file: %s\n",
            m_filename.c_str() );
}

void
ExodusIIMeshWriter::writeMesh( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write ExodusII mesh file taking a tk::UnsMesh object
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  writeHeader( mesh );
  writeNodes( mesh );
  writeElements( mesh );
  writeSidesets( mesh );
  writeNodesets( mesh );
  writeTimeValues( mesh.vartimes() );
  writeNodeVarNames( mesh.nodevarnames() );
  writeElemVarNames( mesh.elemvarnames() );
  writeNodeScalars( mesh.nodevars() );
  writeElemScalars( mesh.elemvars() );
}

void
ExodusIIMeshWriter::writeMesh(
  const std::vector< std::size_t >& tetinp,
  const UnsMesh::Coords& coord,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinp ) const
// *****************************************************************************
//  Write ExodusII mesh file taking inputs to a tk::UnsMesh object
//! \param[in] tetinp Tetrahedron element connectivity
//! \param[in] coord Node coordinates
//! \param[in] bface Boundary face ids for each side set
//! \param[in] triinp Triangle face connectivity (for all faces in bface)
// *****************************************************************************
{
  // Fill element-relative face ids (0,1,2,3) for all side sets with 0 (will
  // use triangles as face elements for side sets)
  std::map< int, std::vector< std::size_t > > faceid;
  for (const auto& s : bface) faceid[s.first].resize( s.second.size(), 0 );

  // Generate file-internal Exodus (triangle face) element ids for all faces of
  // all side sets. bface_exo:: face ids for each side set, triinpoel: triangle
  // face connectivity for all side sets.
  std::map< int, std::vector< std::size_t > > bface_exo;
  std::vector< std::size_t > triinpoel( triinp.size() );
  // Generate/start exodus-file-face ids from max number of tetrahedra because
  // tet elem blocks will be written out first
  std::size_t exo_faceid = tetinp.size() / 4;
  for (const auto& s : bface) {
    auto& b = bface_exo[ s.first ];
    b.resize( s.second.size() );
    std::size_t j = 0;  // side-set-relative face id
    for (auto f : s.second) {   // for all faces on side set s.first
      b[ j++ ] = exo_faceid;    // generate exo-file face id
      auto k = exo_faceid - tetinp.size()/4;
      // copy over triangle connectivity in order
      triinpoel[ k*3+0 ] = triinp[ f*3+0 ];
      triinpoel[ k*3+1 ] = triinp[ f*3+1 ];
      triinpoel[ k*3+2 ] = triinp[ f*3+2 ];
      ++exo_faceid;
    }
  }

  // Write mesh
  writeMesh( tk::UnsMesh( tetinp, coord, bface_exo, triinpoel, faceid ) );
}

void
ExodusIIMeshWriter::writeMesh(
  const std::vector< std::size_t >& tetinp,
  const UnsMesh::Coords& coord,
  const std::map< int, std::vector< std::size_t > >& bnode ) const
// *****************************************************************************
//  Write ExodusII mesh file taking inputs to a tk::UnsMesh object
//! \param[in] tetinp Tetrahedron element connectivity
//! \param[in] coord Node coordinates
//! \param[in] bnode Boundary node ids for each side set
// *****************************************************************************
{
  writeMesh( tk::UnsMesh( tetinp, coord, bnode ) );
}

void
ExodusIIMeshWriter::writeHeader( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write ExodusII header
//! \param[in] mesh Unstructured mesh object
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
                 static_cast< int64_t >( mesh.bnode().size() ),
                 static_cast< int64_t >( mesh.bface().size() )
               ) == 0,
    "Failed to write header to file: " + m_filename );
}

void
ExodusIIMeshWriter::writeHeader( const char* title,
                                 int64_t ndim,
                                 int64_t nnodes,
                                 int64_t nelem,
                                 int64_t nblk,
                                 int64_t node_set,
                                 int64_t side_set ) const
// *****************************************************************************
//  Write ExodusII header
//! \param[in] title ExodusII file header 'title'
//! \param[in] ndim Number of spatial dimensions in ExodusII file
//! \param[in] nnodes Number of mesh nodes in ExodusII file
//! \param[in] nelem Number of mesh elements in ExodusII file
//! \param[in] nblk Number of mesh element blocks in ExodusII file
//! \param[in] node_set Number of node sets in ExodusII file
//! \param[in] side_set Number of side sets in ExodusII file
// *****************************************************************************
{
  ErrChk(
    ex_put_init( m_outFile, title, ndim, nnodes, nelem, nblk, 
                 node_set, side_set) == 0,
    "Failed to write header to file: " + m_filename );
}

void
ExodusIIMeshWriter::writeNodes( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write node coordinates to ExodusII file
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  ErrChk( ex_put_coord( m_outFile, mesh.x().data(), mesh.y().data(),
                        mesh.z().data() ) == 0,
          "Failed to write coordinates to ExodusII file: " + m_filename );
}

void
ExodusIIMeshWriter::writeNodes( const std::vector< tk::real >& x,
                                const std::vector< tk::real >& y,
                                const std::vector< tk::real >& z ) const
// *****************************************************************************
//  Write node coordinates to ExodusII file without Mesh
//! \param[in] x coordinates of mesh nodes
//! \param[in] y coordinates of mesh nodes
//! \param[in] z coordinates of mesh nodes
// *****************************************************************************
{
  ErrChk( ex_put_coord( m_outFile, x.data(), y.data(), z.data() ) == 0,
          "Failed to write coordinates to ExodusII file: " + m_filename );
}

void
ExodusIIMeshWriter::writeElements( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write element connectivity to ExodusII file
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  int elclass = 0;

  // For meshes that have no triangle element block (only tetrahedra), keeping
  // the order as tets first followed by triangles allows keeping the tet ids
  // associated to side sets the same when adding the missing triangle elements
  // by meshconv. Hence this order should be kept as tets first triangles next.

  writeElemBlock( elclass, 4, "TETRAHEDRA", mesh.tetinpoel() );
  writeElemBlock( elclass, 3, "TRIANGLES", mesh.triinpoel() );
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
//! \param[in] inpoel Element connectivity
// *****************************************************************************
{
  if (inpoel.empty()) return;

  // increase number of element classes in file
  ++elclass;

  // Compute number of edges and number of faces for triangles and tetrahedra
  int nedge = 0, nface = 0;
  if (nnpe == 4) {
    nedge = 6;
    nface = 4;
  } else if (nnpe == 3) {
    nedge = 3;
    nface = 1;
  } else Throw( "Write ExodusII element block does not support elements with "
                + std::to_string(nnpe) + " nodes" );

  // Write element block information
  ErrChk(
    ex_put_block( m_outFile,            // exo file handle
                  EX_ELEM_BLOCK,        // block type: elem block
                  elclass,              // element class id
                  eltype.c_str(),       // element block description
                  static_cast< int64_t >( inpoel.size() ) / nnpe, // nof cells
                  nnpe, // number of nodes per element
                  nedge,// number of edges per element
                  nface,// number of faces per element
                  0     // number of attributes per element
                ) == 0,
    "Failed to write " + eltype + " element block to ExodusII file: " +
    m_filename );

  // Write element connectivity with 1-based node ids
  std::vector< int > inp( inpoel.size() );
  std::size_t i = 0;
  for (auto p : inpoel) inp[ i++ ] = static_cast< int >( p+1 );
  ErrChk( ex_put_conn( m_outFile, EX_ELEM_BLOCK, elclass, inp.data(),
                       nullptr, nullptr ) == 0,
          "Failed to write " + eltype + " element connectivity to ExodusII "
          "file: " + m_filename );
}

void
ExodusIIMeshWriter::writeSidesets( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write side sets and their face connectivity to ExodusII file
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  // Write all side sets face list and connectivity in mesh
  for (const auto& s : mesh.bface()) {
    // Write side set parameters
    ErrChk( ex_put_set_param( m_outFile, EX_SIDE_SET, s.first,
                              static_cast<int64_t>(s.second.size()), 0 ) == 0,
      "Failed to write side set parameters to ExodusII file: " + m_filename );

    // ExodusII wants 32-bit integers as IDs of element ids
    std::vector< int > bface( s.second.size() );
    std::size_t i = 0;
    for (auto f : s.second) bface[ i++ ] = static_cast<int>(f)+1;
    // ExodusII wants 32-bit integers as element-relative face IDs
    const auto& fi = tk::cref_find( mesh.faceid(), s.first );
    std::vector< int > faceid( fi.size() );
    i = 0;
    for (auto f : fi) faceid[ i++ ] = static_cast<int>(f)+1;

    // Write side set data: ExodusII-file internal element ids adjacent to side
    // set and face id relative to element indicating which face is aligned with
    // the side set.
    ErrChk( ex_put_set( m_outFile, EX_SIDE_SET, s.first, bface.data(),
                        faceid.data() ) == 0,
      "Failed to write side set face list to ExodusII file: " + m_filename );
  }
}

void
ExodusIIMeshWriter::writeNodesets( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write side sets and their node list to ExodusII file
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  // Write all side set node lists in mesh
  for (const auto& s : mesh.bnode()) {
    // Write side set parameters
    ErrChk( ex_put_set_param( m_outFile, EX_NODE_SET, s.first,
                              static_cast<int64_t>(s.second.size()), 0 ) == 0,
      "Failed to write side set parameters to ExodusII file: " + m_filename );

    // ExodusII wants 32-bit integers as IDs of node ids
    std::vector< int > bnode( s.second.size() );
    std::size_t i = 0;
    for (auto n : s.second) bnode[ i++ ] = static_cast<int>(n)+1;

    // Write side set data
    ErrChk( ex_put_set( m_outFile, EX_NODE_SET, s.first, bnode.data(),
                        nullptr ) == 0,
      "Failed to write side set node list to ExodusII file: " + m_filename );
  }
}

void
ExodusIIMeshWriter::writeTimeStamp( uint64_t it, tk::real time ) const
// *****************************************************************************
//  Write time stamp to ExodusII file
//! \param[in] it Iteration number
//! \param[in] time Time
// *****************************************************************************
{
  ErrChk( ex_put_time( m_outFile, static_cast<int>(it), &time ) == 0,
          "Failed to write time stamp to ExodusII file: " + m_filename );
}

void
ExodusIIMeshWriter::writeTimeValues( const std::vector< tk::real >& tv ) const
// *****************************************************************************
//  Write time values to ExodusII file
//! \param[in] tv Time values for all time steps
// *****************************************************************************
{
   int i = 0;
   for (const auto& v : tv) {
     ErrChk( ex_put_time( m_outFile, ++i, &v ) == 0,
             "Failed to write time value for a time step to ExodusII file: " +
             m_filename );
   }
}

void
ExodusIIMeshWriter::writeNodeVarNames( const std::vector< std::string >& nv )
const
// *****************************************************************************
//  Write the names of nodal output variables to ExodusII file
//! \param[in] nv Nodal variable names
// *****************************************************************************
{
  if (!nv.empty()) {

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
}

void
ExodusIIMeshWriter::writeElemVarNames( const std::vector< std::string >& ev )
const
// *****************************************************************************
//  Write the names of element output variables to ExodusII file
//! \param[in] ev Elem variable names
// *****************************************************************************
{
  if (!ev.empty()) {

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
}

void
ExodusIIMeshWriter::writeNodeScalars(
  const std::vector< std::vector< std::vector< tk::real > > >& var ) const
// *****************************************************************************
//  Write multiple node scalar fields to ExodusII file at multiple time steps
//! \param[in] var Vector of nodal variables to read to: inner vector: nodes,
//!   middle vector: (physics) variable, outer vector: time step
// *****************************************************************************
{
  uint64_t time = 0;
  int varid = 0;

  for (const auto& t : var) {    // for all times
    ++time;
    for (const auto& v : t) {    // for all variables
      writeNodeScalar( time, ++varid, v );
    }
    varid = 0;
  }
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
// *****************************************************************************
{
  if (!var.empty()) {
    ErrChk( ex_put_var( m_outFile,
                        static_cast< int >( it ),
                        EX_NODE_BLOCK,
                        varid,
                        1,
                        static_cast< int64_t >( var.size() ),
                        var.data() ) == 0,
            "Failed to write node scalar to ExodusII file: " + m_filename );
  }
}

void
ExodusIIMeshWriter::writeElemScalars(
  const std::vector< std::vector< std::vector< tk::real > > >& var ) const
// *****************************************************************************
//  Write multiple element scalar fields to ExodusII file at multiple time steps
//! \param[in] var Vector of elemental variables to read to: inner vector:
//!   elements, middle vector: (physics) variable, outer vector: time step
// *****************************************************************************
{
  uint64_t time = 0;
  int varid = 0;

  for (const auto& t : var) {    // for all times
    ++time;
    for (const auto& v : t) {    // for all variables
      writeElemScalar( time, ++varid, v );
    }
    varid = 0;
  }
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
// *****************************************************************************
{
  if (!var.empty()) {
    ErrChk( ex_put_var( m_outFile,
                        static_cast< int >( it ),
                        EX_ELEM_BLOCK,
                        varid,
                        1,
                        static_cast< int64_t >( var.size() ),
                        var.data() ) == 0,
            "Failed to write elem scalar to ExodusII file: " + m_filename );
  }
}
