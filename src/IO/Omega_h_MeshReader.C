// *****************************************************************************
/*!
  \file      src/IO/Omega_h_MeshReader.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Omega_h mesh reader
  \details   Omega_h mesh reader class definition.
*/
// *****************************************************************************

#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>

#include "Macro.h"
#include "Omega_h_MeshReader.h"
#include "Reorder.h"

using tk::Omega_h_MeshReader;

void
Omega_h_MeshReader::readMeshPart(
  std::vector< std::size_t >& ginpoel,
  std::vector< std::size_t >& inpoel,
  std::vector< std::size_t >& gid,
  std::unordered_map< std::size_t, std::size_t >& lid,
  tk::UnsMesh::Coords& coord,
  int n, int m )
// *****************************************************************************
//  Read a part of the mesh (graph and coordinates) from Omega_h file
//! \param[in,out] ginpoel Container to store element connectivity of this PE's
//!   chunk of the mesh (global ids)
//! \param[in,out] inpoel Container to store element connectivity with local
//!   node IDs of this PE's mesh chunk
//! \param[in,out] gid Container to store global node IDs of elements of this
//!   PE's mesh chunk
//! \param[in,out] lid Container to store global->local node IDs of elements of
//!   this PE's mesh chunk
//! \param[in,out] coord Container to store coordinates of mesh nodes of this
//!   PE's mesh chunk
//! \param[in] n Total number of PEs (default n = 1, for a single-CPU read)
//! \param[in] m This PE (default m = 0, for a single-CPU read)
//! \note The last two integer arguments are unused. They are needed because
//!   this function can be used via a polymorphic interface via a base class,
//!   see tk::MeshReader, and other specialized mesh readers, e.g.,
//!   tk::ExodusIIMeshReader, use these arguments. Here we require the Omega_h
//!   input mesh to be pre-partitioned, with Omega_h's osh_part tool, into the
//!   number of partitions that is equal to the number of PEs this fuinction is
//!   called on in parallel.
// *****************************************************************************
{
  IGNORE( n );  // Avoid compiler warning on unused arguments in RELEASE mode
  IGNORE( m );
  Assert( m < n, "Invalid input: PE id must be lower than NumPEs" );
  Assert( ginpoel.empty() && inpoel.empty() && gid.empty() && lid.empty() &&
          coord[0].empty() && coord[1].empty() && coord[2].empty(),
          "Containers to store mesh must be empty" );

  // Create Omega_h library instance
  auto lib = Omega_h::Library( nullptr, nullptr, MPI_COMM_WORLD );

  // Read mesh
  auto mesh = Omega_h::binary::read( m_filename, &lib );

  // Extract connectivity from Omega_h's mesh object
  auto ntets = mesh.nelems();
  ginpoel.resize( static_cast< std::size_t >( ntets ) * 4 );
  auto o_inpoel = mesh.ask_elem_verts();
  auto o_gid = mesh.globals( Omega_h::VERT );
  for (int i=0; i<ntets; ++i)
    for (int j=0; j<4; ++j) {
      auto I = static_cast< std::size_t >( i );
      auto J = static_cast< std::size_t >( j );
      ginpoel[ I*4+J ] = static_cast<std::size_t>( o_gid[ o_inpoel[i*4+j] ] );
    }

  // Extract node coordinates from Omega_h's mesh object
  auto o_coord = mesh.coords();

  // Extract number of vertices from Omega_h's mesh object
  auto nnode = static_cast< std::size_t >( mesh.nverts() );

  // Extract node coordinates from Omega_h's mesh object
  auto& x = coord[0];
  auto& y = coord[1];
  auto& z = coord[2];
  x.resize( nnode );
  y.resize( nnode );
  z.resize( nnode );

  for (std::size_t I=0; I<nnode; ++I) {
    auto i = static_cast< int >( I );
    x[I] = o_coord[ i*3+0 ];
    y[I] = o_coord[ i*3+1 ];
    z[I] = o_coord[ i*3+2 ];
  }

  // Compute local data from global mesh connectivity
  std::tie( inpoel, gid, lid ) = tk::global2local( ginpoel );
}

std::size_t
Omega_h_MeshReader::readSidesetFaces(
  std::map< int, std::vector< std::size_t > >& bface,
  std::map< int, std::vector< int > >& faceid )
// *****************************************************************************
//  Read face list of all side sets from Omega_h file
//! \param[in,out] bface Face-Element lists mapped to side set ids
//! \param[in,out] faceid Side set side lists associated to side set ids
//! \return Total number of boundary faces
// *****************************************************************************
{
IGNORE(bface);
IGNORE(faceid);
  return 0;
}

void
Omega_h_MeshReader::readFaces( std::size_t nbfac,
                               std::vector< std::size_t >& conn ) const
// *****************************************************************************
//  Read face connectivity of a number of boundary faces from Omega_h file
//! \param[in] nbfac Number of boundary faces
//! \param[inout] conn Connectivity vector to push to
//! \details This function reads in the total number of boundary faces,
//!   also called triangle-elements, and their connectivity.
//! \note It is okay to call this function with zero nbfac: it will be no-op.
// *****************************************************************************
{
  // Return if no boundary faces in file
  if (nbfac == 0) return;
IGNORE(conn);
}

std::map< int, std::vector< std::size_t > >
Omega_h_MeshReader::readSidesets()
// *****************************************************************************
//  Read node list of all side sets from Omega_h file
//! \return Node lists mapped to side set ids
// *****************************************************************************
{
  // Node lists mapped to side set ids
  std::map< int, std::vector< std::size_t > > side;

  return side;
}
