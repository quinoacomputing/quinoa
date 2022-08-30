// *****************************************************************************
/*!
  \file      src/IO/Omega_h_MeshReader.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Omega_h mesh reader
  \details   Omega_h mesh reader class definition.
*/
// *****************************************************************************

#include "NoWarning/Omega_h_file.hpp"
#include <Omega_h_library.hpp>

#include "Macro.hpp"
#include "Omega_h_MeshReader.hpp"
#include "Reorder.hpp"

using tk::Omega_h_MeshReader;

void
Omega_h_MeshReader::readMeshPart(
  std::vector< std::size_t >& ginpoel,
  std::vector< std::size_t >& inpoel,
  [[maybe_unused]] std::vector< std::size_t >& triinp,
  std::unordered_map< std::size_t, std::size_t >& lid,
  tk::UnsMesh::Coords& coord,
  std::unordered_map< std::size_t, std::set< std::size_t > >&,
  int numpes,
  [[maybe_unused]] int mype )
// *****************************************************************************
//  Read a part of the mesh (graph and coordinates) from Omega_h file
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
//! \param[in] numpes Total number of PEs (default n = 1, for a single-CPU read)
//! \param[in] mype This PE (default m = 0, for a single-CPU read)
//! \note The last two integer arguments are unused. They are needed because
//!   this function can be used via a polymorphic interface via a base class,
//!   see tk::MeshReader, and other specialized mesh readers, e.g.,
//!   tk::ExodusIIMeshReader, use these arguments. Here we require the Omega_h
//!   input mesh to be pre-partitioned, with Omega_h's osh_part tool, into the
//!   number of partitions that is equal to the number of PEs this fuinction is
//!   called on in parallel.
// *****************************************************************************
{
  Assert( mype < numpes, "Invalid input: PE id must be lower than NumPEs" );
  Assert( ginpoel.empty() && inpoel.empty() && lid.empty() &&
          coord[0].empty() && coord[1].empty() && coord[2].empty(),
          "Containers to store mesh must be empty" );

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wold-style-cast"
  #endif

  // Create Omega_h library instance
  auto lib = Omega_h::Library( nullptr, nullptr, MPI_COMM_WORLD );

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #endif

  // Find out how many partitions the Omega_h mesh was saved with
  auto nparts = Omega_h::binary::read_nparts( m_filename, lib.world() );

  if (numpes < nparts)
    Throw( "The Omega_h mesh reader only supports NumPEs >= nparts, where "
           "nparts is the number of partitions the mesh is partitioned into. "
           "Also note that NumPEs must be a power of 2 if NumPEs > nparts." );

  // Read mesh
  auto mesh = Omega_h::binary::read( m_filename, &lib );

  // Lambda to check if int is a power of two
  auto isPowerOfTwo = []( int x ) { return (x != 0) && ((x & (x - 1)) == 0); };

  if (nparts != numpes) {
    if (!isPowerOfTwo(numpes))
      Throw( "The Omega_h mesh reader only supports NumPEs of power of 2" );
    else
      mesh.balance();
  }

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
  std::vector< std::size_t > gid;
  std::tie( inpoel, gid, lid ) = tk::global2local( ginpoel );
}


std::vector< std::size_t >
Omega_h_MeshReader::triinpoel(
  [[maybe_unused]] std::map< int, std::vector< std::size_t > >& bface,
  [[maybe_unused]] const std::map< int, std::vector< std::size_t > >& faces,
  [[maybe_unused]] const std::vector< std::size_t >& ginpoel,
  [[maybe_unused]] const std::vector< std::size_t >& triinp ) const
// *****************************************************************************
// ...
//! \note Must be preceded by a call to readElemBlockIDs()
// *****************************************************************************
{
  std::vector< std::size_t > bnd_triinpoel;
  return bnd_triinpoel;
}

void
Omega_h_MeshReader::readSidesetFaces(
  [[maybe_unused]] std::map< int, std::vector< std::size_t > >& bface,
  [[maybe_unused]] std::map< int, std::vector< std::size_t > >& faces )
// *****************************************************************************
//  Read side sets from Omega_h file
//! \param[in,out] bface Elem ids of side sets to read into
//! \param[in,out] faces Elem-relative face ids of tets of side sets
// *****************************************************************************
{
}

void
Omega_h_MeshReader::readFaces(
  [[maybe_unused]] std::vector< std::size_t >& conn )
// *****************************************************************************
//  Read face connectivity of a number of boundary faces from Omega_h file
//! \param[in,out] conn Connectivity vector to push to
//! \details This function reads in the total number of boundary faces,
//!   also called triangle-elements, and their connectivity.
// *****************************************************************************
{
}

std::map< int, std::vector< std::size_t > >
Omega_h_MeshReader::readSidesetNodes()
// *****************************************************************************
//  Read node list of all side sets from Omega_h file
//! \return Node lists mapped to side set ids
// *****************************************************************************
{
  // Node lists mapped to side set ids
  std::map< int, std::vector< std::size_t > > side;

  return side;
}
