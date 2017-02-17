/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef TESTFEMESHBOXFIXTURE_HPP
#define TESTFEMESHBOXFIXTURE_HPP

#include <stdio.h>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>
#include <BoxMeshFixture.hpp>

#include <ParallelComm.hpp>

//----------------------------------------------------------------------------

namespace TestFEMesh {

template< class FEMeshType >
struct VerifyUnpack  ;

}

//----------------------------------------------------------------------------

#ifdef HAVE_MPI

namespace TestFEMesh {

template< typename coordinate_scalar_type ,
          unsigned ElemNodeCount ,
          class Device >
void verify_parallel(
  const HybridFEM::FEMesh< coordinate_scalar_type ,
                           ElemNodeCount ,
                           Device > & mesh )
{
  typedef HybridFEM::FEMesh< coordinate_scalar_type, ElemNodeCount, Device > femesh_type ;
  typedef typename femesh_type::node_coords_type node_coords_type ;

  comm::Machine machine = mesh.parallel_data_map.machine ;

  // Communicate node coordinates to verify communication and setup.

  const size_t chunk_size = 3 ;

  KokkosArray::AsyncExchange< coordinate_scalar_type, Device, KokkosArray::ParallelDataMap >
    exchange( mesh.parallel_data_map , chunk_size );

  const size_t send_begin = mesh.parallel_data_map.count_interior ;
  const size_t send_count = mesh.parallel_data_map.count_send ;

  const size_t recv_begin = mesh.parallel_data_map.count_owned ;
  const size_t recv_count = mesh.parallel_data_map.count_receive ;

  typedef KokkosArray::PackArray< node_coords_type > pack_type ;

  pack_type::pack( exchange.buffer(), send_begin, send_count, mesh.node_coords );

  exchange.setup();

  // Launch local-action device kernels

  exchange.send_receive();

  unsigned long local[3] ;
  local[0] = mesh.parallel_data_map.count_owned ;
  local[1] = mesh.parallel_data_map.count_receive ;
  local[2] = TestFEMesh::VerifyUnpack< node_coords_type >::unpack( mesh.node_coords, recv_begin, recv_count, exchange.buffer() );

  unsigned long global[3] = { 0 , 0 , 0 };

  MPI_Allreduce( local , global ,
                 3 , MPI_UNSIGNED_LONG , MPI_SUM , machine.mpi_comm );

  if ( 0 == comm::rank( machine ) ) {
    std::cout << ( global[2] ? "FAILED" : "PASSED" )
              << ": TestFEMesh::verify_parallel "
              << "NP(" << comm::size( machine )
              << ") total_node(" << global[0]
              << ") verified_nodes(" << global[1]
              << ") failed_nodes(" << global[2]
              << ")" << std::endl ;
  }
}

} // namespace TestFEMesh

#else /* ! #ifdef HAVE_MPI */

namespace TestFEMesh {

template< typename coordinate_scalar_type ,
          unsigned ElemNodeCount ,
          class Device >
void verify_parallel(
  const HybridFEM::FEMesh< coordinate_scalar_type ,
                           ElemNodeCount ,
                           Device > & mesh )
{}

} // namespace TestFEMesh

#endif /* ! #ifdef HAVE_MPI */

//----------------------------------------------------------------------------

template< class Device >
void test_box_fixture( comm::Machine machine ,
                       const size_t nodes_nx ,
                       const size_t nodes_ny ,
                       const size_t nodes_nz )
{
  typedef int coordinate_scalar_type ;

  const size_t proc_count = comm::size( machine );
  const size_t proc_local = comm::rank( machine ) ;

  HybridFEM::FEMesh< coordinate_scalar_type , 8 , Device > mesh =
    box_mesh_fixture< coordinate_scalar_type , Device >
      ( proc_count , proc_local , nodes_nx , nodes_ny , nodes_nz );

  mesh.parallel_data_map.machine = machine ;

  TestFEMesh::verify_parallel( mesh );
}

#endif /* #ifndef TESTFEMESHBOXFIXTURE_HPP */


