/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @file
 * @author H. Carter Edwards
 * @date   January 2007
 */

#include <iostream>
#include <vector>
#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>

using namespace phdmesh ;

// Stress test for communication: message sizes

void test_comm_bounds( ParallelMachine comm , std::istream & )
{
  static const char method[] = "phdmesh::test_comm_sizes" ;

  const unsigned u_zero = 0 ;

  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  std::vector<unsigned> send_size( p_size , u_zero );
  std::vector<unsigned> recv_size( p_size , u_zero );

  const unsigned * const ptr_send_size = & send_size[0] ;
  const unsigned * const ptr_recv_size = & recv_size[0] ;

  if ( 1 < p_size ) {

    // Send one point-to-point message doubling the size until it breaks
    // Send one message to next processor,
    // Receive one message from previous processor:

    const unsigned p_send = ( p_rank + 1 )          % p_size ;
    const unsigned p_recv = ( p_rank + p_size - 1 ) % p_size ;

    int error = 0 ;

    const unsigned msg_max = 0x8000000 ; // max to test 100 Mb
    unsigned msg_size = 0x1000 ; // 4k to start

    while ( ! error && msg_size <= msg_max ) {

      send_size[ p_send ] = msg_size ;
      recv_size[ p_recv ] = msg_size ;

      try {
        CommAll all ;

        all.allocate_buffers( comm , p_size , ptr_send_size , ptr_recv_size );

        // Fill send buffer with predefined values

        CommBuffer & send_buf = all.send_buffer( p_send );

        const unsigned n = msg_size / sizeof(unsigned);
        for ( unsigned i = 0 ; i < n ; ++i ) {
          send_buf.pack<unsigned>( msg_size );
        }

        all.communicate();

        // unpack the receive buffer, verify contents

        CommBuffer & recv_buf = all.recv_buffer( p_recv );

        for ( unsigned i = 0 ; i < n ; ++i ) {
          unsigned tmp ;
          recv_buf.unpack<unsigned>( tmp );
          if ( tmp != msg_size ) { error = 1 ; }
        }
      }
      catch ( const std::exception & x ) {
        if ( p_rank == 0 ) {
          std::cout << method << x.what() << std::endl ;
        }
        error = 2 ;
      }
      catch ( ... ) {
        error = 2 ;
      }

      all_reduce( comm , Max<1>( & error ) );

      if ( ! error ) {
        if ( p_rank == 0 ) {
          std::cout << "OK Message size = "
                    << msg_size << std::endl ;
        }
        msg_size <<= 1 ;
      }
    }

    if ( p_rank == 0 ) {
      if ( error == 1 ) {
        std::cout << method << " BAD message at "
                  << msg_size << " bytes" << std::endl ;
      }
    }
  }
}

