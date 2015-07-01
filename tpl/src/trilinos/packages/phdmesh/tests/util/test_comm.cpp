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
 * @author H. Carter Edwards
 */

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>

using namespace phdmesh ;

//----------------------------------------------------------------------

void test_comm_sparse( ParallelMachine comm , std::istream & )
{
  const unsigned zero = 0 ;
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  std::vector<unsigned> send_size( p_size , zero );
  std::vector<unsigned> recv_size( p_size , zero );
  int flag ;

  unsigned num_msg_maximum = 0 ;

  //--------------------------------------------------------------------
  // Test sizing zero send and zero flag

  flag = 0 ;
  flag = comm_sizes( comm , p_size , num_msg_maximum ,
                     & send_size[0] , & recv_size[0] , flag );

  for ( unsigned i = 0 ; i < p_size ; ++i ) {
    if ( recv_size[i] != 0 || flag != 0 || num_msg_maximum ) {
      throw std::runtime_error( std::string( "FAIL SPARSE TEST 1" ) );
    }
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS SPARSE TEST 1 for NP = " << p_size << std::endl ;
    std::cout.flush();
  }
  //--------------------------------------------------------------------
  // Test sizing zero send and non zero flag on one processor

  flag = p_rank == 0  ;
  flag = comm_sizes( comm , p_size , num_msg_maximum ,
                     & send_size[0] , & recv_size[0] , flag );

  for ( unsigned i = 0 ; i < p_size ; ++i ) {
    if ( recv_size[i] != 0 || ! flag || num_msg_maximum ) {
      throw std::runtime_error( std::string( "FAIL SPARSE TEST 2" ) );
    }
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS SPARSE TEST 2 for NP = " << p_size << std::endl ;
    std::cout.flush();
  }
  //--------------------------------------------------------------------
  // Test sizing send one integer to previous and next processor
  {
    const unsigned p_down = ( p_rank + p_size - 1 ) % p_size ;
    const unsigned p_up   = ( p_rank + 1 )      % p_size ;
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      unsigned n = sizeof(int) * ( ( i == p_down || i == p_up ) ? 1 : 0 );
      send_size[i] = n ;
    }

    flag = 0 ;
    flag = comm_sizes( comm , p_size , num_msg_maximum ,
                       & send_size[0] , & recv_size[0] , flag );

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      unsigned n = sizeof(int) * ( ( i == p_down || i == p_up ) ? 1 : 0 );
      if ( recv_size[i] != n || flag || 2 < num_msg_maximum ) {
        throw std::runtime_error( std::string( "FAIL SPARSE TEST 3" ) );
      }
    }
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS SPARSE TEST 3 for NP = " << p_size << std::endl ;
    std::cout.flush();
  }
  //--------------------------------------------------------------------
  // Test sending one integer to previous and next processor
  {
    const unsigned p_down = ( p_rank + p_size - 1 ) % p_size ;
    const unsigned p_up   = ( p_rank + 1 )      % p_size ;
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      unsigned n = sizeof(int) * ( ( i == p_down || i == p_up ) ? 1 : 0 );
      send_size[i] = n ;
    }

    flag = 0 ;
    flag = comm_sizes( comm , p_size , num_msg_maximum ,
                       & send_size[0] , & recv_size[0] , flag );
    if ( 2 < num_msg_maximum ) {
      throw std::runtime_error( std::string( "FAIL SPARSE TEST 4A" ) );
    }

    CommAll cs ;
    flag = cs.allocate_buffers( comm , p_size ,
                                & send_size[0] , & recv_size[0] , flag );

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      const bool send = i == p_down || i == p_up ;
      if ( send ) { cs.send_buffer(i).pack<unsigned>(p_rank); }
    }

    cs.communicate();

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      const bool recv = i == p_down || i == p_up ;
      if ( recv ) {
        if ( cs.recv_buffer(i).remaining() != (int) sizeof(unsigned) ) {
          throw std::runtime_error( std::string( "FAIL SPARSE TEST 4B" ) );
        }
        unsigned tmp ;
        cs.recv_buffer(i).unpack<unsigned>(tmp);
        if ( tmp != i || flag ) {
          throw std::runtime_error( std::string( "FAIL SPARSE TEST 4C" ) );
        }
      }
    }
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      if ( cs.recv_buffer(i).remaining() ) {
        throw std::runtime_error( std::string( "FAIL SPARSE TEST 4C" ) );
      }
    }
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS SPARSE TEST 4 for NP = " << p_size << std::endl ;
    std::cout.flush();
  }
  //--------------------------------------------------------------------
  // Test sending 'my_rank + 1' integers to previous and next processor
  {
    const unsigned p_down = ( p_rank + p_size - 1 ) % p_size ;
    const unsigned p_up   = ( p_rank + 1 )      % p_size ;

    CommAll cs( comm );

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      unsigned n = ( i == p_down || i == p_up ) ? ( p_rank + 1 ) : 0 ;
      for ( unsigned j = 0 ; j < n ; ++j ) {
        cs.send_buffer(i).pack<unsigned>( j );
      }
    }

    flag = 0 ;
    flag = cs.allocate_buffers( p_size , false , flag );

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      unsigned n = ( i == p_down || i == p_up ) ? ( p_rank + 1 ) : 0 ;
      for ( unsigned j = 0 ; j < n ; ++j ) {
        cs.send_buffer(i).pack<unsigned>( j );
      }
    }

    cs.communicate();

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      unsigned n = ( i == p_down || i == p_up ) ? ( i + 1 ) : 0 ;
      if ( cs.recv_buffer(i).remaining() != (int)( n * sizeof(unsigned) ) ) {
        throw std::runtime_error( std::string( "FAIL SPARSE TEST 5A" ) );
      }
      for ( unsigned j = 0 ; j < n ; ++j ) {
        unsigned tmp ;
        cs.recv_buffer(i).unpack<unsigned>( tmp );
        if ( tmp != j ) {
          throw std::runtime_error( std::string( "FAIL SPARSE TEST 5B" ) );
        }
      }
    }
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      if ( cs.recv_buffer(i).remaining() ) {
        throw std::runtime_error( std::string( "FAIL SPARSE TEST 5C" ) );
      }
    }
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS SPARSE TEST 5 for NP = " << p_size << std::endl ;
    std::cout.flush();
  }
  //--------------------------------------------------------------------
  // Test sending 'my_rank + 1' integers to '1 + my_rank % 3' other processors
  {
    CommAll cs( comm );

    unsigned n_send = 1 + ( p_rank % 3 );
    unsigned n_size = p_rank + 1 ;

    for ( unsigned i = 0 ; i < n_send ; ++i ) {
      unsigned dst = ( p_rank + 1 + i ) % p_size ;
      for ( unsigned j = 0 ; j < n_size ; ++j ) {
        cs.send_buffer(dst).pack<unsigned>( j );
      }
    }

    flag = 0 ;
    flag = cs.allocate_buffers( p_size , false , flag );

    for ( unsigned i = 0 ; i < n_send ; ++i ) {
      unsigned dst = ( p_rank + 1 + i ) % p_size ;
      for ( unsigned j = 0 ; j < n_size ; ++j ) {
        cs.send_buffer(dst).pack<unsigned>( j );
      }
    }

    cs.communicate();

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      // How many messages is 'i' send?
      unsigned i_n_send = 1 + ( i % 3 );
      unsigned i_n_size = i + 1 ;
      // Should 'i' have sent to me?
      for ( unsigned j = 0 ; j < i_n_send ; ++j ) {
        unsigned i_dst = ( i + 1 + j ) % p_size ;
        if ( p_rank == i_dst ) {
          if ( cs.recv_buffer(i).remaining() !=
               (int)( i_n_size * sizeof(unsigned) ) ) {
            throw std::runtime_error( std::string( "FAIL SPARSE TEST 6A" ) );
          }
          for ( unsigned k = 0 ; k < i_n_size ; ++k ) {
            unsigned tmp ;
            cs.recv_buffer(i).unpack<unsigned>( tmp );
            if ( tmp != k ) {
              throw std::runtime_error( std::string( "FAIL SPARSE TEST 6B" ) );
            }
          }
        }
      }
    }
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      if ( cs.recv_buffer(i).remaining() ) {
        throw std::runtime_error( std::string( "FAIL SPARSE TEST 6C" ) );
      }
    }
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS SPARSE TEST 6 for NP = " << p_size << std::endl ;
    std::cout.flush();
  }
  //--------------------------------------------------------------------
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void test_comm_dense( ParallelMachine comm , std::istream & )
{
  const unsigned zero = 0 ;
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  std::vector<unsigned> send_size( p_size , zero );
  std::vector<unsigned> recv_size( p_size , zero );
  int flag ;

  //--------------------------------------------------------------------
  // Test sizing zero send and zero flag

  flag = 0 ;
  flag = comm_dense_sizes( comm , & send_size[0] , & recv_size[0] , flag );

  for ( unsigned i = 0 ; i < p_size ; ++i ) {
    if ( recv_size[i] != 0 || flag != 0 ) {
      throw std::runtime_error( std::string( "FAIL DENSE TEST 1" ) );
    }
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS DENSE TEST 1 for NP = " << p_size << std::endl ;
    std::cout.flush();
  }
  //--------------------------------------------------------------------
  // Test sizing zero send and non zero flag on one processor

  flag = p_rank == 0  ;
  flag = comm_dense_sizes( comm , & send_size[0] , & recv_size[0] , flag );

  for ( unsigned i = 0 ; i < p_size ; ++i ) {
    if ( recv_size[i] != 0 || ! flag ) {
      throw std::runtime_error( std::string( "FAIL DENSE TEST 2" ) );
    }
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS DENSE TEST 2 for NP = " << p_size << std::endl ;
    std::cout.flush();
  }
  //--------------------------------------------------------------------
  // Test sizing send one integer to previous and next processor
  {
    const unsigned p_down = ( p_rank + p_size - 1 ) % p_size ;
    const unsigned p_up   = ( p_rank + 1 )      % p_size ;
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      unsigned n = sizeof(int) * ( ( i == p_down || i == p_up ) ? 1 : 0 );
      send_size[i] = n ;
    }

    flag = 0 ;
    flag = comm_dense_sizes( comm , & send_size[0] , & recv_size[0] , flag );

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      unsigned n = sizeof(int) * ( ( i == p_down || i == p_up ) ? 1 : 0 );
      if ( recv_size[i] != n || flag ) {
        throw std::runtime_error( std::string( "FAIL DENSE TEST 3" ) );
      }
    }
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS DENSE TEST 3 for NP = " << p_size << std::endl ;
    std::cout.flush();
  }
  //--------------------------------------------------------------------
  // Test sending one integer to previous and next processor
  {
    const unsigned p_down = ( p_rank + p_size - 1 ) % p_size ;
    const unsigned p_up   = ( p_rank + 1 )      % p_size ;
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      unsigned n = sizeof(int) * ( ( i == p_down || i == p_up ) ? 1 : 0 );
      send_size[i] = n ;
    }

    flag = 0 ;
    CommAll cs ;
    flag = cs.allocate_buffers( comm , 0 ,
                                & send_size[0] , & recv_size[0] , flag );

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      const bool send = i == p_down || i == p_up ;
      if ( send ) { cs.send_buffer(i).pack<unsigned>(p_rank); }
    }

    cs.communicate();

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      const bool recv = i == p_down || i == p_up ;
      if ( recv ) {
        if ( cs.recv_buffer(i).remaining() != (int) sizeof(unsigned) ) {
          throw std::runtime_error( std::string( "FAIL DENSE TEST 4A" ) );
        }
        unsigned tmp ;
        cs.recv_buffer(i).unpack<unsigned>(tmp);
        if ( tmp != i || flag ) {
          throw std::runtime_error( std::string( "FAIL DENSE TEST 4B" ) );
        }
      }
    }
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      if ( cs.recv_buffer(i).remaining() ) {
        throw std::runtime_error( std::string( "FAIL DENSE TEST 4C" ) );
      }
    }
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS DENSE TEST 4 for NP = " << p_size << std::endl ;
    std::cout.flush();
  }
  //--------------------------------------------------------------------
  // Test sending 'my_rank + 1' integers to previous and next processor
  {
    const unsigned p_down = ( p_rank + p_size - 1 ) % p_size ;
    const unsigned p_up   = ( p_rank + 1 )      % p_size ;

    CommAll cs( comm );

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      unsigned n = ( i == p_down || i == p_up ) ? ( p_rank + 1 ) : 0 ;
      for ( unsigned j = 0 ; j < n ; ++j ) {
        cs.send_buffer(i).pack<unsigned>( j );
      }
    }

    flag = 0 ;
    flag = cs.allocate_buffers( 0 , false , flag );

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      unsigned n = ( i == p_down || i == p_up ) ? ( p_rank + 1 ) : 0 ;
      for ( unsigned j = 0 ; j < n ; ++j ) {
        cs.send_buffer(i).pack<unsigned>( j );
      }
    }

    cs.communicate();

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      unsigned n = ( i == p_down || i == p_up ) ? ( i + 1 ) : 0 ;
      if ( cs.recv_buffer(i).remaining() != (int)( n * sizeof(unsigned) ) ) {
        throw std::runtime_error( std::string( "FAIL DENSE TEST 5A" ) );
      }
      for ( unsigned j = 0 ; j < n ; ++j ) {
        unsigned tmp ;
        cs.recv_buffer(i).unpack<unsigned>( tmp );
        if ( tmp != j ) {
          throw std::runtime_error( std::string( "FAIL DENSE TEST 5B" ) );
        }
      }
    }
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      if ( cs.recv_buffer(i).remaining() ) {
        throw std::runtime_error( std::string( "FAIL DENSE TEST 5C" ) );
      }
    }
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS DENSE TEST 5 for NP = " << p_size << std::endl ;
    std::cout.flush();
  }
  //--------------------------------------------------------------------
  // Test sending 'my_rank + 1' integers to '1 + my_rank % 3' other processors
  {
    CommAll cs( comm );

    unsigned n_send = 1 + ( p_rank % 3 );
    unsigned n_size = p_rank + 1 ;

    for ( unsigned i = 0 ; i < n_send ; ++i ) {
      unsigned dst = ( p_rank + 1 + i ) % p_size ;
      for ( unsigned j = 0 ; j < n_size ; ++j ) {
        cs.send_buffer(dst).pack<unsigned>( j );
      }
    }

    flag = 0 ;
    flag = cs.allocate_buffers( 0 , false , flag );

    for ( unsigned i = 0 ; i < n_send ; ++i ) {
      unsigned dst = ( p_rank + 1 + i ) % p_size ;
      for ( unsigned j = 0 ; j < n_size ; ++j ) {
        cs.send_buffer(dst).pack<unsigned>( j );
      }
    }

    cs.communicate();

    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      // How many messages is 'i' send?
      unsigned i_n_send = 1 + ( i % 3 );
      unsigned i_n_size = i + 1 ;
      // Should 'i' have sent to me?
      for ( unsigned j = 0 ; j < i_n_send ; ++j ) {
        unsigned i_dst = ( i + 1 + j ) % p_size ;
        if ( p_rank == i_dst ) {
          if ( cs.recv_buffer(i).remaining() !=
               (int)( i_n_size * sizeof(unsigned) ) ) {
            throw std::runtime_error( std::string( "FAIL DENSE TEST 6A" ) );
          }
          for ( unsigned k = 0 ; k < i_n_size ; ++k ) {
            unsigned tmp ;
            cs.recv_buffer(i).unpack<unsigned>( tmp );
            if ( tmp != k ) {
              throw std::runtime_error( std::string( "FAIL DENSE TEST 6B" ) );
            }
          }
        }
      }
    }
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      if ( cs.recv_buffer(i).remaining() ) {
        throw std::runtime_error( std::string( "FAIL DENSE TEST 6C" ) );
      }
    }
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS DENSE TEST 6 for NP = " << p_size << std::endl ;
    std::cout.flush();
  }
  //--------------------------------------------------------------------
}

