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
#include <util/NamedValue.hpp>

using namespace phdmesh ;

//----------------------------------------------------------------------

void test_comm_all( ParallelMachine comm , std::istream & s )
{
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  //--------------------------------------------------------------------
  // Test parameters

  NamedValue<unsigned> num_cycle("cycle");
  NamedValue<unsigned> dense_flag("dense");
  NamedValue<unsigned> max_msg_length("length");
  NamedValue< std::vector<int> > neighbor_template("neighbors");

  num_cycle.value = 10 ;
  dense_flag.value = 0 ;
  max_msg_length.value = 1000 ;

  NamedValueSet input_values ;
  input_values.insert( & dense_flag );
  input_values.insert( & max_msg_length );
  input_values.insert( & neighbor_template );

  s >> input_values ;

  std::vector<unsigned> neighbors ;

  if ( p_rank == 0 ) { std::cout << input_values << std::endl ; }

  for ( std::vector<int>::iterator
        i =  neighbor_template.value.begin() ;
        i != neighbor_template.value.end() ; ++i ) {
    int p = ( *i + p_rank ) % p_size ;
    if ( p < 0 ) { p += p_size ; }
    unsigned n = (unsigned) p ;
    neighbors.push_back(n);
  }

  //--------------------------------------------------------------------

  const unsigned zero = 0 ;

  std::vector<unsigned> send_size( p_size , zero );
  std::vector<unsigned> recv_size( p_size , zero );

  for ( unsigned i = 0 ; i < neighbors.size() ; ++i ) {
    send_size[ neighbors[i] ] = max_msg_length.value * sizeof(unsigned);
  }

  //--------------------------------------------------------------------

  unsigned num_msg_maximum ;
  int flag = 0 ;
  flag = comm_sizes( comm , p_size , num_msg_maximum ,
                       & send_size[0] , & recv_size[0] , flag );

  const unsigned p_msg_max = dense_flag.value ? 0 : p_size ;
  CommAll cs ;
  flag = cs.allocate_buffers( comm , p_msg_max ,
                              & send_size[0] , & recv_size[0] , flag );

  for ( unsigned i = 0 ; i < neighbors.size() ; ++i ) {
    CommBuffer & b = cs.send_buffer( neighbors[i] );
    for ( unsigned k = 0 ; k < max_msg_length.value ; ++k ) {
      b.pack<unsigned>( p_rank );
    }
  }

  double t = wall_time();

  for ( unsigned i = 0 ; i < num_cycle.value ; ++i ) {
    cs.communicate();
  }
  double dt_max = wall_dtime( t );
  double dt_min = dt_max ;

  all_reduce( comm , Min<1>( & dt_min ) , Max<1>( & dt_max ) );

  dt_max /= num_cycle.value ;
  dt_min /= num_cycle.value ;

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & b = cs.recv_buffer(p);

    while ( b.remaining() ) {
      unsigned tmp ; b.unpack<unsigned>( tmp );
      if ( tmp != p ) {
        throw std::runtime_error( std::string("FAIL comm_all CONTENT") );
      }
    }
  }

  if ( p_rank == 0 ) {
    std::cout << "comm_all PASS ( min , max ) = ( "
              << dt_min << " , " << dt_max << " )"
              << std::endl << std::endl ;
    std::cout.flush();
  }
}

