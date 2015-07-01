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

#include <limits>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <util/Parallel.hpp>
#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>
#include <util/OctTree.hpp>
#include <util/OctTreeOps.hpp>

using namespace phdmesh ;

void test_oct_tree_keys()
{
  std::cout << "OctTreeSize<0>::value = "
            << (unsigned) OctTreeSize<0>::value << std::endl ;
  std::cout << "OctTreeSize<1>::value = "
            << (unsigned) OctTreeSize<1>::value << std::endl ;
  std::cout << "OctTreeSize<2>::value = "
            << (unsigned) OctTreeSize<2>::value << std::endl ;
  std::cout << "OctTreeSize<3>::value = "
            << (unsigned) OctTreeSize<3>::value << std::endl ;
  std::cout << "OctTreeSize<4>::value = "
            << (unsigned) OctTreeSize<4>::value << std::endl ;
  std::cout << "OctTreeSize<5>::value = "
            << (unsigned) OctTreeSize<5>::value << std::endl ;
  std::cout << "OctTreeSize<6>::value = "
            << (unsigned) OctTreeSize<6>::value << std::endl ;
  std::cout << "OctTreeSize<7>::value = "
            << (unsigned) OctTreeSize<7>::value << std::endl ;
  std::cout << "OctTreeSize<8>::value = "
            << (unsigned) OctTreeSize<8>::value << std::endl ;
  std::cout << "OctTreeSize<9>::value = "
            << (unsigned) OctTreeSize<9>::value << std::endl ;
  std::cout << "OctTreeSize<10>::value = "
            << (unsigned) OctTreeSize<10>::value << std::endl ;

  std::cout << std::endl ;

  OctTreeKey key ;
  unsigned count = 0 ;

  enum { Depth = 4 };

  try {
    if ( oct_tree_offset( Depth , key ) != count ) throw count ;
    if ( key.depth() != 0 ) throw count ;
    ++count ;

    for ( unsigned i = 0 ; i < 8 ; ++i ) {
      key.set_index<1>( i + 1 );

      if ( oct_tree_offset( Depth , key ) != count ) throw count ;
      if ( key.depth() != 1 ) throw count ;
      ++count ;

      for ( unsigned j = 0 ; j < 8 ; ++j ) {
        key.set_index<2>( j + 1 );

        if ( oct_tree_offset( Depth , key ) != count ) throw count ;
        if ( key.depth() != 2 ) throw count ;
        ++count ;

        for ( unsigned k = 0 ; k < 8 ; ++k ) {
          key.set_index<3>( k + 1 );

          if ( oct_tree_offset( Depth , key ) != count ) throw count ;
          if ( key.depth() != 3 ) throw count ;
          ++count ;

          for ( unsigned l = 0 ; l < 8 ; ++l ) {
            key.set_index<4>( l + 1 );

            if ( oct_tree_offset( Depth , key ) != count ) throw count ;
            if ( key.depth() != 4 ) throw count ;
            ++count ;
          }
          key.clear_index<4>();
        }
        key.clear_index<3>();
      }
      key.clear_index<2>();
    }
    if ( count != OctTreeSize<Depth>::value ) throw count ;
  }
  catch( unsigned ix ) {
    std::cout << "oct_tree_offset" << key
              << " = " << oct_tree_offset( Depth , key )
              << " != " << ix << std::endl ;
  }

  try {
    OctTreeKey k ;
    for ( unsigned i = 1 ; i <= 16 ; ++i ) {
      k.set_index( i , 1 );
      if ( k.depth() != i ) throw i ;
    }
  }
  catch( unsigned e ) {
    std::cout << "k.depth() != " << e << std::endl ;
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

OctTreeKey to_key( const unsigned val )
{
  enum { bits = std::numeric_limits<unsigned>::digits };
  enum { mask = 0x07 };
  OctTreeKey k ;
  for ( unsigned d = 1 ; d <= 10 ; ++d ) {
    const unsigned shift = bits - d * 3 ;
    const unsigned index = ( val >> shift ) & mask ;
    k.set_index( d , index + 1 );
  }
  return k ;
}

unsigned to_uint( const OctTreeKey & k )
{
  enum { bits = std::numeric_limits<unsigned>::digits };
  enum { mask = 0x07 };
  unsigned val = 0 ;
  for ( unsigned d = 1 ; d <= 10 ; ++d ) {
    const unsigned index = k.index(d) ;
    if ( index ) {
      const unsigned shift = bits - d * 3 ;
      val |= ( index - 1 ) << shift ;
    }
  }
  return val ;
}

void proc_weights( ParallelMachine comm ,
                   const std::vector<OctTreeKey> & cut_keys ,
                   const std::vector<OctTreeKey> & keys ,
                   const std::vector<float>      & weights ,
                         std::vector<float>      & w_local ,
                         std::vector<float>      & w_global )
{
  const float fzero = 0.0 ;
  const unsigned ncuts = cut_keys.size();

  w_local.assign(  ncuts + 1 , fzero );
  w_global.assign( ncuts + 1 , fzero );

  const std::vector<OctTreeKey>::const_iterator ib = cut_keys.begin();
  const std::vector<OctTreeKey>::const_iterator ie = cut_keys.end();

  for ( unsigned i = 0 ; i < keys.size() ; ++i ) {
    const std::vector<OctTreeKey>::const_iterator ik =
      std::lower_bound( ib , ie , keys[i] );
    const unsigned p = ( ik - ib ) - 1 ;
    w_local[p] += weights[i] ;
    w_local[ ncuts ] += weights[i] ;
  }

  all_reduce_sum( comm , & w_local[0] , & w_global[0] , ncuts + 1 );
}

void print_cuts( const std::vector< OctTreeKey > & cut_keys ,
                 const std::vector< float >      & weights )
{
  const unsigned p_size = cut_keys.size();
  const double total = weights[ p_size ];
  const double mean  = total / ((double) p_size );

  for ( unsigned i = 0 ; i < p_size ; ++i ) {
    std::cout << "  P" << i << " cut begin = " ;
    std::cout << std::hex << to_uint( cut_keys[i] ) << std::dec ;
    std::cout << " , weight/mean = " ;
    std::cout << ((float)( weights[i] / mean ));
    std::cout << std::endl ;
  }
}

}

void test_oct_tree_comm_part( ParallelMachine comm , std::istream & is )
{
  typedef OctTreeKey OctTreeKey ;

  unsigned number = 10000 ;

  if ( is.good() ) { is >> number ; }

  const unsigned max = ~((unsigned) 0);
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );
  const unsigned global_number = p_size * number ;

  std::vector<float>      weights( number );
  std::vector<OctTreeKey> cut_keys( p_size );
  std::vector<OctTreeKey> keys( number );

  std::vector<float> local_weights( p_size + 1 );
  std::vector<float> global_weights( p_size + 1 );
  //--------------------------------------------------------------------
  // Test uniform distibution of keys and constant weights

  for ( unsigned i = 0 ; i < number ; ++i ) {
    // const unsigned ig = i + p_rank * number ;
    const unsigned ig = i * p_size + p_rank ;
    const double   r  = (0.5 + (double) ig) / ((double) global_number);
    keys[i] = to_key((unsigned)( max * r ) );
    weights[i] = 1.0 ;
  }

  oct_tree_partition_course(
     comm , number , & keys[0] , & weights[0] , & cut_keys[0] );

  proc_weights( comm, cut_keys, keys, weights, local_weights, global_weights );

  if ( p_rank == 0 ) {
    std::cout << "OctTree Uniform distribution test course" << std::endl ;
    print_cuts( cut_keys , global_weights );
  }

  oct_tree_partition_fine(
     comm , number , & keys[0] , & weights[0] , & cut_keys[0] );

  proc_weights( comm, cut_keys, keys, weights, local_weights, global_weights );

  if ( p_rank == 0 ) {
    std::cout << "OctTree Uniform distribution test fine" << std::endl ;
    print_cuts( cut_keys , global_weights );
    std::cout << std::endl ;
  }

  //--------------------------------------------------------------------
  // Test clustered distribution of weights:

  for ( unsigned i = 0 ; i < number ; ++i ) {
    // const unsigned ig = i + p_rank * number ;
    const unsigned ig = i * p_size + p_rank ;
    const double   r  = (0.5 + (double) ig) / ((double) global_number);
    keys[i] = to_key((unsigned)( max * ( r * 0.5 ) ) );
    weights[i] = 1.0 ;
  }

  oct_tree_partition_course(
     comm , number , & keys[0] , & weights[0] , & cut_keys[0] );

  proc_weights( comm, cut_keys, keys, weights, local_weights, global_weights );

  if ( p_rank == 0 ) {
    std::cout << "OctTree Half distribution test course" << std::endl ;
    print_cuts( cut_keys , global_weights );
  }

  oct_tree_partition_fine(
     comm , number , & keys[0] , & weights[0] , & cut_keys[0] );

  proc_weights( comm, cut_keys, keys, weights, local_weights, global_weights );

  if ( p_rank == 0 ) {
    std::cout << "OctTree Half distribution test fine" << std::endl ;
    print_cuts( cut_keys , global_weights );
    std::cout << std::endl ;
  }

  //--------------------------------------------------------------------
  // Test nonuniform distribution of weights:

  for ( unsigned i = 0 ; i < number ; ++i ) {
    // const unsigned ig = i + p_rank * number ;
    const unsigned ig = i * p_size + p_rank ;
    const double   r  = (0.5 + (double) ig) / ((double) global_number);
    keys[i] = to_key((unsigned)( max * r ) );
    weights[i] = ((float)(1 + ig ));
  }

  oct_tree_partition_course(
     comm , number , & keys[0] , & weights[0] , & cut_keys[0] );

  proc_weights( comm, cut_keys, keys, weights, local_weights, global_weights );

  if ( p_rank == 0 ) {
    std::cout << "OctTree Nonlinear distribution test course" << std::endl ;
    print_cuts( cut_keys , global_weights );
  }

  oct_tree_partition_fine(
     comm , number , & keys[0] , & weights[0] , & cut_keys[0] );

  proc_weights( comm, cut_keys, keys, weights, local_weights, global_weights );

  if ( p_rank == 0 ) {
    std::cout << "OctTree Nonuniform distribution test fine" << std::endl ;
    print_cuts( cut_keys , global_weights );
    std::cout << std::endl ;
  }
}

void test_oct_tree_part_course( ParallelMachine comm , std::istream & is )
{
  typedef OctTreeKey OctTreeKey ;

  unsigned number = 10000 ;
  unsigned ncuts  = 16 ;

  if ( is.good() ) { is >> number ; is >> ncuts ; }

  const unsigned max = ~((unsigned) 0);
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );
  const unsigned global_number = p_size * number ;

  std::vector<float>      weights( number );
  std::vector<OctTreeKey> cut_keys( ncuts );
  std::vector<OctTreeKey> keys( number );

  std::vector<float> local_weights( p_size + 1 );
  std::vector<float> global_weights( p_size + 1 );
  //--------------------------------------------------------------------
  // Test uniform distibution of keys and constant weights

  for ( unsigned i = 0 ; i < number ; ++i ) {
    // const unsigned ig = i + p_rank * number ;
    const unsigned ig = i * p_size + p_rank ;
    const double   r  = (0.5 + (double) ig) / ((double) global_number);
    keys[i] = to_key((unsigned)( max * r ) );
    weights[i] = 1.0 ;
  }

  oct_tree_partition_course(
     comm , number , & keys[0] , & weights[0] , & cut_keys[0] , ncuts );

  proc_weights( comm, cut_keys, keys, weights, local_weights, global_weights );

  if ( p_rank == 0 ) {
    std::cout << "OctTree Uniform distribution test course" << std::endl ;
    print_cuts( cut_keys , global_weights );
  }

  //--------------------------------------------------------------------
  // Test clustered distribution of weights:

  for ( unsigned i = 0 ; i < number ; ++i ) {
    // const unsigned ig = i + p_rank * number ;
    const unsigned ig = i * p_size + p_rank ;
    const double   r  = (0.5 + (double) ig) / ((double) global_number);
    keys[i] = to_key((unsigned)( max * ( r * 0.5 ) ) );
    weights[i] = 1.0 ;
  }

  oct_tree_partition_course(
     comm , number , & keys[0] , & weights[0] , & cut_keys[0] , ncuts );

  proc_weights( comm, cut_keys, keys, weights, local_weights, global_weights );

  if ( p_rank == 0 ) {
    std::cout << "OctTree Half distribution test course" << std::endl ;
    print_cuts( cut_keys , global_weights );
  }

  //--------------------------------------------------------------------
  // Test nonuniform distribution of weights:

  for ( unsigned i = 0 ; i < number ; ++i ) {
    // const unsigned ig = i + p_rank * number ;
    const unsigned ig = i * p_size + p_rank ;
    const double   r  = (0.5 + (double) ig) / ((double) global_number);
    keys[i] = to_key((unsigned)( max * r ) );
    weights[i] = ((float)(1 + ig ));
  }

  oct_tree_partition_course(
     comm , number , & keys[0] , & weights[0] , & cut_keys[0] , ncuts );

  proc_weights( comm, cut_keys, keys, weights, local_weights, global_weights );

  if ( p_rank == 0 ) {
    std::cout << "OctTree Nonlinear distribution test course" << std::endl ;
    print_cuts( cut_keys , global_weights );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void fill_weights(
  const bool uniform ,
  OctTreeKey k , const unsigned d , float * const w )
{
  const unsigned d1 = k.depth() + 1 ;

  if ( k.depth() < d ) {
    for ( unsigned i = 1 ; i <= 8 ; ++i ) {
      k.set_index( d1 , i );
      fill_weights( uniform , k , d , w );
    }
  }
  else {
    const unsigned ord = 1 + 2 * oct_tree_offset( d , k );
    if ( uniform ) {
      w[ord] = 64 ;
    }
    else {
      w[ord] = (float) ord ;
    }
  }
}

template<unsigned Depth>
void oct_tree_part(
  const unsigned np , const double tolerance , const bool uniform )
{
  enum { tree_size = OctTreeSize<Depth>::value };
  enum { tree_size_2 = tree_size * 2 };

  float weights[ tree_size_2 ];
  for ( unsigned i = 0 ; i < tree_size_2 ; ++i ) { weights[i] = 0 ; }

  std::vector<OctTreeKey> cuts( np );

  // Fill lowest level with uniform weights

  fill_weights( uniform , OctTreeKey() , Depth , weights );

  oct_tree_partition_private( 0 , np , Depth ,
                              tolerance , weights , np , & cuts[0] );

  const double max = ((double) (~((unsigned)0))) + 1 ;
  const double inc = max / np ;

  std::cout << "Test oct_tree_part<" << ((unsigned)Depth)
            << ">( " << np << " , " << tolerance << " )"
            << " tree_size = " << ((unsigned)tree_size)
            << std::endl ;
  for ( unsigned i = 0 ; i < np ; ++i ) {
    std::cout << "  P" << i << " begin = " ;
    std::cout << std::hex << to_uint( cuts[i] ) << std::dec ;
    std::cout << " versus " ;
    std::cout << std::hex << ((unsigned)( inc * i )) << std::dec ;
    std::cout << std::endl ;
  }
  std::cout << std::endl ;
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void test_global_box( ParallelMachine comm , std::istream & is )
{
  unsigned number = 100 ;

  if ( is.good() ) { is >> number ; }

  const float global_zero[3] = { -1 , -2 , -3 };
  const float global_size[3] = {  3 ,  7 , 11 };

  // const unsigned p_size = phdmesh::parallel_machine_size( comm );
  const unsigned p_rank = phdmesh::parallel_machine_rank( comm );

  std::vector<IdentProcBox> boxes( number );

  for ( unsigned i = 0 ; i < number ; ++i ) {
    IdentProcBox & obj = boxes[i] ;
    const unsigned n = p_rank * number ;
    const double r  = ((double) n + i )     / ((double) number );
    const double r1 = ((double) n + i + 1 ) / ((double) number );
    obj.ident = p_rank * number + i ;
    obj.proc  = p_rank ;
    obj.part  = 0 ;
    obj.box[0] = (float)( global_zero[0] + global_size[0] * r );
    obj.box[1] = (float)( global_zero[1] + global_size[1] * r );
    obj.box[2] = (float)( global_zero[2] + global_size[2] * r );
    obj.box[3] = (float)( global_zero[0] + global_size[0] * r1 );
    obj.box[4] = (float)( global_zero[1] + global_size[1] * r1 );
    obj.box[5] = (float)( global_zero[2] + global_size[2] * r1 );
  }

  float global_box[6] ;

  box_global_bounds( comm , number , & boxes[0] , 0 , NULL , global_box );

  if ( p_rank == 0 ) {
    std::ios::fmtflags savefmt = std::cout.flags();

    std::cout.precision(15);
    std::cout << "Test global_box bounds" << std::endl ;
    std::cout << "  lower = " << (double) global_box[0]
              << " , "        << (double) global_box[1]
              << " , "        << (double) global_box[2]
              << std::endl ;
    std::cout << "  upper = " << (double) global_box[3]
              << " , "        << (double) global_box[4]
              << " , "        << (double) global_box[5]
              << std::endl
              << std::endl ;

    std::cout.flags(savefmt);
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void test_oct_tree_global_search( ParallelMachine comm , std::istream & is )
{
  unsigned number = 10000 ;

  if ( is.good() ) { is >> number ; }

  const float global_zero[3] = { -1 , -2 , -3 };
  const float global_size[3] = {  3 ,  7 , 11 };

  const unsigned p_size = phdmesh::parallel_machine_size( comm );
  const unsigned p_rank = phdmesh::parallel_machine_rank( comm );

  if ( p_rank == 0 ) {
    std::cout << "Test proximity" << std::endl ;
  }
  
  std::vector<IdentProcBox> boxes( number );

  for ( unsigned i = 0 ; i < number ; ++i ) {
    IdentProcBox & obj = boxes[i] ;
    const unsigned n = p_rank * number ;
    const double r  = ((double) n + i )     / ((double) number );
    const double r1 = ((double) n + i + 1 ) / ((double) number );
    obj.ident = p_rank * number + i ;
    obj.proc  = p_rank ;
    obj.part  = 0 ;

    obj.box[0] = (float)( global_zero[0] + global_size[0] * r );
    obj.box[1] = (float)( global_zero[1] + global_size[1] * r );
    obj.box[2] = (float)( global_zero[2] + global_size[2] * r );
    obj.box[3] = (float)( global_zero[0] + global_size[0] * r1 );
    obj.box[4] = (float)( global_zero[1] + global_size[1] * r1 );
    obj.box[5] = (float)( global_zero[2] + global_size[2] * r1 );
  }

  std::vector< std::pair<IdentProc,IdentProc> > prox ;

  float global_box[6];

  box_global_bounds( comm , number , & boxes[0] , 0 , NULL , global_box );

  oct_tree_proximity_search( comm ,
                             global_box ,
                             number , & boxes[0] ,
                             number , & boxes[0] ,
                             NULL , prox );

  // Take turns outputting

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    parallel_machine_barrier( comm );
    if ( p_rank == p ) {
      std::vector< std::pair<IdentProc,IdentProc> >::iterator j ;
      for ( j = prox.begin() ; j != prox.end() ; ++j ) {
        std::cout << "  P" << p_rank
                  << ": ( ( " << j->first.ident
                  << " , P" << j->first.proc
                  << " ) , ( " << j->second.ident
                  
                  << " ) )" << std::endl ;
      }
      std::cout << std::endl ;
      std::cout.flush();
    }
    parallel_machine_barrier( comm );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void test_oct_tree_global_search_time(
  phdmesh::ParallelMachine comm ,
  std::istream & )
{
  parallel_machine_barrier( comm );

  const float global_zero[3] = { -1 , -2 , -3 };
  const float global_size[3] = {  3 ,  7 , 11 };

  const unsigned p_size = phdmesh::parallel_machine_size( comm );
  const unsigned p_rank = phdmesh::parallel_machine_rank( comm );

  if ( p_rank == 0 ) {
    std::cout << "Time proximity" << std::endl ;
  }

  const unsigned nloop[4] = { 1000 , 100 , 10 , 1 };
  const unsigned count[4] = { 100 , 1000 , 10000 , 100000 };

  float global_box[6] ;

  unsigned stats[6];
  
  for ( unsigned j = 0 ; j < 4 ; ++j ) {
    const unsigned number = count[j] ;
    const unsigned loops  = nloop[j] ;

    std::vector<IdentProcBox> boxes( number );

    unsigned num_local = 0 ;

    for ( unsigned i = 0 ; i < number ; ++i ) {

      if ( i % p_size == p_rank ) {

        IdentProcBox & obj = boxes[ num_local ] ; ++num_local ;

        const double r  = ((double) number + i )     / ((double) number );
        const double r1 = ((double) number + i + 1 ) / ((double) number );

        obj.ident = i ;
        obj.proc  = p_rank ;
        obj.part  = 0 ;
        obj.box[0] = (float)( global_zero[0] + global_size[0] * r );
        obj.box[1] = (float)( global_zero[1] + global_size[1] * r );
        obj.box[2] = (float)( global_zero[2] + global_size[2] * r );
        obj.box[3] = (float)( global_zero[0] + global_size[0] * r1 );
        obj.box[4] = (float)( global_zero[1] + global_size[1] * r1 );
        obj.box[5] = (float)( global_zero[2] + global_size[2] * r1 );
      }
    }

    // Timing loop

    const double time_start = wall_time();

    for ( unsigned k = 0 ; k < loops ; ++k ) {
      std::vector< std::pair<IdentProc,IdentProc> > prox ;

      box_global_bounds( comm, num_local, & boxes[0], 0, NULL, global_box );

      oct_tree_proximity_search( comm ,
                                 global_box ,
                                 num_local , & boxes[0] ,
                                 num_local , & boxes[0] ,
                                 NULL , prox , stats );
    }
    const double time_end = wall_time();

    const double total = ( time_end - time_start ) / loops ;

    if ( p_rank == 0 ) {
      std::cout << "  Time = " << total
                << " for " << number << " boxes "
                << " on " << p_size << " processors"
                << " ; Search "
                << " cell-min = " << stats[0]
                << " cell-max = " << stats[1]
                << " cell-avg = " << stats[2]
                << std::endl ;
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void test_oct_tree( ParallelMachine , std::istream & )
{
  const double tolerance = 1.0e-8 ;
  test_oct_tree_keys();
  oct_tree_part<4>(7,tolerance,true);
  oct_tree_part<5>(7,tolerance,true);
  oct_tree_part<4>(8,tolerance,true);
  oct_tree_part<5>(8,tolerance,true);
  oct_tree_part<4>(13,tolerance,true);
  oct_tree_part<5>(13,tolerance,true);
  oct_tree_part<6>(13,tolerance,true);
}

