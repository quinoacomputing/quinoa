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

#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <util/TPI.hpp>
#include <util/Parallel.hpp>
#include <util/ParallelReduce.hpp>
#include <util/ParallelInputStream.hpp>

#include <txblas/Reduction.hpp>

using namespace phdmesh ;

namespace {

struct FillWork {
  double   x_mag ;
  double * x_beg ;
  unsigned x_length ;
};

void task_rand_fill( void * arg , TPI::ThreadPool pool )
{
  const double r_max = RAND_MAX ;
  int p_rank , p_size ; TPI::Rank( pool , p_rank , p_size );

  FillWork & w = * reinterpret_cast<FillWork*>( arg );

  const unsigned p_next = p_rank + 1 ;
  const unsigned       n = w.x_length ;
        double * const xe = w.x_beg + ( n * p_next ) / p_size ;
        double *       x  = w.x_beg + ( n * p_rank ) / p_size ;

  const double mag = w.x_mag ;

  unsigned seed = p_rank ;

  for ( ; xe != x ; ++x ) {
    const double r = rand_r( & seed );
    *x = mag * 2.0 * ( ( r / r_max ) - 0.5 );
  }
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void test_accuracy_txblas( ParallelMachine comm ,
                           const unsigned M ,
                           const double mag )
{
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  const unsigned m_rem   = M % p_size ;
  const unsigned m_local = M / p_size + ( p_rank < m_rem ? 1 : 0 );

  std::vector<double> values( m_local );

  {
    FillWork data ;
    data.x_mag = mag ;
    data.x_beg = & values[0] ;
    data.x_length = m_local ;
    TPI::Run( & task_rand_fill , & data );
  }

  const double init = 1 ;
  double a = 0 ;
  Summation z ;
  if ( p_rank == 0 ) {
    a = init ;
    z += init ;
  }

  std::vector<double>::iterator i ;


  txdsum_add_array( z.xdval , m_local , & values[0] );

  for ( i = values.begin() ; i != values.end() ; ++i ) { a += *i ; }

  all_reduce( comm , Sum<1>( & z ) , Sum<1>( & a ) );

  for ( i = values.begin() ; i != values.end() ; ++i ) { *i = - *i ; }

  txdsum_add_array( z.xdval , m_local , & values[0] );

  for ( i = values.begin() ; i != values.end() ; ++i ) { a += *i ; }

  all_reduce( comm , Sum<1>( & z ) , Sum<1>( & a ) );


  const double z_val = z.value();

  if ( p_rank == 0 ) {
    std::cout << "ACCURACY txdsum_add_array[ N = "
              << M << " , mag = " << mag << " ] " << std::endl
              << "  init_value[ " << init
              << " ] , txdsum_add[ " << z_val
              << " ] , simple_sum[ " << a << " ]" << std::endl ;
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void test_accuracy( phdmesh::ParallelMachine comm ,
                    std::istream & is )
{
  unsigned num = 100000000 ;
  double mag = 1e10 ;
  if ( is.good() ) { is >> num ; }
  if ( is.good() ) { is >> mag ; }
  test_accuracy_txblas( comm , num , mag );
}

//----------------------------------------------------------------------

void test_reduce( ParallelMachine comm , std::istream & is )
{
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  unsigned M = 1000000 ;
  double mag = 1e10 ;

  if ( is.good() ) { is >> M ; }
  if ( is.good() ) { is >> mag ; }

  //--------------------------------------------------------------------
  {
    srand( p_rank ); // each processor do a different series

    const double r_max = RAND_MAX ;

    std::vector<double> values( M );

    for ( std::vector<double>::iterator
          i = values.begin() ; i != values.end() ; ++i ) {
      const double r = rand() ;
      const double x = 2.0 * ( ( r / r_max ) - 0.5 );
      *i = x * mag ;
    }

    std::vector<double> values_neg( values );

    for ( std::vector<double>::iterator
          i = values_neg.begin() ; i != values_neg.end() ; ++i ) {
      *i = - *i ;
    }

    const double a0 = 1 ;
    double a = a0 ;
    double z[4] = { 0 , 0 , 0 , 0 };
    xdsum_add_value( z , a0 );

    for ( std::vector<double>::iterator
          i = values.begin() ; i != values.end() ; ++i ) {
      a += *i ;
    }

    txdsum_add_array( z , M , & values[0] );

    for ( std::vector<double>::iterator
          i = values_neg.begin() ; i != values_neg.end() ; ++i ) {
      a += *i ;
    }

    txdsum_add_array( z , M , & values_neg[0] );

    xdsum_get_value( z , z );

    const int local_flag = z[0] < a0 || a0 < z[0] ;
    int flag = local_flag ;

    all_reduce( comm , Max<1>( & flag ) );

    if ( flag ) {
      if ( local_flag ) {
        std::cout << "P" << p_rank << " SUMMATION ERROR = "
                  << ( z[0] - a0 ) << std::endl ;
      }
      throw std::runtime_error( std::string( "FAIL SUMMATION TEST 1" ) );
    }

    if ( p_rank == 0 ) {
      std::cout << "PASS SUMMATION TEST: P0 has "
                << z[0] << " vs. " << a << std::endl ;
    }
  }
  //--------------------------------------------------------------------

  double dval[5] = { 0 , 1 , 2 , 3 , 0 };
  long   ival[3] = { 0 , 0 , 0 };
  Summation aval ;
  unsigned flag  = 0 ;

  ival[1] = p_rank + 1 ;
  ival[2] = - ((long)( p_rank + 1 ));
  dval[4] = p_rank + 1 ;
  aval += dval[1] ;

  all_reduce( comm , Sum<5>( dval ) ,
                     Max<3>( ival ) ,
                     Sum<1>( & aval ) ,
                     BitOr<1>( & flag ) );

  const double sum_ranks = ( p_size * ( p_size + 1 ) ) / 2 ;

  if ( (unsigned) dval[0] != 0 ||
       (unsigned) dval[1] != p_size ||
       (unsigned) dval[2] != p_size * 2 ||
       (unsigned) dval[3] != p_size * 3 ||
       (unsigned) dval[4] != (unsigned) sum_ranks ||
       ival[0] != 0 ||
       ival[1] != ( (long) p_size ) ||
       ival[2] != -1 ||
       (unsigned) ( aval.value() ) != p_size ||
       flag != 0 ) {
    std::cout << "P" << p_rank << " "
              << dval[0] << " " << dval[1] << " "
              << dval[2] << " " << dval[3] << " "
              << dval[4] << " , "
              << ival[0] << " " << ival[1] << " "
              << ival[2] << " , "
              << aval.value() << " , "
              << flag << std::endl ;
    throw std::runtime_error( std::string( "FAIL REDUCE TEST 1" ) );
  }

  parallel_machine_barrier( comm );

  if ( p_rank == 0 ) {
    std::cout << "PASS REDUCE TEST 1 for NP = " << p_size << " : "
              << dval[0] << " " << dval[1] << " "
              << dval[2] << " " << dval[3] << " "
              << dval[4] << " , "
              << ival[0] << " " << ival[1] << " "
              << ival[2] << " , " << flag << std::endl ;
    std::cout.flush();
  }
}

