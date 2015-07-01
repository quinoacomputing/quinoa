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
 * @date   August 2007
 */

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include <util/ParallelComm.hpp>

#include <txblas/CR4Matrix.hpp>
#include <txblas/cr4_mxv.h>

namespace phdmesh {

void simple_partition( ParallelMachine comm ,
                       const unsigned nglobal ,
                       std::vector<unsigned> & partition )
{
  const unsigned p_size = parallel_machine_size( comm );

  // Partition rows evenly among processors

  const unsigned min_row_per_proc = nglobal / p_size ;
  const unsigned remaining_rows   = nglobal % p_size ;

  partition.resize( p_size + 1 );

  partition[0] = 0 ;
  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    const unsigned n = min_row_per_proc + ( p < remaining_rows ? 1 : 0 );
    partition[p+1] = partition[p] + n ;
  }
}

//----------------------------------------------------------------------

namespace {

#if defined( HAVE_MPI )

#define PARALLEL_DATATYPE_DOUBLE   MPI_DOUBLE
#define PARALLEL_DATATYPE_UNSIGNED MPI_UNSIGNED

void all_to_all( ParallelMachine comm ,
                 ParallelDatatype type ,
                 const bool sparse ,
                 void * send_buf , const int * send_disp ,
                 void * recv_buf , const int * recv_disp )
{
  static const char method[] = "txblas_SparseMatrix.cpp(all_to_all)" ;

  std::ostringstream msg ;

  const int mpi_tag = 0 ;
  const unsigned p_rank = parallel_machine_rank( comm );
  const unsigned p_size = parallel_machine_size( comm );

  if ( sparse ) {
    int type_size = 0 ;
    MPI_Type_size( type , & type_size );

    unsigned num_recv = 0 ;
    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      if ( p != p_rank && recv_disp[p] < recv_disp[p+1] ) { ++num_recv ; }
    }

    // Post receives for specific processors with specific sizes

    MPI_Request request_null = MPI_REQUEST_NULL ;
    std::vector<MPI_Request> request( num_recv , request_null );
    std::vector<MPI_Status>  status(  num_recv );

    int result = MPI_SUCCESS ;
    unsigned count = 0 ;

    for ( unsigned p = 0 ; result == MPI_SUCCESS && p < p_size ; ++p ) {
      const unsigned recv_size = recv_disp[p+1] - recv_disp[p] ;
      if ( p != p_rank && recv_size ) {
        void * ptr = ((char *) recv_buf) + type_size * recv_disp[p] ;
        result = MPI_Irecv( ptr , recv_size , type ,
                            p , mpi_tag , comm , & request[count] );
        ++count ;
      }
    }

    if ( MPI_SUCCESS != result ) {
      msg << method << " LOCAL[" << p_rank << "] ERROR: "
          << result << " == MPI_Irecv , " ;
    }

    //------------------------------
    // Sync to allow ready sends and for a potential error

    int local_error = MPI_SUCCESS == result ? 0 : 1 ;
    int global_error = 0 ;

    result = MPI_Allreduce( & local_error , & global_error ,
                            1 , MPI_INT , MPI_SUM , comm );

    if ( MPI_SUCCESS != result ) {
      msg << method << " GLOBAL ERROR: " << result << " == MPI_Allreduce" ;
    }
    else if ( global_error ) {
      result = MPI_ERR_UNKNOWN ;
    }
    else {
      //------------------------------
      // Ready-send the buffers, rotate the send processor
      // in a simple attempt to smooth out the communication traffic.

      for ( unsigned p = 0 ; MPI_SUCCESS == result && p < p_size ; ++p ) {
        const int p_dst = ( p + p_rank ) % p_size ;
        const int send_size = send_disp[ p_dst + 1 ] - send_disp[ p_dst ] ;
        if ( ( p_dst != (int) p_rank ) && send_size ) {
          void * ptr = ((char *) send_buf) + type_size * send_disp[ p_dst ] ;
          result = MPI_Rsend( ptr, send_size, type, p_dst, mpi_tag, comm );
        }
      }

      if ( MPI_SUCCESS != result ) {
        msg << method << " LOCAL ERROR: " << result << " == MPI_Rsend , " ;
      }
      else {
        MPI_Request * const p_request = & request[0] ;
        MPI_Status  * const p_status  = & status[0] ;

        result = MPI_Waitall( num_recv , p_request , p_status );
      }
    }
  }
  else {
    std::vector<int> send_size( p_size );
    std::vector<int> recv_size( p_size );

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      send_size[p] = send_disp[p+1] - send_disp[p] ;
      recv_size[p] = recv_disp[p+1] - recv_disp[p] ;
    }
    send_size[p_rank] = 0 ;
    recv_size[p_rank] = 0 ;

    int * sdisp = const_cast<int*>( send_disp );
    int * rdisp = const_cast<int*>( recv_disp );

    MPI_Alltoallv( send_buf , & send_size[0] , sdisp , type ,
                   recv_buf , & recv_size[0] , rdisp , type , comm );

  }
}

#else

#define PARALLEL_DATATYPE_DOUBLE   0
#define PARALLEL_DATATYPE_UNSIGNED 0

void all_to_all( ParallelMachine ,
                 ParallelDatatype ,
                 const bool ,
                 void * , const int * ,
                 void * , const int * )
{}

#endif


void ordered_insert( std::vector<unsigned> & vec , unsigned val )
{
  const std::vector<unsigned>::iterator e = vec.end();
        std::vector<unsigned>::iterator i = vec.begin();

  i = std::lower_bound( i , e , val );

  if ( e == i || val != *i ) { vec.insert( i , val ); }
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

CR_Matrix::~CR_Matrix() {}

//----------------------------------------------------------------------

void CR_Matrix::multiply(
  const double * const x ,
        double * const y ) const
{
  if ( 1 < m_comm_size ) {

    std::vector<double> x_work( m_work_disp[ m_comm_size ] );

    { // Remote gathered portion:
      std::vector<double> x_send( m_send_disp[ m_comm_size ] );

      for ( unsigned i = 0 ; i < m_send_map.size() ; ++i ) {
        x_send[i] = x[ m_send_map[i] ];
      }

      // Skips local-to-local
      all_to_all( m_comm , PARALLEL_DATATYPE_DOUBLE , m_sparse ,
                  & x_send[0] , & m_send_disp[0] ,
                  & x_work[0] , & m_work_disp[0] );
    }

    { // Local portion:
      double * d = & x_work[0] + m_work_disp[ m_comm_rank ];
      const double * s = x ;
      const double * const e = x + m_row_size ;
      while ( e != s ) { *d++ = *s++ ; }
    }

    txblas_cr_mxv( m_row_size ,
                   & m_prefix[0] ,
                   & m_coli[0] ,
                   & m_coef[0] ,
                   & x_work[0] , y );
  }
  else {
    txblas_cr_mxv( m_row_size ,
                   & m_prefix[0] ,
                   & m_coli[0] ,
                   & m_coef[0] ,
                   x , y );
  }
}

//----------------------------------------------------------------------


CR_Matrix::CR_Matrix(
 ParallelMachine arg_comm ,
 const std::vector<unsigned> & arg_partition ,
       std::vector<unsigned> & arg_prefix ,
       std::vector<unsigned> & arg_coli ,
       std::vector<double>   & arg_coef )
  : m_comm( arg_comm ),
    m_comm_size( parallel_machine_size( arg_comm ) ),
    m_comm_rank( parallel_machine_rank( arg_comm ) ),
    m_sparse( false ),
    m_work_disp(),
    m_send_disp(),
    m_send_map(),
    m_row_size( 0 ),
    m_prefix(),
    m_coli(),
    m_coef()
{
  static const char method[] = "phdmesh::CR_Matrix::CR_Matrix" ;

  if ( arg_prefix.empty() ) { return ; }

  //------------------------------------

  if ( arg_coli.size() != arg_prefix.back() ||
       arg_coef.size() != arg_prefix.back() ) {
    std::ostringstream msg ;
    msg << method << " ERROR" ;
    msg << " arg_coli.size() = " << arg_coli.size() ;
    msg << " arg_coef.size() = " << arg_coef.size() ;
    msg << " !=  arg_prefix.back() = " << arg_prefix.back() ;
    throw std::invalid_argument( msg.str() );
  }

  swap( m_prefix , arg_prefix );
  swap( m_coli , arg_coli );
  swap( m_coef , arg_coef );

  m_row_size = m_prefix.size() - 1 ;

  if ( 1 == m_comm_size ) { return ; }

  //------------------------------------

  if ( arg_partition.size() != 1 + m_comm_size ) {
    std::ostringstream msg ;
    msg << method << " ERROR" ;
    msg << " comm_size = " << m_comm_size ;
    msg << " + 1  !=  arg_partition.size() = " << arg_partition.size() ;
    throw std::invalid_argument( msg.str() );
  }

  const unsigned row_first = arg_partition[ m_comm_rank ];
  const unsigned row_end   = arg_partition[ m_comm_rank + 1 ] ;

  if ( m_row_size != ( row_end - row_first ) ) {
    std::ostringstream msg ;
    msg << method << " ERROR" ;
    msg << " arg_prefix'row_size = " << m_row_size ;
    msg << " !=  arg_partition'row_size = " << ( row_end - row_first );
    throw std::invalid_argument( msg.str() );
  }

  //------------------------------------

  m_send_disp.resize( m_comm_size + 1 );
  m_work_disp.resize( m_comm_size + 1 );

  // Generate a vector of off-processor column identifiers

  std::vector<unsigned> work_col_ident ;

  {
    const std::vector<unsigned>::iterator j = m_coli.end();
          std::vector<unsigned>::iterator b = m_coli.begin();
          std::vector<unsigned>::iterator i ;

    for ( i = b ; j != i ; ++i ) { 
      const unsigned global_col = *i ;
      if ( global_col < row_first || row_end <= global_col ) {
        ordered_insert( work_col_ident , global_col );
      }
    }
  }

  //------------------------------------
  // Map column global identifiers to local work offsets

  {
    const std::vector<unsigned>::iterator b = work_col_ident.begin();
    const std::vector<unsigned>::iterator e = work_col_ident.end();
          std::vector<unsigned>::iterator j ;

    j = std::lower_bound( b , e , row_end );

    const unsigned local_row_end = j - b ;

    for ( std::vector<unsigned>::iterator
          i = m_coli.begin() ; i != m_coli.end() ; ++i ) {
      const unsigned global_col = *i ;

      j = std::lower_bound( b, e, global_col );

      unsigned local_col = j - b ;

      if ( row_end <= global_col ) { local_col += local_row_end ; }

      *i = local_col ;
    }
  }

  //------------------------------------
  // Displacement prefix for work vector

  {
    std::vector<unsigned>::const_iterator i = work_col_ident.begin() ;

    m_work_disp[0] = 0 ;

    for ( unsigned p = 0 ; p < m_comm_size ; ++p ) {
      const unsigned p_row_end = arg_partition[p+1] ;
      unsigned count = 0 ;
      for ( ; i != work_col_ident.end() && *i < p_row_end ; ++i ) {
        ++count ;
      }

      m_work_disp[p+1] = m_work_disp[p] + count ;
    }
  }

  //------------------------------------
  // Set up communications to gather work subvector

  {
    std::vector<unsigned> send_col_size( m_comm_size );
    std::vector<unsigned> recv_col_size( m_comm_size );

    for ( unsigned p = 0 ; p < m_comm_size ; ++p ) {
      send_col_size[p] = m_work_disp[p+1] - m_work_disp[p] ;
    }

    if ( send_col_size[ m_comm_rank ] ) {
      std::ostringstream msg ;
      msg << method << " ERROR with communication sizing logic" ;
      throw std::logic_error( msg.str() );
    }

    unsigned num_msg_maximum = 0 ;

    comm_sizes( m_comm , m_comm_size / 4 , num_msg_maximum ,
                & send_col_size[0] , & recv_col_size[0] );

    m_sparse = num_msg_maximum < ( m_comm_size / 4 );

    m_send_disp[0] = 0 ;
    for ( unsigned p = 0 ; p < m_comm_size ; ++p ) {
      m_send_disp[p+1] = m_send_disp[p] + recv_col_size[p] ;
    }
  }

  const unsigned send_map_size = m_send_disp[ m_comm_size ];

  m_send_map.resize( send_map_size );

  all_to_all( m_comm , PARALLEL_DATATYPE_UNSIGNED , m_sparse ,
              & work_col_ident[0] , & m_work_disp[0],
              & m_send_map[0] ,     & m_send_disp[0] );

  //------------------------------------
  // Remap the 'm_work_disp' for receiving coefficients into the
  // work vector: [ lower_row_recv , local_row , upper_row_recv ]

  for ( unsigned p = m_comm_rank ; p < m_comm_size ; ++p ) {
    m_work_disp[p+1] += m_row_size ;
  }

  //------------------------------------
  // Map the send_map from global to local indices,
  // also sanity check it.

  for ( unsigned i = 0 ; i < send_map_size ; ++i ) {

    if ( m_send_map[i] < (int) row_first ||
                         (int) row_end <= m_send_map[i] ) {
      std::ostringstream msg ;
      msg << method << " ERROR Received index " ;
      msg << m_send_map[i] ;
      msg << " out of range [ " ;
      msg << row_first ;
      msg << " : " ;
      msg << row_end ;
      msg << " )" ;
      throw std::runtime_error( msg.str() );
    }

    m_send_map[i] -= row_first ;
  }
}

//----------------------------------------------------------------------

}


