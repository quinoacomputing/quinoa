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

#ifndef SPARSELINEARSYSTEM_HPP
#define SPARSELINEARSYSTEM_HPP

#include <cmath>
#include <impl/KokkosArray_Timer.hpp>
#include <KokkosArray_Value.hpp>
#include <KokkosArray_MultiVector.hpp>
#include <KokkosArray_CrsArray.hpp>

namespace KokkosArray {
namespace Impl {

template< class Scalar , class DeviceType , class > struct Dot ;
template< class Scalar , class DeviceType > struct WAXPBY ;
template< class Scalar , class DeviceType > struct AXPBY ;
template< class Scalar , class DeviceType > struct FILL ;
template< class AType , class XType , class YType > struct Multiply ;

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

//----------------------------------------------------------------------------

template< typename ScalarType , class Device >
struct CrsMatrix {
  typedef Device      device_type ;
  typedef ScalarType  value_type ;

  typedef CrsArray< int , device_type , int >     graph_type ;
  typedef MultiVector< value_type , device_type > coefficients_type ;

  graph_type         graph ;
  coefficients_type  coefficients ;
};

//----------------------------------------------------------------------------

template< typename Scalar , class Device >
void waxpby( const ParallelDataMap & data_map ,
             const double alpha ,
             const MultiVector< Scalar , Device > & x ,
             const double beta ,
             const MultiVector< Scalar , Device > & y ,
             const MultiVector< Scalar , Device > & w )
{
  if ( y == w ) {
    Impl::AXPBY<Scalar,Device>::apply( data_map.count_owned , alpha , x , beta , y );
  }
  else {
    Impl::WAXPBY<Scalar,Device>::apply( data_map.count_owned , alpha , x , beta , y , w );
  }
}

template< typename Scalar , class Device >
void fill( const double alpha ,
           const MultiVector< Scalar , Device > & w )
{
  Impl::FILL<Scalar,Device>::apply( w.length() , alpha , w );
}

//----------------------------------------------------------------------------

template< typename Scalar , class Device >
double dot( const ParallelDataMap & data_map ,
            const MultiVector< Scalar , Device > & x ,
            const MultiVector< Scalar , Device > & y ,
            const Value< double , Device > & temp_result )
{
  typedef Impl::Dot< Scalar , Device , Impl::unsigned_<2> > op_type ;

  Impl::Dot< Scalar , Device , Impl::unsigned_<2> >
      ::apply( data_map.count_owned , x , y , temp_result );

  double global_result = 0 ;

  deep_copy( global_result , temp_result );

#if defined( HAVE_MPI )

  double local_result = global_result ;

  MPI_Allreduce( & local_result , & global_result , 1 , MPI_DOUBLE , MPI_SUM ,
                 data_map.machine.mpi_comm );

#endif

  return global_result ;
}

template< typename Scalar , class Device >
double dot( const ParallelDataMap & data_map ,
            const MultiVector< Scalar , Device > & x ,
            const MultiVector< Scalar , Device > & y )
{
  typedef Value< double , Device > value_type ;
  value_type temp_result = create_value< value_type >();
  return dot( data_map , x , y , temp_result );
}

template< typename Scalar , class Device >
double dot( const ParallelDataMap & data_map ,
            const MultiVector< Scalar , Device > & x ,
            const Value< double , Device > & temp_result )
{
  Impl::Dot< Scalar , Device , Impl::unsigned_<1> >
      ::apply( data_map.count_owned , x , temp_result );

  double global_result = 0 ;

  deep_copy( global_result , temp_result );

#if defined( HAVE_MPI )

  double local_result = global_result ;

  MPI_Allreduce( & local_result , & global_result , 1 , MPI_DOUBLE , MPI_SUM ,
                 data_map.machine.mpi_comm );

#endif

  return global_result ;
}

template< typename Scalar , class Device >
double dot( const ParallelDataMap & data_map ,
            const MultiVector< Scalar , Device > & x )
{
  typedef Value< double , Device > value_type ;
  value_type temp_result = create_value< value_type >();
  return dot( data_map , x , temp_result );
}

//----------------------------------------------------------------------------

template< typename AScalarType ,
          typename VScalarType ,
          class Device >
class Operator {
  typedef CrsMatrix<AScalarType,Device>    matrix_type ;
  typedef MultiVector<VScalarType,Device>  vector_type ;

private:
  const CrsMatrix<AScalarType,Device> A ;

#if defined( HAVE_MPI )

  ParallelDataMap                                         data_map ;
  AsyncExchange< VScalarType , Device , ParallelDataMap > exchange ;

public:

  Operator( const ParallelDataMap                  & arg_data_map ,
            const CrsMatrix<AScalarType,Device>    & arg_A )
    : A( arg_A )
    , data_map( arg_data_map )
    , exchange( arg_data_map , 1 )
    {}

#else /* ! defined( HAVE_MPI ) */

public:

  Operator( const ParallelDataMap                  & ,
            const CrsMatrix<AScalarType,Device>    & arg_A )
    : A( arg_A ) {}

#endif

  void apply( const MultiVector<VScalarType,Device>  & x ,
              const MultiVector<VScalarType,Device>  & y )
  {
#if defined( HAVE_MPI )
    // Gather off-processor data for 'x'

    PackArray< vector_type >::pack( exchange.buffer() ,
                                    data_map.count_interior ,
                                    data_map.count_send , x );

    exchange.setup();

    // If interior & boundary matrices then could launch interior multiply

    exchange.send_receive();

    UnpackArray< vector_type >::unpack( x , exchange.buffer() ,
                                        data_map.count_owned ,
                                        data_map.count_receive );
    const typename Device::size_type nrow = data_map.count_owned ;
    const typename Device::size_type ncol = data_map.count_owned +
                                            data_map.count_receive ;
#else /* ! defined( HAVE_MPI ) */

    const typename Device::size_type nrow = A.graph.row_map.length();
    const typename Device::size_type ncol = A.graph.row_map.length();

#endif


    Impl::Multiply<matrix_type,vector_type,vector_type>
      ::apply( A , nrow , ncol , x , y );
  }
};

//----------------------------------------------------------------------------

template< typename AScalarType , typename VScalarType , class Device >
void cgsolve(
  const ParallelDataMap                 data_map ,
  const CrsMatrix<AScalarType,Device>   A ,
  const MultiVector<VScalarType,Device> b ,
  const MultiVector<VScalarType,Device> x ,
  size_t & iteration ,
  double & normr ,
  double & iter_time ,
  const size_t maximum_iteration = 200 ,
  const double tolerance = std::numeric_limits<VScalarType>::epsilon() )
{
  typedef MultiVector<VScalarType,Device> vector_type ;
  typedef Value      <VScalarType,Device> value_type ;

  const size_t count_owned = data_map.count_owned ;
  const size_t count_total = data_map.count_owned + data_map.count_receive ;

  Operator<AScalarType,VScalarType,Device> matrix_operator( data_map , A );

  vector_type r  = create_multivector< vector_type >( "cg::r" , count_owned );
  vector_type p  = create_multivector< vector_type >( "cg::p" , count_total );
  vector_type Ap = create_multivector< vector_type >( "cg::Ap", count_owned );
  value_type dot_temp = create_value< value_type >( "cg::dot_temp" );

  /* p  = x      */ deep_copy( p , x , count_owned );
  /* Ap = A * p  */ matrix_operator.apply( p , Ap );
  /* r  = b - Ap */ waxpby( data_map , 1.0 , b , -1.0 , Ap , r );
  /* p  = r      */ deep_copy( p , r , count_owned );

  double old_rdot = dot( data_map , r , dot_temp );

  normr     = sqrt( old_rdot );
  iteration = 0 ;

  KokkosArray::Impl::Timer wall_clock ;

  while ( tolerance < normr && iteration < maximum_iteration ) {

    /* Ap = A * p  */ matrix_operator.apply( p , Ap );

    const double pAp_dot = dot( data_map , p , Ap , dot_temp );
    const double alpha   = old_rdot / pAp_dot ;

    /* x += alpha * p ;  */ waxpby( data_map,  alpha, p , 1.0 , x , x );
    /* r -= alpha * Ap ; */ waxpby( data_map, -alpha, Ap, 1.0 , r , r );

    const double r_dot = dot( data_map , r , dot_temp );
    const double beta  = r_dot / old_rdot ;

    /* p = r + beta * p ; */ waxpby( data_map , 1.0 , r , beta , p , p );

    normr = sqrt( old_rdot = r_dot );
    ++iteration ;
  }

  iter_time = wall_clock.seconds();
}

//----------------------------------------------------------------------------

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef SPARSELINEARSYSTEM_HPP */

