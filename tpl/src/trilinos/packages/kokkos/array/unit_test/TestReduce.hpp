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

#include <gtest/gtest.h>

#ifndef KOKKOS_MACRO_DEVICE
#error "KOKKOS_MACRO_DEVICE undefined"
#endif

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <KokkosArray_Value.hpp>
#include <KokkosArray_ParallelReduce.hpp>

#include <impl/KokkosArray_Preprocessing_macros.hpp>

/*--------------------------------------------------------------------------*/

namespace Test {

template< typename ScalarType , class DeviceType >
class ReduceFunctor ;

template< typename ScalarType >
class ReduceFunctor< ScalarType , KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  struct value_type {
    ScalarType value[3] ;
  };

  const size_type nwork ;

  ReduceFunctor( const size_type & arg_nwork ) : nwork( arg_nwork ) {}

  ReduceFunctor( const ReduceFunctor & rhs )
    : nwork( rhs.nwork ) {}

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void init( value_type & dst )
  {
    dst.value[0] = 0 ;
    dst.value[1] = 0 ;
    dst.value[2] = 0 ;
  }

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & dst ,
                    const volatile value_type & src )
  {
    dst.value[0] += src.value[0] ;
    dst.value[1] += src.value[1] ;
    dst.value[2] += src.value[2] ;
  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( device_type :: size_type iwork , value_type & dst ) const
  {
    dst.value[0] += 1 ;
    dst.value[1] += iwork + 1 ;
    dst.value[2] += nwork - iwork ;
  }

#endif /* defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

};

} // namespace Test

namespace {

template< typename ScalarType , class DeviceType >
class TestReduce ;

template< typename ScalarType >
class TestReduce< ScalarType , KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  typedef Test::ReduceFunctor< ScalarType , device_type > functor_type ;

  typedef typename functor_type::value_type value_type ;

  //------------------------------------

  TestReduce( const size_type & nwork )
  { run_test(nwork); }

  void run_test( const size_type & nwork )
  {
    value_type result ;

    typedef KokkosArray::Value< value_type , device_type > result_type ;

    result_type device_result = KokkosArray::create_value< result_type >();

    KokkosArray::parallel_reduce( nwork , functor_type( nwork ) , device_result );

    KokkosArray::deep_copy( result , device_result );

    const unsigned long nw   = nwork ;
    const unsigned long nsum = nw % 2 ? nw * (( nw + 1 )/2 )
                                      : (nw/2) * ( nw + 1 );

    ASSERT_EQ( result.value[0], (ScalarType) nw);
    ASSERT_EQ( result.value[1], (ScalarType) nsum);
    ASSERT_EQ( result.value[2], (ScalarType) nsum);

  }
};

}

/*--------------------------------------------------------------------------*/

