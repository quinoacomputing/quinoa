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
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   November 2006
 */

#ifndef util_ParallelReduce_hpp
#define util_ParallelReduce_hpp

#include <cstddef>
#include <iosfwd>
#include <string>
#include <util/Parallel.hpp>
#include <util/SimpleArrayOps.hpp>

//------------------------------------------------------------------------

namespace phdmesh {

/** Write string from any or all processors
 *  to the ostream on the root processor.
 */
void all_write_string( ParallelMachine ,
                       std::ostream & ,
                       const std::string & );

void all_reduce_sum( ParallelMachine ,
                     const double * local , double * global , unsigned count );

void all_reduce_sum( ParallelMachine ,
                     const float * local , float * global , unsigned count );

void all_reduce_sum( ParallelMachine ,
                     const int * local , int * global , unsigned count );

void all_reduce_bor( ParallelMachine ,
                     const unsigned * local ,
                     unsigned * global , unsigned count );

/** Aggregated parallel in-place reduce-to-all-processors operations.
 *
 *  example:
 *    ParallelMachine comm = ... ;
 *    double a[5] ;
 *    int    b[3] ;
 *    all_reduce( comm , Sum<5>( a ) , Max<3>( b ) , ... );
 *
 *  Reduction options include Sum, Prod, Max, Min, BitOr, and BitAnd
 */

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace phdmesh {

extern "C" {
typedef void (*ParallelReduceOp)
  ( void * inv , void * outv , int * , ParallelDatatype * );
}

void all_reduce_internal( ParallelMachine  arg_comm ,
                          ParallelReduceOp arg_op ,
                          void           * arg_in ,
                          void           * arg_out ,
                          unsigned         arg_len );

namespace {

// Blank namespace so that this class produces local symbols,
// avoiding complaints from a linker of multiple-define symbols.

struct ReduceEnd {
  struct BufferType {};
  void copyin(  BufferType & ) const {}
  void copyout( BufferType & ) const {}
  static void op( BufferType & , BufferType & ) {}
};

// Workhorse class for aggregating reduction operations.

template < class Oper , class Next = ReduceEnd >
struct Reduce {
  typedef typename Oper::type Type ;
  enum { N = Oper::N };

  struct BufferType {
    Type                      m_value[N];
    typename Next::BufferType m_next ;
  };

  Next   m_next ;
  Type * m_ptr ;

  Next & set( const Oper & arg ) { m_ptr = arg.ptr ; return m_next ; }

  void reduce( ParallelMachine comm ) const ;

  void copyin( BufferType & b ) const
    { Copy<N>( b.m_value , m_ptr ); m_next.copyin( b.m_next ); }

  void copyout( BufferType & b ) const
    { Copy<N>( m_ptr , b.m_value ); m_next.copyout( b.m_next ); }

  static void op( BufferType & dst , BufferType & src )
    { Oper::op(dst.m_value,src.m_value); Next::op(dst.m_next,src.m_next); }

  static void void_op( void*inv, void*inoutv, int*, ParallelDatatype*);
};

template <class Oper, class Next>
void Reduce<Oper,Next>::void_op( void*inv, void*inoutv,int*,ParallelDatatype*)
{
  op( * reinterpret_cast<BufferType*>( inoutv ) ,
      * reinterpret_cast<BufferType*>( inv ) );
}

template <class Oper, class Next>
void Reduce<Oper,Next>::reduce( ParallelMachine comm ) const
{
  ParallelReduceOp f = reinterpret_cast<ParallelReduceOp>( & void_op );
  BufferType inbuf , outbuf ;
  copyin( inbuf );
  all_reduce_internal( comm , f , & inbuf , & outbuf , sizeof(BufferType) );
  copyout( outbuf );
}

} // namespace
} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace phdmesh {

template < class Op1 >
inline
void all_reduce( ParallelMachine comm , const Op1 & op1 )
{
  Reduce< Op1 > work ;
  work.set( op1 );
  work.reduce( comm );
}

template < class Op1 , class Op2 >
inline
void all_reduce( ParallelMachine comm , const Op1 & op1 ,
                                        const Op2 & op2 )
{
  Reduce< Op1 ,
  Reduce< Op2 > > work ;
  work.set( op1 ).set( op2 );
  work.reduce( comm );
}

template < class Op1 , class Op2 , class Op3 >
inline
void all_reduce( ParallelMachine comm , const Op1 & op1 ,
                                        const Op2 & op2 ,
                                        const Op3 & op3 )
{
  Reduce< Op1 ,
  Reduce< Op2 ,
  Reduce< Op3 > > > work ;
  work.set( op1 ).set( op2 ).set( op3 );
  work.reduce( comm );
}

template < class Op1 , class Op2 , class Op3 , class Op4 >
inline
void all_reduce( ParallelMachine comm , const Op1 & op1 ,
                                        const Op2 & op2 ,
                                        const Op3 & op3 ,
                                        const Op4 & op4 )
{
  Reduce< Op1 ,
  Reduce< Op2 ,
  Reduce< Op3 ,
  Reduce< Op4 > > > > work ;
  work.set( op1 ).set( op2 ).set( op3 ).set( op4 );
  work.reduce( comm );
}

template < class Op1 , class Op2 , class Op3 , class Op4 ,
           class Op5 >
inline
void all_reduce( ParallelMachine comm , const Op1 & op1 ,
                                        const Op2 & op2 ,
                                        const Op3 & op3 ,
                                        const Op4 & op4 ,
                                        const Op5 & op5 )
{
  Reduce< Op1 ,
  Reduce< Op2 ,
  Reduce< Op3 ,
  Reduce< Op4 ,
  Reduce< Op5 > > > > > work ;
  work.set( op1 ).set( op2 ).set( op3 ).set( op4 ).set( op5 );
  work.reduce( comm );
}

template < class Op1 , class Op2 , class Op3 , class Op4 ,
           class Op5 , class Op6 >
inline
void all_reduce( ParallelMachine comm , const Op1 & op1 ,
                                        const Op2 & op2 ,
                                        const Op3 & op3 ,
                                        const Op4 & op4 ,
                                        const Op5 & op5 ,
                                        const Op6 & op6 )
{
  Reduce< Op1 ,
  Reduce< Op2 ,
  Reduce< Op3 ,
  Reduce< Op4 ,
  Reduce< Op5 ,
  Reduce< Op6 > > > > > > work ;
  work.set( op1 ).set( op2 ).set( op3 ).set( op4 ).set( op5 ).set( op6 );
  work.reduce( comm );
}

}

//----------------------------------------------------------------------

#endif

