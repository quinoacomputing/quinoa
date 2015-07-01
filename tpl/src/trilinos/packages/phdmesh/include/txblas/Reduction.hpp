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
 * @date   May 2007
 */

#ifndef txblas_Reduction_hpp
#define txblas_Reduction_hpp

#include <math.h>
#include <txblas/reduction.h>
#include <util/ParallelReduce.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** Summation of numerical values in double-double precision
 *  and mitigation of catastrophic cancelation for mixed
 *  positive and negative values.
 *
 *  Given double precision addition error bound of epsilon and
 *  the summation of N terms the error should be O(N*epsilon^2)
 *  as opposed to O(N*epsilon) for a conventional summation.
 *
 *  This should provide a permutation-insensitive result for
 *  the summation of a quadrillion double precision values
 *  N ~ O( 1,000,000,000,000,000 ), for epsilon < 1e-15.
 */

class Summation {
public:

  double xdval[4] ;

  inline
  ~Summation() {}

  inline
  Summation() { xdval[0] = 0 ; xdval[1] = 0 ; xdval[2] = 0 ; xdval[3] = 0 ; }

  inline
  Summation( const Summation & y )
    { xdval[0] = y.xdval[0] ; xdval[1] = y.xdval[1] ;
      xdval[2] = y.xdval[2] ; xdval[3] = y.xdval[3] ; }

  inline
  Summation & operator = ( const Summation & y )
    { xdval[0] = y.xdval[0] ; xdval[1] = y.xdval[1] ;
      xdval[2] = y.xdval[2] ; xdval[3] = y.xdval[3] ; return *this ; }

  inline
  Summation & operator += ( const double a )
    { xdsum_add_value( xdval , a ); return *this ; }

  inline
  Summation & operator += ( const Summation & a )
    { xdsum_add_dsum( xdval , a.xdval ); return *this ; }

  inline
  double value( double * remainder = NULL ) const
    {
      double x[2] ; xdsum_get_value( x , xdval );
      if ( remainder ) { *remainder = x[1] ; }
      return x[0] ;
    }
};

//----------------------------------------------------------------------

inline
double sum( unsigned n , const double * x )
{ Summation s ; txdsum_add_array( s.xdval , n , x ); return s.value(); }

inline
double dot( unsigned n , const double * x )
{ double d[2] ; txddot1( d , n , x ); return d[0] ; }

inline
double dot( unsigned n , const double * x , const double * y )
{ Summation s ; txddot( s.xdval , n , x , y ); return s.value(); }

inline
double norm2( unsigned n , const double * x )
{ return sqrt( dot( n , x ) ); }

inline
double norm1( unsigned n , const double * x )
{ double d[2] ; txdnorm1( d , n , x ); return d[0] ; }

//----------------------------------------------------------------------

inline
double sum( ParallelMachine comm , unsigned n , const double * x )
{
  Summation s ;
  txdsum_add_array( s.xdval , n , x );
  all_reduce( comm , Sum<1>( & s ) );
  return s.value();
}

inline
double dot( ParallelMachine comm , unsigned n , const double * x )
{
  Summation s ;
  txddot1( s.xdval , n , x );
  all_reduce( comm , Sum<1>( & s ) );
  return s.xdval[0] ;
}

inline
double dot( ParallelMachine comm , unsigned n , const double * x , const double * y )
{
  Summation s ;
  txddot( s.xdval , n , x , y );
  all_reduce( comm , Sum<1>( & s ) );
  return s.value();
}

inline
double norm2( ParallelMachine comm , unsigned n , const double * x )
{
  return sqrt( dot( comm , n , x ) );
}

inline
double norm1( ParallelMachine comm , unsigned n , const double * x )
{
  Summation s ;
  txdnorm1( s.xdval , n , x );
  all_reduce( comm , Sum<1>( & s ) );
  return s.xdval[0] ;
}

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

