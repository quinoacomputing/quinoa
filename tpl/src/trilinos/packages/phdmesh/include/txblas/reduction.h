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
 * @date   July 2007
 */

#ifndef txblas_reduction_h
#define txblas_reduction_h

#if defined( __cplusplus )
extern "C" {
#endif

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
 *
 *  Each 'value' consists of a large magnitude value and remainder
 *  such that:  large == large + remainder.  Summations consist of
 *  a positive value-pair and negative value-pair.  The negative
 *  value pair is accumulated as a positive value.
 *
 *  Example usage for a permutation-insensitive result:
 *    double dot( unsigned n , const double * x , const double * y )
 *    {
 *      double sum[4] ;
 *      phdmesh_dot( sum , n , x , y );
 *      phdmesh_dsum_value( sum , sum );
 *      return sum[0] ;
 *    }
 */

/** Add-in a double value.
 *    s4[0] = positive large magnitude value
 *    s4[1] = positive remainder
 *    s4[2] = negative large magnitude value, negated
 *    s4[3] = negative remainder, negated
 */
void xdsum_add_value( double * s4 , double );

/** Add-in a summation value, equivalent to:
 *    phdmesh_dsum_add( s4 ,   a4[0] );
 *    phdmesh_dsum_add( s4 ,   a4[1] );
 *    phdmesh_dsum_add( s4 , - a4[2] );
 *    phdmesh_dsum_add( s4 , - a4[3] );
 */
void xdsum_add_dsum( double * s4 , const double * a4 );

/** Extract the value and remainder.
 *    v2[0] = s4'positive - s4'negative, large magnitude
 *    v2[1] = s4'positive - s4'negative, remainder
 */
void xdsum_get_value( double * v2 , const double * s4 );

/** Summation of a set of numbers, equivalent to:
 *    for ( i = 0 ; i < n ; ++i ) {
 *      xdsum_add_value( s4 , x[i] );
 *    }
 */
void xdsum_add_array( double * s4 , unsigned n , const double * x );

/** Task-parallel summation of a set of numbers. */
void txdsum_add_array( double * s4 , unsigned n , const double * x );

void xdnorm1(  double * s2 , unsigned , const double * );
void txdnorm1( double * s2 , unsigned , const double * );

/** Dot product with itself, equivalent to:
 *    for ( i = 0 ; i < n ; ++i ) {
 *      phdmesh_dsum_sum( s2 , x[i] * x[i] );
 *    }
 *  Note that the contributions are always positive,
 *  thus the negative part is untouched.
 */
void xddot1(  double * s2 , unsigned , const double * );
void txddot1( double * s2 , unsigned , const double * );

/** Dot product, equivalent to:
 *    for ( i = 0 ; i < n ; ++i ) {
 *      phdmesh_dsum_sum( s4 , x[i] * y[i] );
 *    }
 */
void xddot(   double * s4 , unsigned , const double * , const double * );

void txddot( double * s4 , unsigned , const double * , const double * );


double ddot( unsigned , const double * , const double * );

void tddot( double * s , unsigned , const double * , const double * );

#if defined( __cplusplus )
} /* extern "C" */
#endif

#endif

