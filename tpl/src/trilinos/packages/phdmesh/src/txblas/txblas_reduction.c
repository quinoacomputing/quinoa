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

#include <stddef.h>
#include <math.h>
#include <util/TPI.h>
#include <txblas/reduction.h>
#include <stdlib.h>

#define BLOCKING_SIZE 0

/*--------------------------------------------------------------------*/


#define SUM_ADD( v , a ) \
  { \
    const double VpA = v[0] + a ; \
    v[1] += a < v[0] ? ( a - ( VpA - v[0] ) ) : ( v[0] - ( VpA - a ) ); \
    v[0]  = VpA + v[1] ; \
    v[1] += VpA - v[0] ; \
  }

#define SUM_SUB( v , a ) \
  { \
    const double VmA = v[0] - a ; \
    v[1] += a < fabs(v[0]) ? ( -a - ( VmA - v[0] ) ) : ( v[0] - ( VmA + a ) ); \
    v[0]  = VmA + v[1] ; \
    v[1] += VmA - v[0] ; \
  }

/*--------------------------------------------------------------------*/

void xdsum_add_value( double * v , double a )
{
  if ( a < 0 ) { a = -a ; v += 2 ; }
  SUM_ADD( v , a );
}

void xdsum_add_dsum( double * v , const double * const a )
{
  SUM_ADD( v , a[0] );
  SUM_ADD( v , a[1] );
  v += 2 ;
  SUM_ADD( v , a[2] );
  SUM_ADD( v , a[3] );
}

void xdsum_get_value( double * const y , const double * const v )
{
  y[0] = v[0] ;
  y[1] = v[1] ;
  SUM_SUB( y , v[2] );
  SUM_SUB( y , v[3] );
}

/*--------------------------------------------------------------------*/

struct TaskX {
        double * x_sum ;
  const double * x_beg ;
        unsigned number ;
};

struct TaskXY {
        double * xy_sum ;
  const double * x_beg ;
  const double * y_beg ;
        unsigned number ;
        unsigned block ;
};

/*--------------------------------------------------------------------*/

static
void add_array( double * const s , const double * x , const double * const xe )
{
  while ( xe != x ) {
    double * p;
    double   a = *x ; ++x ;
    p = s ;
    if ( a < 0 ) { a = -a ; p += 2 ; }
    SUM_ADD( p , a );
  }
}

void xdsum_add_array( double * s , unsigned n , const double * x )
{ add_array( s , x , x + n ); }

static void task_sum_work( void * arg , TPI_ThreadPool pool )
{
  int p_size , p_rank ;

  if ( ! TPI_Rank( pool , & p_rank , & p_size ) ) {

    struct TaskX * const t  = (struct TaskX *) arg ;

    const unsigned p_next = p_rank + 1 ;
    const unsigned n = t->number ;
    const double * const xb = t->x_beg + ( n * p_rank ) / p_size ;
    const double * const xe = t->x_beg + ( n * p_next ) / p_size ;
          double * const v  = t->x_sum ;

    double partial[4] = { 0 , 0 , 0 , 0 };

    add_array( partial , xb , xe );

    TPI_Lock( pool , 0 );

    xdsum_add_dsum( v , partial );

    TPI_Unlock( pool , 0 );
  }
}

void txdsum_add_array( double * s , unsigned n , const double * x )
{
  struct TaskX data ;
  data.x_sum  = s ;
  data.x_beg  = x ;
  data.number = n ;
  TPI_Set_lock_size( 1 );
  TPI_Run( & task_sum_work , & data , 0 );
}

/*--------------------------------------------------------------------*/

static void norm1( double * s , const double * x , const double * const xe )
{
  while ( xe != x ) { const double a = fabs(*x); ++x ; SUM_ADD( s , a ); }
}

void xdnorm1( double * s2 , unsigned n , const double * x )
{ norm1( s2 , x , x + n ); }

static void task_norm1_work( void * arg , TPI_ThreadPool pool )
{
  int p_size , p_rank ;

  if ( ! TPI_Rank( pool , & p_rank , & p_size ) ) {

    struct TaskX * const t  = (struct TaskX *) arg ;

    const unsigned p_next = p_rank + 1 ;
    const unsigned n = t->number ;
    const double * const xb = t->x_beg + ( n * p_rank ) / p_size ;
    const double * const xe = t->x_beg + ( n * p_next ) / p_size ;
          double * const v  = t->x_sum ;

    double partial[2] = { 0 , 0 };

    norm1( partial , xb , xe );

    TPI_Lock( pool , 0 );

    SUM_ADD( v , partial[0] );
    SUM_ADD( v , partial[1] );

    TPI_Unlock( pool , 0 );
  }
}

void txdnorm1( double * s , unsigned n , const double * x )
{
  struct TaskX data ;
  data.x_sum  = s ;
  data.x_beg  = x ;
  data.number = n ;
  TPI_Set_lock_size( 1 );
  TPI_Run( & task_norm1_work , & data , 0 );
}

/*--------------------------------------------------------------------*/

static void dot1_unroll( double * s , const double * x , const size_t n )
{
  enum { NB = 8 };

  const double * const x_blk = x + n % NB ;
  for ( ; x_blk != x ; ++x ) {
    double a = *x ; a *= a ; SUM_ADD( s , a );
  }

  {const double * const x_end = x + n ;
  for ( ; x_end != x ; x += NB ) {
    double a0 = x[0] ;
    double a1 = x[1] ;
    double a2 = x[2] ;
    double a3 = x[3] ;
    double a4 = x[4] ;
    double a5 = x[5] ;
    double a6 = x[6] ;
    double a7 = x[7] ;
    a0 *= a0 ;
    a1 *= a1 ;
    a2 *= a2 ;
    a3 *= a3 ;
    a4 *= a4 ;
    a5 *= a5 ;
    a6 *= a6 ;
    a7 *= a7 ;
    SUM_ADD( s , a0 );
    SUM_ADD( s , a1 );
    SUM_ADD( s , a2 );
    SUM_ADD( s , a3 );
    SUM_ADD( s , a4 );
    SUM_ADD( s , a5 );
    SUM_ADD( s , a6 );
    SUM_ADD( s , a7 );
  }
  }
}

void xddot1( double * s2 , unsigned n , const double * x )
{ dot1_unroll( s2 , x , n ); }

static void task_xddot_x_work( void * arg , TPI_ThreadPool pool )
{
  int p_size , p_rank ;

  if ( ! TPI_Rank( pool , & p_rank , & p_size ) ) {

    double partial[2] = { 0 , 0 };
    struct TaskX * const t  = (struct TaskX *) arg ;

    {
      const unsigned p_next   = p_rank + 1 ;
      const unsigned n_global = t->number ;
      const unsigned n_begin  = ( ( n_global * p_rank ) / p_size );
      const unsigned n_local  = ( ( n_global * p_next ) / p_size ) - n_begin ;

      dot1_unroll( partial , t->x_beg + n_begin , n_local );
    }

    {
      TPI_Lock(pool,0);
      {
      double * const v = t->x_sum ;
      SUM_ADD( v , partial[0] );
      SUM_ADD( v , partial[1] );
      TPI_Unlock(pool,0);
      }
    }
  }
}

void txddot1( double * s , unsigned n , const double * x )
{
  struct TaskX data ;
  data.x_sum  = s ;
  data.x_beg  = x ;
  data.number = n ;
  TPI_Set_lock_size( 1 );
  TPI_Run( & task_xddot_x_work , & data , 0 );
}

/*--------------------------------------------------------------------*/

double ddot( unsigned n , const double * x , const double * y )
{
  enum { STRIDE = 8 };

  const double * const x_end = x + n ;
  const double * const x_blk = x_end - n % STRIDE ;

  double result = 0 ;

  for ( ; x < x_blk ; x += STRIDE , y += STRIDE ) {
    result += x[0] * y[0] +
              x[1] * y[1] +
              x[2] * y[2] +
              x[3] * y[3] +
              x[4] * y[4] +
              x[5] * y[5] +
              x[6] * y[6] +
              x[7] * y[7] ;
  }

  for ( ; x < x_end ; ++x , ++y ) {
    result += *x * *y ;
  }

  return result ;
}

static void task_ddot_xy_work_blocking( void * arg , TPI_ThreadPool pool )
{
  int p_size , p_rank ;

  if ( ! TPI_Rank( pool , & p_rank , & p_size ) ) {

    struct TaskXY * const t = (struct TaskXY *) arg ;

    const unsigned block_size   = t->block ;
    const unsigned block_start  = block_size * p_rank ;
    const unsigned block_stride = block_size * p_size ;

    const double * const x_end = t->x_beg + t->number ;
    const double * x = t->x_beg + block_start ;
    const double * y = t->x_beg + block_start ;

    double local = 0.0 ;

    for ( ; x < x_end ; x += block_stride , y += block_stride ) {
      const unsigned n = x_end - x ;
      local += ddot( ( block_size < n ? block_size : n ) , x , y );
    }

    t->xy_sum[ p_rank ] = local ;
  }
}

static void task_ddot_xy_work( void * arg , TPI_ThreadPool pool )
{
  int p_size , p_rank ;

  if ( ! TPI_Rank( pool , & p_rank , & p_size ) ) {

    struct TaskXY * const t = (struct TaskXY *) arg ;

    const unsigned n_total = t->number ;
    const unsigned n_begin = ( n_total * ( p_rank     ) ) / p_size ;
    const unsigned n_end   = ( n_total * ( p_rank + 1 ) ) / p_size ;
    const unsigned n_local = ( n_end - n_begin );

    const double * x = t->x_beg + n_begin ;
    const double * y = t->y_beg + n_begin ;

    t->xy_sum[ p_rank ] = ddot( n_local , x , y );
  }
}

void tddot( double * s , unsigned n , const double * x , const double * y )
{
  int p_size ;
  if ( ! TPI_Size( & p_size ) ) {
    double* tmp = malloc( p_size * sizeof(double));
    struct TaskXY data = { tmp , x , y , n , BLOCKING_SIZE };
    int i ;
    for ( i = 0 ; i < p_size ; ++i ) { tmp[i] = 0 ; }
    if ( data.block ) {
      TPI_Run( & task_ddot_xy_work_blocking , & data , 0 );
    }
    else {
      TPI_Run( & task_ddot_xy_work , & data , 0 );
    }
    for ( i = 1 ; i < p_size ; ++i ) { tmp[0] += tmp[i] ; }
    *s = tmp[0] ;
    free(tmp);
  }
}

/*--------------------------------------------------------------------*/

void xddot( double * s4 , unsigned n , const double * x , const double * y )
{
  enum { STRIDE = 8 };

  const double * const x_end = x + n ;
  const double * const x_blk = x_end - n % STRIDE ;

  for ( ; x < x_blk ; x += STRIDE , y += STRIDE ) {
    xdsum_add_value( s4 , x[0] * y[0] );
    xdsum_add_value( s4 , x[1] * y[1] );
    xdsum_add_value( s4 , x[2] * y[2] );
    xdsum_add_value( s4 , x[3] * y[3] );
    xdsum_add_value( s4 , x[4] * y[4] );
    xdsum_add_value( s4 , x[5] * y[5] );
    xdsum_add_value( s4 , x[6] * y[6] );
    xdsum_add_value( s4 , x[7] * y[7] );
  }

  for ( ; x < x_end ; ++x , ++y ) {
    xdsum_add_value( s4 , *x * *y );
  }
}

/*--------------------------------------------------------------------*/

static void task_xddot_xy_work_blocking( void * arg , TPI_ThreadPool pool )
{
  int p_size , p_rank ;

  if ( ! TPI_Rank( pool , & p_rank , & p_size ) ) {

    struct TaskXY * const t = (struct TaskXY *) arg ;

    const unsigned block_size   = t->block ;
    const unsigned block_start  = block_size * p_rank ;
    const unsigned block_stride = block_size * p_size ;

    const double * const x_end = t->x_beg + t->number ;
    const double * x = t->x_beg + block_start ;
    const double * y = t->x_beg + block_start ;

    double s_local[4] = { 0 , 0 , 0 , 0 };

    for ( ; x < x_end ; x += block_stride , y += block_stride ) {
      const unsigned n = x_end - x ;
      xddot( s_local , ( block_size < n ? block_size : n ) , x , y );
    }

    {
    double * const xy_sum = t->xy_sum + 4 * p_rank ;

    xy_sum[0] = s_local[0] ;
    xy_sum[1] = s_local[1] ;
    xy_sum[2] = s_local[2] ;
    xy_sum[3] = s_local[3] ;
    }
  }
}

static void task_xddot_xy_work( void * arg , TPI_ThreadPool pool )
{
  int p_size , p_rank ;

  if ( ! TPI_Rank( pool , & p_rank , & p_size ) ) {

    struct TaskXY * const t = (struct TaskXY *) arg ;

    const unsigned n_total = t->number ;
    const unsigned n_begin = ( n_total * ( p_rank     ) ) / p_size ;
    const unsigned n_end   = ( n_total * ( p_rank + 1 ) ) / p_size ;
    const unsigned n_local = ( n_end - n_begin );

    const double * const x = t->x_beg + n_begin ;
    const double * const y = t->y_beg + n_begin ;

    double s_local[4] = { 0 , 0 , 0 , 0 };

    xddot( s_local , n_local , x , y );

    {
    double * const xy_sum = t->xy_sum + 4 * p_rank ;

    xy_sum[0] = s_local[0] ;
    xy_sum[1] = s_local[1] ;
    xy_sum[2] = s_local[2] ;
    xy_sum[3] = s_local[3] ;
    }
  }
}


void txddot( double * s , unsigned n , const double * x , const double * y )
{
  int p_size ;
  if ( ! TPI_Size( & p_size ) ) {
    double* tmp;
    const int ntmp = 4 * p_size ;
    tmp = malloc(ntmp * sizeof(double));
    {struct TaskXY data = { tmp , x , y , n , BLOCKING_SIZE };
    int i ;
    for ( i = 0 ; i < ntmp ; ++i ) { tmp[i] = 0 ; }
    if ( data.block ) {
      TPI_Run( & task_xddot_xy_work_blocking , & data , 0 );
    }
    else {
      TPI_Run( & task_xddot_xy_work , & data , 0 );
    }
    for ( i = 0 ; i < p_size ; ++i ) {
      xdsum_add_dsum( s , tmp + 4 * i );
    }
    }
    free(tmp);
  }
}

/*--------------------------------------------------------------------*/

