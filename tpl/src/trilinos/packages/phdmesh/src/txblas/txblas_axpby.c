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
#include "stdlib.h"

struct TaskXY {
  const double   alpha ;
  const double   beta ;
  const double * x_beg ;
        double * y_beg ;
        unsigned number ;
        unsigned * iter ;
};

static void task_axpby_work( void * , TPI_ThreadPool );
static void task_axpby_work_block( void * , TPI_ThreadPool );
static void task_axpby_work_steal( void * , TPI_ThreadPool );

void tdaxpby( unsigned n ,
              double a , const double * x ,
              double b , double * y ,
              int block )
{
  int p_size ;
  TPI_Size( & p_size );

  {
    unsigned *tmp = malloc( p_size );
    struct TaskXY data = { a , b , x , y , n , tmp };
    int i ;
    for ( i = 0 ; i < p_size ; ++i ) { tmp[i] = i ; }
    if ( 0 < block ) {
      TPI_Run( & task_axpby_work_block , & data , 0 );
    }
    else if ( block < 0 ) {
      TPI_Set_lock_size( p_size );
      TPI_Run( & task_axpby_work_steal , & data , 0 );
    }
    else {
      TPI_Run( & task_axpby_work , & data , 0 );
    }
    free(tmp);
  }
}

/*--------------------------------------------------------------------*/

enum { UNROLL = 8 };

static void daxpby_work( unsigned n ,
                         const double a , const double * x ,
                         const double b , double * y )
{
  const double * const xe = x + n ;
  const double * const xb = xe - n % UNROLL ;

  for ( ; x < xb ; x += UNROLL , y += UNROLL ) {
    y[0] = a * x[0] + b * y[0] ;
    y[1] = a * x[1] + b * y[1] ;
    y[2] = a * x[2] + b * y[2] ;
    y[3] = a * x[3] + b * y[3] ;
    y[4] = a * x[4] + b * y[4] ;
    y[5] = a * x[5] + b * y[5] ;
    y[6] = a * x[6] + b * y[6] ;
    y[7] = a * x[7] + b * y[7] ;
  }

  for ( ; x < xe ; ++x , ++y ) {
    *y = a * *x + b * *y ;
  }
}

static void task_axpby_work( void * arg , TPI_ThreadPool pool )
{
  int p_size ;
  int p_rank ;

  if ( ! TPI_Rank( pool , & p_rank , & p_size ) ) {

    struct TaskXY * const t = (struct TaskXY *) arg ;

    const int      n_rem = t->number % p_size ;
    const unsigned n_num = t->number / p_size ;
    const unsigned n_beg = p_rank * n_num + ( p_rank < n_rem ? p_rank : n_rem );
    const unsigned n_len = n_num + ( p_rank < n_rem ? 1 : 0 );

    daxpby_work( n_len , t->alpha , t->x_beg + n_beg ,
                         t->beta  , t->y_beg + n_beg );
  }
}

static void task_axpby_work_block( void * arg , TPI_ThreadPool pool )
{
  enum { BLOCK = UNROLL * 1024 };
  int p_size ;
  int p_rank ;

  if ( ! TPI_Rank( pool , & p_rank , & p_size ) ) {

    struct TaskXY * const t = (struct TaskXY *) arg ;

    const unsigned inc = BLOCK * p_size ;
    const unsigned num = t->number ;

    if ( 1 < p_size && 2 * inc < num ) { /* More than two blocks of work */
      const double a = t->alpha ;
      const double b = t->beta ;
      const double * const x = t->x_beg ;
            double * const y = t->y_beg ;

      int len ;

      unsigned i ;
      for ( i = BLOCK * p_rank ; 0 < ( len = num - i ) ; i += inc ) {
        daxpby_work( ( BLOCK < len ? BLOCK : len ) , a , x + i , b , y + i );
      }
    }
    else { /* Even partitioning */
      const int      n_rem = num % p_size ;
      const unsigned n_num = num / p_size ;
      const unsigned n_beg = p_rank*n_num + (p_rank < n_rem ? p_rank : n_rem);
      const unsigned n_len = n_num + ( p_rank < n_rem ? 1 : 0 );

      daxpby_work( n_len , t->alpha , t->x_beg + n_beg ,
                           t->beta ,  t->y_beg + n_beg );
    }
  }
}

static void task_axpby_work_steal( void * arg , TPI_ThreadPool pool )
{
  enum { BLOCK = UNROLL * 128 };
  int p_size ;
  int p_rank ;

  if ( ! TPI_Rank( pool , & p_rank , & p_size ) ) {

    struct TaskXY * const t = (struct TaskXY *) arg ;

    const double   a = t->alpha ;
    const double   b = t->beta ;
    const unsigned n = t->number ;
    const double * const x = t->x_beg ;
          double * const y = t->y_beg ;

    unsigned * const all_iter = t->iter ;
    unsigned * const my_iter  = all_iter + p_size ;

    {
      unsigned i ;
      for ( i = 0 ; i < n ; ) {
        TPI_Lock( pool , p_rank );
        i = *my_iter * BLOCK ; *my_iter += p_size ;
        TPI_Unlock( pool , p_rank );
        if ( i < n ) {
          const unsigned len = BLOCK < n - i ? BLOCK : n - i ;
          daxpby_work( len, a, x + i, b, y + i );
        }
      }
    }

    /* Finished my work, steal work from someone else */

    {
    int working ; 
    int p = 0 ;
    for ( working = 1 ; working ; ) {
      working = 0 ;
      for ( p = 0 ; p < p_size ; ++p ) {
        if ( all_iter[p] * BLOCK < n ) {
          if ( ! TPI_Trylock( pool , p ) ) {
            const unsigned i = all_iter[p] * BLOCK ;
            all_iter[p] += p_size ;
            TPI_Unlock( pool , p );
            if ( i < n ) {
              const unsigned len = BLOCK < n - i ? BLOCK : n - i ;
              daxpby_work( len, a, x + i, b, y + i );
            }
          }
          working = 1 ;
        }
      }
    }
    }
  }
}


