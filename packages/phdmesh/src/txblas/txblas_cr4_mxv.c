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
#include <util/TPI.h>
#include <txblas/cr4_mxv.h>

/*--------------------------------------------------------------------*/

typedef struct txblasTask_cr_MatrixStruct {
  unsigned         number_row ;
  const unsigned * pc_begin ;
  const unsigned * ia_begin ;
  const double   * a_begin ;
  const double   * x_begin ;
        double   * y_begin ;
} txblasTask_cr_Matrix ;

/*--------------------------------------------------------------------*/

static void txblas_task_cr_mxv( void * data , TPI_ThreadPool pool )
{
  int p_size , p_rank ;

  if ( ! TPI_Rank( pool , & p_rank , & p_size ) ) {

    txblasTask_cr_Matrix * const t = (txblasTask_cr_Matrix*) data ;

    const unsigned beg_row = ( t->number_row * ( p_rank     ) ) / p_size ;
    const unsigned end_row = ( t->number_row * ( p_rank + 1 ) ) / p_size ;

    const unsigned * const pc_end = t->pc_begin + end_row ;
    const unsigned * const ia_beg = t->ia_begin ;
    const double   * const a_beg  = t->a_begin ;
    const double   * const x_beg  = t->x_begin ;
          double   *       y      = t->y_begin + beg_row ;

    const unsigned * pc = t->pc_begin + beg_row ;
    const unsigned * ia = ia_beg + *pc ;
    const double   * a  = a_beg  + *pc ;

    while ( pc < pc_end ) {
      double ytmp = 0 ;

      const unsigned * const ia_end = ia_beg + *++pc ;

      {
        enum { STRIDE = 4 };

        const unsigned * const ia_blk = ia_end - ( ia_end - ia ) % STRIDE ;

        for ( ; ia < ia_blk ; ia += STRIDE , a += STRIDE ) {
          ytmp += a[0] * x_beg[ ia[0] ] +
                  a[1] * x_beg[ ia[1] ] +
                  a[2] * x_beg[ ia[2] ] +
                  a[3] * x_beg[ ia[3] ] ;
        }
      }

      for ( ; ia < ia_end ; ++ia , ++a ) {
        ytmp += *a * x_beg[ *ia ];
      }

      *y++ = ytmp ;
    }
  }
}


void txblas_cr_mxv(
  const unsigned nr  /* Number rows */ ,
  const unsigned pc[] ,
  const unsigned ia[] ,
  const double   a[] ,
  const double   x[] ,  /* Input vector */
        double   y[] )  /* Output vector */
{
  txblasTask_cr_Matrix data = { nr , pc , ia , a , x , y };
  TPI_Run( & txblas_task_cr_mxv , & data , 0 );
}

/*--------------------------------------------------------------------*/
