/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define max(x,y) (( x > y ) ? x : y)     /* max function  */
double smsrres (int m, int n, 
	      double *val, int *indx, 
	      double *xlocal, double *x, double *b)
{
    int i, j, jbgn, jend, ione = 1;
    double sum, norm_tmp = 0.0, norm_b = 0.0;
    double scaled_res_norm, res_norm, *tmp, max_norm = 0.0;


/*     Computes the residual

                      res = || b - A*x ||

       where x and b are vectors and A is a sparse matrix stored
       in MSR format. */

/*     -------------------------- 
       First executable statement 
       -------------------------- */

    /* Create tmp workspace */
    tmp = (double *) calloc(m,sizeof(double));

/* .....initialize soln */

    for (i = 0; i < m; i++)
	tmp[i] = b[i] - val[i] * xlocal[i];

/* .....do a series of SPDOTs (sparse dot products) */

    for (i = 0; i <m ; i++) 
      {
	jbgn = indx[i];
	jend = indx[i + 1];
	sum = 0.0;
	
	for (j = jbgn; j < jend; j++)
	  sum += val[j] * x[indx[j]];

	tmp[i] -= sum;
	max_norm = max(fabs(tmp[i]),max_norm);
	norm_tmp += tmp[i]*tmp[i];
	norm_b += b[i]*b[i];
      }
   
    res_norm = sqrt(norm_tmp);
    printf("\n\nMax norm of residual        = %12.4g\n",max_norm);
    printf(    "Two norm of residual        = %12.4g\n",res_norm);
    if (norm_b > 1.0E-7) 
      {
	   scaled_res_norm = res_norm/sqrt(norm_b);
	   printf(    "Scaled two norm of residual = %12.4g\n",scaled_res_norm);
      }
    free((void *) tmp);

    return(scaled_res_norm);

} /* smsrres */

