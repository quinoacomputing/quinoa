//  ------------------------------------------------------------------------------------------------------------
//
//  Copyright 2007 Jozsef Bakosi
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  ------------------------------------------------------------------------------------------------------------
//
//  Functions that deal with 3x3 matrices
//  for more info see main.cc
//

#include <stdio.h>
#include <math.h>
#include "matrix3.h"


static double m2[9], m3[9], m4[9], m5[9];






double m3det( double *m )
//
// returns the determinant of 3x3 matrix m
//
{
  return( m[0]*(m[4]*m[8]-m[5]*m[7]) + m[1]*(m[6]*m[5]-m[3]*m[8]) + m[2]*(m[3]*m[7]-m[6]*m[4]) );
}




void m3inv( double *m, double *invm )
//
// calculates the inverse matrix (invm) for a 3x3 matrix (m)
//
{
  double det;

  
  det = m3det(m);

  invm[0] = (m[4]*m[8]-m[7]*m[5])/det;
  invm[1] = (m[2]*m[7]-m[8]*m[1])/det;
  invm[2] = (m[1]*m[5]-m[4]*m[2])/det;
  invm[3] = (m[5]*m[6]-m[8]*m[3])/det;
  invm[4] = (m[0]*m[8]-m[6]*m[2])/det;
  invm[5] = (m[2]*m[3]-m[5]*m[0])/det;
  invm[6] = (m[3]*m[7]-m[6]*m[4])/det;
  invm[7] = (m[1]*m[6]-m[7]*m[0])/det;
  invm[8] = (m[0]*m[4]-m[3]*m[1])/det;
}





void m3mult( double *A, double *B, double *mult )
//
// calculates the product of two 3x3 matrices,
// mult = A*B
//
{
  mult[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
  mult[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
  mult[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8];
  mult[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
  mult[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
  mult[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8];
  mult[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6];
  mult[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7];
  mult[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
}






void m3exp( double *m, double *expm )
//
// calculates the matrix exponential for a 3x3 matrix
// using the definition of the matrix exponential (power series)
//
{
  m3mult( m, m, m2 );	// m2 = m*m
  m3mult( m2, m, m3 );	// m3 = m2*m
  m3mult( m3, m, m4 );	// m4 = m3*m
  m3mult( m4, m, m5 );	// m5 = m4*m

  expm[0] = 1.0 + m[0] + 0.5*m2[0] + m3[0]/6.0 + m4[0]/24.0 + m5[0]/120.0;
  expm[1] = m[1] + 0.5*m2[1] + m3[1]/6.0 + m4[1]/24.0 + m5[1]/120.0;
  expm[2] = m[2] + 0.5*m2[2] + m3[2]/6.0 + m4[2]/24.0 + m5[2]/120.0;
  expm[3] = m[3] + 0.5*m2[3] + m3[3]/6.0 + m4[3]/24.0 + m5[3]/120.0;
  expm[4] = 1.0 + m[4] + 0.5*m2[4] + m3[4]/6.0 + m4[4]/24.0 + m5[4]/120.0;
  expm[5] = m[5] + 0.5*m2[5] + m3[5]/6.0 + m4[5]/24.0 + m5[5]/120.0;
  expm[6] = m[6] + 0.5*m2[6] + m3[6]/6.0 + m4[6]/24.0 + m5[6]/120.0;
  expm[7] = m[7] + 0.5*m2[7] + m3[7]/6.0 + m4[7]/24.0 + m5[7]/120.0;
  expm[8] = 1.0 + m[8] + 0.5*m2[8] + m3[8]/6.0 + m4[8]/24.0 + m5[8]/120.0;
}





void m3trans( double *m, double *trans )
//
// calculates the sum of a matrix and its transpose,
// trans = m+m^T for a 3x3 matrix
//
{
  trans[0] = 2.0*m[0]; 
  trans[1] = trans[3] = m[1] + m[3];
  trans[2] = trans[6] = m[2] + m[6];
  trans[4] = 2.0*m[4];
  trans[5] = trans[7] = m[5] + m[7];
  trans[8] = 2.0*m[8];
}





void m3cholesky( double *m )
//
// overwrites the 3x3 matrix m with its Cholesky
// decomposition's lower triangular part (including the diagonal,
// also puts zeros in the upper diagonal)
//
{
  m[1] = m[2] = m[5] = 0.0;
  m[0] = sqrt(m[0]);
  m[3] /= m[0];
  m[6] /= m[0];
  m[4] = sqrt(m[4]-m[3]*m[3]);
  m[7] = (m[7]-m[6]*m[3])/m[4];
  m[8] = sqrt(m[8]-m[6]*m[6]-m[7]*m[7]);
}




void m3out( double *m, char *name, FILE *ofile )
//
// outputs 3x3 matrix m to file 'ofile' with name 'name'
// ofile should be open
//
{
  fprintf( ofile, "%s = [%g, %g, %g; %g, %g, %g; %g, %g, %g]\n", name,
		  m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8] );
}
