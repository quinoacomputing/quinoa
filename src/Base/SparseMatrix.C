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
//  Functions dealing with sparse matrices. For mor info, see main.cc.
//
//


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include "omp.h"
#endif
#include "macros.h"
#include "sparsemat.h"









void create_sparsemat( sparsemat *m, int DOF, int size, int *psup1, int *psup2 )
//
// allocates data structures for a size x size sparse symmetric matrix with DOF degrees of freedom,
// ie. the real size will be (size x DOF) x (size x DOF) and symmetric, storing only the upper triangle
//
// see also header-file for structure definition
//
// also allocates vectors for the conjugate gradients solve: r, p, z, q, d, u
// also allocates linked lists for level scheduling for SGS preconditioner
//
{
  int i, j, k, l, n, e, itmp;
  int *depth;


  // put in size
  m->size = size;

  // put in real size
  m->rsize = size*DOF;

  // put in number of unknowns/point
  m->dof = DOF;

  // allocate vectors for conjugate gradients solve
  if ( !(m->r = (double*)calloc(m->rsize,sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(m->p = (double*)calloc(m->rsize,sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(m->z = (double*)calloc(m->rsize,sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(m->q = (double*)calloc(m->rsize,sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(m->d = (double*)calloc(m->rsize,sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(m->u = (double*)calloc(m->rsize,sizeof(double))) ) ERR("Can't allocate memory!");

  // allocate array for row indices
  if ( !(m->ia = (int*)calloc(size*DOF+1,sizeof(int))) ) ERR("Can't allocate memory!");
  // allocate array for for storing the nonzeros in each row
  if ( !(m->rnz = (int*)calloc(size,sizeof(int))) ) ERR("Can't allocate memory!");

  // calculate number of nonzeros in each block row (rnz[]),
  // total number of nonzeros (nnz) and fill row indices (ia[])
  for ( m->ia[0]=1, m->nnz=i=0; i < size; i++ )
  {
    // add up and store nonzeros of row i (only upper triangular part, matrix is symmetric)
    for ( m->rnz[i] = 1, j = psup2[i]+1; j <= psup2[i+1]; j++ )
      m->rnz[i]++;

    // add up total number of nonzeros
    m->nnz += m->rnz[i]*DOF;
    
    // fill up rowindex
    for ( k = 0; k < DOF; k++ )
      m->ia[i*DOF+k+1] = m->ia[i*DOF+k] + m->rnz[i];
  }

  // allocate array for nonzero matrix values
  if ( !(m->a = (double*)calloc(m->nnz,sizeof(double))) ) ERR("Can't allocate memory!");
  // allocate array for column indices
  if ( !(m->ja = (int*)calloc(m->nnz,sizeof(int))) ) ERR("Can't allocate memory!");

  // allocate temporary array to store depths for each row
  if ( !(depth = (int*)calloc(m->rsize,sizeof(int))) ) ERR("Can't allocate memory!");
  
  // fill column indices and compute depths of rows (for level scheduling)
  m->nlev = 0;
  for ( i = 0; i < size; i++ )					// loop over all points
    for ( k = 0; k < DOF; k++ )					// loop over all degrees of freedom in a point
    {
      itmp = i*DOF+k;
      m->ja[m->ia[itmp]-1] = itmp+1;				// put in column index of main diagonal
      for ( l=0, n=1, j=psup2[i]+1; j <= psup2[i+1]; j++ )	// loop over all points surrounding point i
      {
	m->ja[m->ia[itmp]-1+(n++)] = e = psup1[j]*DOF+k+1;	// put in column index of an off-diagonal
        if ( (i > psup1[j]) && (depth[e-1] > l) )		// find maximum depth of previous rows needed
	  l = depth[e-1];					// to compute row i (consider lower triangle)
      }
      depth[itmp] = 1 + l;					// store depth of row i
      if ( depth[itmp] > m->nlev ) m->nlev = depth[itmp];	// find maximum level
    }

  // (bubble-)sort column indices
  for ( i = 0; i < size; i++ )					// loop over all points
    for ( k = 0; k < DOF; k++ )					// loop over all degrees of freedom in a point
      for ( j = psup2[i]+1; j <= psup2[i+1]; j++ )		// loop over all points surrounding point i
         for ( l = 1; l < m->rnz[i]; l++ )			// sort column indices of row i
            for ( e = 0; e < (m->rnz[i]-l); e++ )
              if ( m->ja[m->ia[i*DOF+k]-1+e] > m->ja[m->ia[i*DOF+k]+e] )
	        SWAP(m->ja[m->ia[i*DOF+k]-1+e], m->ja[m->ia[i*DOF+k]+e], itmp);

  // allocate arrays for linked lists: lev1, lev2
  // lev1: permutation of rows according to consecutive levels
  // lev2: pointers to the beginning of the ith level in lev1
  if ( !(m->lev1 = (int*)calloc(m->rsize,sizeof(int))) ) ERR("Can't allocate memory!");
  if ( !(m->lev2 = (int*)calloc(m->nlev+1,sizeof(int))) ) ERR("Can't allocate memory!");

  // generate linked lists: lev1, lev2 (for level scheduling)
  m->lev2[0] = 0;
  for ( n=j=0; j < m->nlev; j++ )				// loop over all depths
  {
    for ( i = 0; i < m->rsize; i++ )				// loop over all rows
      if ( depth[i] == j+1 ) m->lev1[n++] = i;			// collect row indices for depth j+1
    m->lev2[j+1] = n;						// store pointer to beginning of depth j+1
  }

  // free temporary array
  free( depth );
}







void destroy_sparsemat( sparsemat *m )
//
// frees memory allocated for data structures of sparse matrix m
//
{
  free( m->rnz );
  free( m->ia );
  free( m->ja );
  free( m->a );
  free( m->r );
  free( m->p );
  free( m->z );
  free( m->q );
  free( m->d );
  free( m->u );
  free( m->lev1 );
  free( m->lev2 );
}










void addmr( sparsemat *M, int row, int column, int i, double value )
//
// adds a value to sparse matrix M into a specified position
// using relative addressing
//
// row :    block row
// column : block column
// i :      position in block
// value :  value to add
//
{
  int n, j, idx, rmdof;


  rmdof = row * M->dof;

  for (n=0, j=M->ia[rmdof+i]-1; j <= M->ia[rmdof+i+1]-2; j++, n++)
    if (column*M->dof+i+1 == M->ja[j])
    {
      idx = n;
      j = M->ia[rmdof+i+1]-1;
    }

  #ifdef _OPENMP
  #pragma omp atomic
  #endif
  M->a[M->ia[rmdof+i]-1+idx] += value;
}






void insmr( sparsemat *M, int row, int column, int i, double value )
//
// inserts a value to sparse matrix M into a specified position
// using relative addressing
//
// row :    block row
// column : block column
// i :      position in block
// value :  value to insert
//
{
  int n, j, idx, rmdof;


  rmdof = row * M->dof;
  
  for (n=0, j=M->ia[rmdof+i]-1; j <= M->ia[rmdof+i+1]-2; j++, n++)
    if (column*M->dof+i+1 == M->ja[j])
    {
      idx = n;
      j = M->ia[rmdof+i+1]-1;
    }

  M->a[M->ia[rmdof+i]-1+idx] = value;
}







void addma( sparsemat *M, int row, int column, double value )
//
// adds a value to sparse matrix M into a specified position
// using absolute addressing
//
// row :    block row
// column : block column
// value :  value to add
//
{
  int n, j, idx;


  for (n=0, j=M->ia[row]-1; j <= M->ia[row+1]-2; j++, n++)
    if (column+1 == M->ja[j])
    {
      idx = n;
      j = M->ia[row+1]-1;
    }

  #ifdef _OPENMP
  #pragma omp atomic
  #endif
  M->a[M->ia[row]-1+idx] += value;
}







void insma( sparsemat *M, int row, int column, double value )
//
// inserts a value to sparse matrix M into a specified position
// using absolute addressing
//
// row :    block row
// column : block column
// value :  value to insert
//
{
  int n, j, idx;


  for (n=0, j=M->ia[row]-1; j <= M->ia[row+1]-2; j++, n++)
    if (column+1 == M->ja[j])
    {
      idx = n;
      j = M->ia[row+1]-1;
    }

  M->a[M->ia[row]-1+idx] = value;
}







double getmr( sparsemat *M, int row, int column, int i )
//
// obtains a value from sparse matrix M from a specified position
// using relative addressing
//
// row :    block row
// column : block column
// i :      position in block
// 
// returns value from matrix
//
{
  int n, j, idx, rmdof;


  rmdof = row * M->dof;

  for (n=0, j=M->ia[rmdof+i]-1; j <= M->ia[rmdof+i+1]-2; j++, n++)
    if (column*M->dof+i+1 == M->ja[j])
    {
      idx = n;
      j = M->ia[rmdof+i+1]-1;
    }

  return( M->a[M->ia[rmdof+i]-1+idx] );
}







double getma( sparsemat *M, int row, int column )
//
// obtains a value from sparse matrix M from a specified position
// using absolute addressing
//
// row :    block row
// column : block column
// 
// returns value from matrix
//
{
  int n, j, idx;


  for (n=0, j=M->ia[row]-1; j <= M->ia[row+1]-2; j++, n++)
    if (column+1 == M->ja[j])
    {
      idx = n;
      j = M->ia[row+1]-1;
    }

  return( M->a[M->ia[row]-1+idx] );
}






void dpzero( double *ptr, int size, int nthreads )
//
// parallel zero of an array of doubles
//
{
  int i, myid;

 
  // compute chunk size
  i = size / nthreads;
  
  // zero remaining portion
  memset( ptr + nthreads*i, 0, (size % nthreads)*sizeof(double) );

  #ifdef _OPENMP
  #pragma omp parallel private(myid)
  #endif
  {
    #ifdef _OPENMP
    myid = omp_get_thread_num();
    #else
    myid = 0;
    #endif
    
    // each processor zeros its own portion of the matrix
    memset( ptr + myid*i, 0, i*sizeof(double) );
  }
}







void ipzero( int *ptr, int size, int nthreads )
//
// parallel zero of an array of integers
//
{
  int i, myid;

 
  // compute chunk size
  i = size / nthreads;
  
  // zero remaining portion
  memset( ptr + nthreads*i, 0, (size % nthreads)*sizeof(int) );

  #ifdef _OPENMP
  #pragma omp parallel private(myid)
  #endif
  {
    #ifdef _OPENMP
    myid = omp_get_thread_num();
    #else
    myid = 0;
    #endif
    
    // each processor zeros its own portion of the matrix
    memset( ptr + myid*i, 0, i*sizeof(int) );
  }
}





/*
void ipzero_u( int *ptr, int size, int nthreads )
//
// parallel zero of an array of integers, except the portion belonging to CPU 0
//
{
  int i, myid;

 
  // compute chunk size
  i = size / nthreads;
  
  // zero remaining portion
  memset( ptr + nthreads*i, 0, (size % nthreads)*sizeof(int) );

  #ifdef _OPENMP
  #pragma omp parallel private(myid)
  #endif
  {
    #ifdef _OPENMP
    myid = omp_get_thread_num();
    #else
    myid = 0;
    #endif
    
    // each processor zeros its own portion of the matrix
    memset( ptr + myid*i, 0, i*sizeof(int) );
  }
}
*/







void printmat_as_stored( sparsemat *M )
//
// print out sparse matrix as stored
//
{
  int i;


  printf("int size = %d\n",M->size);
  printf("int rsize = %d\n",M->rsize);
  printf("int dof = %d\n",M->dof);
  printf("int nnz = %d\n",M->nnz);

  printf("int rnz[] = { ");
  for ( i = 0; i < M->size-1; i++ ) printf("%d, ",M->rnz[i]);
  printf("%d };\n",M->rnz[i]);
  
  printf("int ia[] = { ");
  for ( i = 0; i < M->rsize; i++ ) printf("%d, ",M->ia[i]);
  printf("%d };\n",M->ia[i]);
  
  printf("int ja[] = { ");
  for ( i = 0; i < M->nnz-1; i++ ) printf("%d, ",M->ja[i]);
  printf("%d };\n",M->ja[i]);
  
  printf("double a[] = { ");
  for ( i = 0; i < M->nnz-1; i++ ) printf("%.3g, ",M->a[i]);
  printf("%.3g };\n",M->a[i]);
}







void printmat2file_as_stored( sparsemat *M, char *filename, double *rhs )
//
// print out sparse matrix and a vector (rhs) into file as stored
//
// vector rhs should be the same size as M->rsize (no error checkin is done here)
//
{
  int i;
  FILE *fout;


  if (!(fout = fopen(filename,"w"))) ERR("cannot open outputfile for matrix");

  fprintf(fout,"%d\t%d\n",M->rsize,M->nnz);
  
  for ( i = 0; i < M->rsize+1; i++ ) fprintf(fout,"%d\t",M->ia[i]);
  fprintf(fout,"\n");
  
  for ( i = 0; i < M->nnz; i++ ) fprintf(fout,"%d\t",M->ja[i]);
  fprintf(fout,"\n");
  
  for ( i = 0; i < M->nnz; i++ ) fprintf(fout,"%g\t",M->a[i]);
  fprintf(fout,"\n");
  
  for ( i = 0; i < M->rsize; i++ ) fprintf(fout,"%g\t",rhs[i]);
  fprintf(fout,"\n");

  fclose( fout );
}








void printmat_as_structure( sparsemat *M )
//
// print out nonzero structure of sparse matrix
//
{
  int i, j, n;


  for ( i = 0; i < M->rsize; i++ )
  {
    for ( j = 1; j < M->ja[M->ia[i]-1]; j++ ) printf(". ");		// leading zeros
    for ( n = M->ia[i]-1; n < M->ia[i+1]-1; n++ )
    {
      if (n>M->ia[i]-1)
        for ( j = M->ja[n-1]; j < M->ja[n]-1; j++ ) printf(". ");	// zeros between nonzeros
      printf("o ");							// nonzero
    }
    for ( j = M->ja[M->ia[i+1]-2]; j < M->rsize; j++ ) printf(". ");	// trailing zeros
    printf("\n");
  }
}







void printmat_as_matrix( sparsemat *M )
//
// print out sparse matrix as a real matrix
//
{
  int i, j, n;


  for ( i = 0; i < M->rsize; i++ )
  {
    for ( j = 1; j < M->ja[M->ia[i]-1]; j++ ) printf("0\t");
    for ( n = M->ia[i]-1; n < M->ia[i+1]-1; n++ )
    {
      if (n>M->ia[i]-1)
        for ( j = M->ja[n-1]; j < M->ja[n]-1; j++ ) printf("0\t");
      printf("%.3g\t",M->a[n]);
    }
    for ( j = M->ja[M->ia[i+1]-2]; j < M->rsize; j++ ) printf("0\t");
    printf("\n");
  }
}






void printmat_as_matlab( sparsemat *M )
//
// print out sparse matrix as a matlab matrix
//
{ int i, j, n;


  printf("A = [ ");
  for ( i = 0; i < M->rsize; i++ )
  {
    for ( j = 1; j < M->ja[M->ia[i]-1]; j++ ) printf("0 ");
    for ( n = M->ia[i]-1; n < M->ia[i+1]-1; n++ )
    {
      if (n>M->ia[i]-1)
        for ( j = M->ja[n-1]; j < M->ja[n]-1; j++ ) printf("0 ");
      printf("%.3g ",M->a[n]);
    }
    for ( j = M->ja[M->ia[i+1]-2]; j < M->rsize; j++ ) printf("0 ");
    printf(";\n");
  }
  printf("]\n");
}
