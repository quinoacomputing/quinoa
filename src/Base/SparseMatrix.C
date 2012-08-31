//******************************************************************************
/*!
  \file      src/Base/SparseMatrix.C
  \author    J. Bakosi
  \date      Thu 30 Aug 2012 10:54:54 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Sparse matrix definition
  \details   Sparse matrix base class definition
*/
//******************************************************************************

#include <cstdio>
#include <cstdlib>
#include <string.h>
#ifdef _OPENMP
#include "omp.h"
#endif
#include "Macros.h"
#include "SparseMatrix.h"

using namespace Quinoa;

SparseMatrix::SparseMatrix()
//******************************************************************************
//  Constructor
//! \author    J. Bakosi
//******************************************************************************
{
}

SparseMatrix::~SparseMatrix()
//******************************************************************************
//  Destructor
//! \author    J. Bakosi
//******************************************************************************
{
}

// void destroy_sparsemat(sparsemat *m)
// //
// // frees memory allocated for data structures of sparse matrix m
// //
// {
//   free( m->rnz );
//   free( m->ia );
//   free( m->ja );
//   free( m->a );
//   free( m->r );
//   free( m->p );
//   free( m->z );
//   free( m->q );
//   free( m->d );
//   free( m->u );
//   free( m->lev1 );
//   free( m->lev2 );
// }
// 
// void addmr(sparsemat *M, int row, int column, int i, double value)
// // -----------------------------------------------------------------------------
// // Routine: addmr - Add a value to sparse matrix M into a specified position
// //                  using relative addressing
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // row :    block row
// // column : block column
// // i :      position in block
// // value :  value to add
// // -----------------------------------------------------------------------------
// {
//   int n, j, idx, rmdof;
// 
//   rmdof = row * M->dof;
// 
//   for (n=0, j=M->ia[rmdof+i]-1; j <= M->ia[rmdof+i+1]-2; j++, n++)
//     if (column*M->dof+i+1 == M->ja[j])
//     {
//       idx = n;
//       j = M->ia[rmdof+i+1]-1;
//     }
// 
//   #ifdef _OPENMP
//   #pragma omp atomic
//   #endif
//   M->a[M->ia[rmdof+i]-1+idx] += value;
// }
// 
// void insmr(sparsemat *M, int row, int column, int i, double value)
// // -----------------------------------------------------------------------------
// // Routine: insmr - Insert a value to sparse matrix M into a specified position
// //                  using relative addressing
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // row :    block row
// // column : block column
// // i :      position in block
// // value :  value to insert
// // -----------------------------------------------------------------------------
// {
//   int n, j, idx, rmdof;
// 
//   rmdof = row * M->dof;
//   
//   for (n=0, j=M->ia[rmdof+i]-1; j <= M->ia[rmdof+i+1]-2; j++, n++)
//     if (column*M->dof+i+1 == M->ja[j])
//     {
//       idx = n;
//       j = M->ia[rmdof+i+1]-1;
//     }
// 
//   M->a[M->ia[rmdof+i]-1+idx] = value;
// }
// 
// void addma(sparsemat *M, int row, int column, double value)
// // -----------------------------------------------------------------------------
// // Routine: addma - Add a value to sparse matrix M into a specified position
// //                  using absolute addressing
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // row :    block row
// // column : block column
// // value :  value to add
// // -----------------------------------------------------------------------------
// {
//   int n, j, idx;
// 
//   for (n=0, j=M->ia[row]-1; j <= M->ia[row+1]-2; j++, n++)
//     if (column+1 == M->ja[j])
//     {
//       idx = n;
//       j = M->ia[row+1]-1;
//     }
// 
//   #ifdef _OPENMP
//   #pragma omp atomic
//   #endif
//   M->a[M->ia[row]-1+idx] += value;
// }
// 
// void insma(sparsemat *M, int row, int column, double value)
// // -----------------------------------------------------------------------------
// // Routine: insma - Insert a value to sparse matrix M into a specified position
// //                  using absolute addressing
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // row :    block row
// // column : block column
// // value :  value to insert
// // -----------------------------------------------------------------------------
// {
//   int n, j, idx;
// 
//   for (n=0, j=M->ia[row]-1; j <= M->ia[row+1]-2; j++, n++)
//     if (column+1 == M->ja[j])
//     {
//       idx = n;
//       j = M->ia[row+1]-1;
//     }
// 
//   M->a[M->ia[row]-1+idx] = value;
// }
// 
// double getmr( sparsemat *M, int row, int column, int i )
// // -----------------------------------------------------------------------------
// // Routine: getmr - Obtain a value from sparse matrix M from a specified position
// //                  using relative addressing
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // row :    block row
// // column : block column
// // i :      position in block
// // returns value from matrix
// // -----------------------------------------------------------------------------
// {
//   int n, j, idx, rmdof;
// 
//   rmdof = row * M->dof;
// 
//   for (n=0, j=M->ia[rmdof+i]-1; j <= M->ia[rmdof+i+1]-2; j++, n++)
//     if (column*M->dof+i+1 == M->ja[j])
//     {
//       idx = n;
//       j = M->ia[rmdof+i+1]-1;
//     }
// 
//   return( M->a[M->ia[rmdof+i]-1+idx] );
// }
// 
// double getma( sparsemat *M, int row, int column )
// // -----------------------------------------------------------------------------
// // Routine: getma - Obtain a value from sparse matrix M from a specified position
// //                  using absolute addressing
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // row :    block row
// // column : block column
// // returns value from matrix
// // -----------------------------------------------------------------------------
// {
//   int n, j, idx;
// 
//   for (n=0, j=M->ia[row]-1; j <= M->ia[row+1]-2; j++, n++)
//     if (column+1 == M->ja[j])
//     {
//       idx = n;
//       j = M->ia[row+1]-1;
//     }
// 
//   return( M->a[M->ia[row]-1+idx] );
// }
// 
// void dpzero( double *ptr, int size, int nthreads )
// // -----------------------------------------------------------------------------
// // Routine: dpzero - Parallel zero of an array of doubles
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   int i, myid;
//  
//   // compute chunk size
//   i = size / nthreads;
//   
//   // zero remaining portion
//   memset( ptr + nthreads*i, 0, (size % nthreads)*sizeof(double) );
// 
//   #ifdef _OPENMP
//   #pragma omp parallel private(myid)
//   #endif
//   {
//     #ifdef _OPENMP
//     myid = omp_get_thread_num();
//     #else
//     myid = 0;
//     #endif
//     
//     // each processor zeros its own portion of the matrix
//     memset( ptr + myid*i, 0, i*sizeof(double) );
//   }
// }
// 
// void ipzero( int *ptr, int size, int nthreads )
// // -----------------------------------------------------------------------------
// // Routine: ipzero - Parallel zero of an array of integers
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   int i, myid;
//  
//   // compute chunk size
//   i = size / nthreads;
//   
//   // zero remaining portion
//   memset( ptr + nthreads*i, 0, (size % nthreads)*sizeof(int) );
// 
//   #ifdef _OPENMP
//   #pragma omp parallel private(myid)
//   #endif
//   {
//     #ifdef _OPENMP
//     myid = omp_get_thread_num();
//     #else
//     myid = 0;
//     #endif
//     
//     // each processor zeros its own portion of the matrix
//     memset( ptr + myid*i, 0, i*sizeof(int) );
//   }
// }
// 
// void ipzero_u( int *ptr, int size, int nthreads )
// // -----------------------------------------------------------------------------
// // Routine: ipzero_u - Parallel zero of an array of integers,
// //                     except the portion belonging to CPU 0
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   int i, myid;
//  
//   // compute chunk size
//   i = size / nthreads;
//   
//   // zero remaining portion
//   memset( ptr + nthreads*i, 0, (size % nthreads)*sizeof(int) );
// 
//   #ifdef _OPENMP
//   #pragma omp parallel private(myid)
//   #endif
//   {
//     #ifdef _OPENMP
//     myid = omp_get_thread_num();
//     #else
//     myid = 0;
//     #endif
//     
//     // each processor zeros its own portion of the matrix
//     memset( ptr + myid*i, 0, i*sizeof(int) );
//   }
// }
// 
// void printmat_as_stored( sparsemat *M )
// // -----------------------------------------------------------------------------
// // Routine: ipzero_u - Print out sparse matrix as stored
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   int i;
// 
//   printf("int size = %d\n",M->size);
//   printf("int rsize = %d\n",M->rsize);
//   printf("int dof = %d\n",M->dof);
//   printf("int nnz = %d\n",M->nnz);
// 
//   printf("int rnz[] = { ");
//   for ( i = 0; i < M->size-1; i++ ) printf("%d, ",M->rnz[i]);
//   printf("%d };\n",M->rnz[i]);
//   
//   printf("int ia[] = { ");
//   for ( i = 0; i < M->rsize; i++ ) printf("%d, ",M->ia[i]);
//   printf("%d };\n",M->ia[i]);
//   
//   printf("int ja[] = { ");
//   for ( i = 0; i < M->nnz-1; i++ ) printf("%d, ",M->ja[i]);
//   printf("%d };\n",M->ja[i]);
//   
//   printf("double a[] = { ");
//   for ( i = 0; i < M->nnz-1; i++ ) printf("%.3g, ",M->a[i]);
//   printf("%.3g };\n",M->a[i]);
// }
// 
// void printmat2file_as_stored( sparsemat *M, char *filename, double *rhs )
// // -----------------------------------------------------------------------------
// // Routine: printmat2file_as_stores - Print out sparse matrix and a vector (rhs)
// //                                    into file as stored
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // Vector rhs should be the same size as M->rsize (no error checking is done here)
// {
//   int i;
//   FILE *fout;
// 
//   if (!(fout = fopen(filename,"w"))) ERR("cannot open outputfile for matrix");
// 
//   fprintf(fout,"%d\t%d\n",M->rsize,M->nnz);
//   
//   for ( i = 0; i < M->rsize+1; i++ ) fprintf(fout,"%d\t",M->ia[i]);
//   fprintf(fout,"\n");
//   
//   for ( i = 0; i < M->nnz; i++ ) fprintf(fout,"%d\t",M->ja[i]);
//   fprintf(fout,"\n");
//   
//   for ( i = 0; i < M->nnz; i++ ) fprintf(fout,"%g\t",M->a[i]);
//   fprintf(fout,"\n");
//   
//   for ( i = 0; i < M->rsize; i++ ) fprintf(fout,"%g\t",rhs[i]);
//   fprintf(fout,"\n");
// 
//   fclose( fout );
// }
// 
// void printmat_as_structure( sparsemat *M )
// // -----------------------------------------------------------------------------
// // Routine: printmat_as_structure - Print out nonzero structure of sparse matrix
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   for (int i=0; i<M->rsize; i++) {
// 
//     for (int j=1; j<M->ja[M->ia[i]-1]; j++)
//        printf(". ");// leading zeros
// 
//     for (int n=M->ia[i]-1; n<M->ia[i+1]-1; n++) {
//       if (n>M->ia[i]-1)
//         for (int j=M->ja[n-1]; j<M->ja[n]-1; j++)
//            printf(". "); // zeros between nonzeros
//       printf("o "); // nonzero
//     }
// 
//     for (int j=M->ja[M->ia[i+1]-2]; j<M->rsize; j++)
//        printf(". "); // trailing zeros
// 
//     printf("\n");
//   }
// }
// 
// void printmat_as_matrix( sparsemat *M )
// // -----------------------------------------------------------------------------
// // Routine: printmat_as_matrix - Print out sparse matrix as a real matrix
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   for (int i=0; i<M->rsize; i++) {
//     for (int j=1; j<M->ja[M->ia[i]-1]; j++)
//        printf("0\t");
// 
//     for (int n=M->ia[i]-1; n<M->ia[i+1]-1; n++) {
//       if (n>M->ia[i]-1)
//         for (int j=M->ja[n-1]; j<M->ja[n]-1; j++)
//           printf("0\t");
//       printf("%.3g\t",M->a[n]);
//     }
// 
//     for (int j=M->ja[M->ia[i+1]-2]; j<M->rsize; j++)
//       printf("0\t");
// 
//     printf("\n");
//   }
// }
// 
// void printmat_as_matlab( sparsemat *M )
// // -----------------------------------------------------------------------------
// // Routine: printmat_as_matrix - Print out sparse matrix as a matlab matrix
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   printf("A = [ ");
//   for (int i=0; i<M->rsize; i++) {
//     for (int j=1; j<M->ja[M->ia[i]-1]; j++)
//        printf("0 ");
// 
//     for (int n=M->ia[i]-1; n<M->ia[i+1]-1; n++) {
//       if (n>M->ia[i]-1)
//         for (int j=M->ja[n-1]; j<M->ja[n]-1; j++)
//           printf("0 ");
//       printf("%.3g ",M->a[n]);
//     }
// 
//     for (int j=M->ja[M->ia[i+1]-2]; j<M->rsize; j++)
//       printf("0 ");
// 
//     printf(";\n");
//   }
//   printf("]\n");
// }
