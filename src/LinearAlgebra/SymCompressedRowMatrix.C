//******************************************************************************
/*!
  \file      src/Base/SymCompressedRowMatrix.C
  \author    J. Bakosi
  \date      Thu Sep  6 16:37:06 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Symmetric compressed row sparse matrix definition
  \details   Derived sparse matrix class for symmetric compressed sparse row
             (CSR) storage format, with only the upper triangle stored,
             including the main diagonal.
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
#include "SymCompressedRowMatrix.h"

using namespace Quinoa;

SymCompressedRowMatrix::SymCompressedRowMatrix(Memory* memory,
                                               string name,
                                               Int size,
                                               Int dof,
                                               Int *psup1,
                                               Int *psup2) :
  SparseMatrix(name, size, dof)
//******************************************************************************
//  Constructor
//! \details Creates a size x size compressed row sparse matrix with dof degrees
//!          of freedom, ie. the real size will be (size x dof) x (size x dof)
//!          and symmetric, storing only the upper triangle.
//!
//! \param[in]  memory       Pointer to MemoryStore object
//! \param[in]  size         Size of matrix
//! \param[in]  dof          Number of degrees of freedom
//! \param[in]  psup1,psup2  Linked lists storing points surrounding points,
//!                          i.e. the graph of the nonzero structure
//!
//! \author    J. Bakosi
//******************************************************************************
{
  // Store memory store pointer
  m_memory = memory;

  // Allocate array for for storing the nonzeros in each row
  m_rnz = memory->newZeroEntry(size, INT_VAL, SCALAR_VAR, name+"_rnz");
  // Get its raw pointer right away
  Int* rnz = memory->getPtr<Int>(m_rnz);

  // Allocate array for row indices
  m_ia = memory->newZeroEntry(size*dof+1, INT_VAL, SCALAR_VAR, name+"_ia");
  // Get and store its raw pointer right away
  m_pia = memory->getPtr<Int>(m_ia);

  // calculate number of nonzeros in each block row (rnz[]),
  // total number of nonzeros (nnz) and fill row indices (m_pia[])
  m_nnz = 0;
  m_pia[0] = 1;
  for (int i=0; i<size; i++) {
    // add up and store nonzeros of row i
    // (only upper triangular part, matrix is symmetric)
    rnz[i] = 1;
    for ( int j=psup2[i]+1; j<=psup2[i+1]; j++) {
      rnz[i]++;
    }

    // add up total number of nonzeros
    m_nnz += rnz[i]*dof;
    
    // fill up rowindex
    for (int k=0; k<dof; k++)
      m_pia[i*dof+k+1] = m_pia[i*dof+k] + rnz[i];
  }

  // Allocate array for column indices
  m_ja = memory->newZeroEntry(m_nnz, INT_VAL, SCALAR_VAR, name+"_ja");
  // Get and store its raw pointer right away
  m_pja = memory->getPtr<Int>(m_ja);

  // Allocate array for nonzero matrix values
  m_a = memory->newZeroEntry(m_nnz, REAL_VAL, SCALAR_VAR, name+"_a");
  // Get and store its raw pointer right away
  m_pa = m_memory->getPtr<Real>(m_a);

  // Fill column indices
  for (int i=0; i<size; i++) { // loop over all points
    for (int k=0; k<dof; k++) { // loop over all degrees of freedom in a point
      int itmp = i*dof+k;
      m_pja[m_pia[itmp]-1] = itmp+1;    // put in column index of main diagonal
      // loop over all points surrounding point i
      int j = psup2[i]+1;
      for (int n=1; j<=psup2[i+1]; j++) {
        int e = psup1[j]*dof+k+1;
        // put in column index of an off-diagonal
	m_pja[m_pia[itmp]-1+(n++)] = e;
      }
    }
  }

  // (bubble-)Sort column indices
  // loop over all points
  for (int i=0; i<size; i++) {
    // loop over all degrees of freedom in a point
    for (int k=0; k<dof; k++) {
      // loop over all points surrounding point i
      for (int j=psup2[i]+1; j<=psup2[i+1]; j++) {
        // sort column indices of row i
        for (int l=1; l<rnz[i]; l++) {
          for (int e=0; e<rnz[i]-l; e++) {
            if (m_pja[m_pia[i*dof+k]-1+e] > m_pja[m_pia[i*dof+k]+e]) {
              int itmp;
	      SWAP(m_pja[m_pia[i*dof+k]-1+e], m_pja[m_pia[i*dof+k]+e], itmp);
            }
          }
        }
      }
    }
  }
}

SymCompressedRowMatrix::~SymCompressedRowMatrix()
//******************************************************************************
//  Destructor
//! \details Free all memory allocated by SymCompressedRowMatrix in the
//!          constructor when leaving scope
//! \author  J. Bakosi
//******************************************************************************
{
  m_memory->freeEntry(m_rnz);
  m_memory->freeEntry(m_ia);
  m_memory->freeEntry(m_ja);
  m_memory->freeEntry(m_a);
}

void
SymCompressedRowMatrix::add(Int row, Int column, Int i, Real value)
//******************************************************************************
//  Add value to matrix in specified position using relative indexing
//!
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \param[in]  i       relative position in block
//! \param[in]  value   value to add
//!
//! \author  J. Bakosi
//******************************************************************************
{
  Int idx;
  Int rmdof = row*m_dof;

  for (Int n=0, j=m_pia[rmdof+i]-1; j<=m_pia[rmdof+i+1]-2; j++, n++)
    if (column*m_dof+i+1 == m_pja[j]) {
      idx = n;
      j = m_pia[rmdof+i+1]-1;
    }

  #ifdef _OPENMP
  #pragma omp atomic
  #endif
  m_pa[m_pia[rmdof+i]-1+idx] += value;
}

void
SymCompressedRowMatrix::add(Int row, Int column, Real value)
//******************************************************************************
//  Add value to matrix in specified position using absolute indexing
//!
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \param[in]  value   value to add
//!
//! \author  J. Bakosi
//******************************************************************************
{
  Int idx;

  for (Int n=0, j=m_pia[row]-1; j<=m_pia[row+1]-2; j++, n++)
    if (column+1 == m_pja[j]) {
      idx = n;
      j = m_pia[row+1]-1;
    }

  #ifdef _OPENMP
  #pragma omp atomic
  #endif
  m_pa[m_pia[row]-1+idx] += value;
}

void
SymCompressedRowMatrix::ins(Int row, Int column, Int i, Real value)
//******************************************************************************
//  Insert value to matrix in specified position using relative indexing
//!
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \param[in]  i       relative position in block
//! \param[in]  value   value to add
//!
//! \author  J. Bakosi
//******************************************************************************
{
  Int idx;
  Int rmdof = row*m_dof;
  
  for (Int n=0, j=m_pia[rmdof+i]-1; j<=m_pia[rmdof+i+1]-2; j++, n++)
    if (column*m_dof+i+1 == m_pja[j]) {
      idx = n;
      j = m_pia[rmdof+i+1]-1;
    }

  m_pa[m_pia[rmdof+i]-1+idx] = value;
}

void
SymCompressedRowMatrix::ins(Int row, Int column, Real value)
//******************************************************************************
//  Insert value to matrix in specified position using absolute indexing
//!
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \param[in]  value   value to add
//!
//! \author  J. Bakosi
//******************************************************************************
{
  Int idx;

  for (Int n=0, j=m_pia[row]-1; j<=m_pia[row+1]-2; j++, n++)
    if (column+1 == m_pja[j]) {
      idx = n;
      j = m_pia[row+1]-1;
    }

  m_pa[m_pia[row]-1+idx] = value;
}

Real
SymCompressedRowMatrix::get(Int row, Int column, Int i)
//******************************************************************************
//  Get value from matrix from specified position using relative indexing
//!
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \param[in]  i       relative position in block
//! \return             matrix value
//!
//! \author  J. Bakosi
//******************************************************************************
{
  Int idx;
  Int rmdof = row*m_dof;

  for (Int n=0, j=m_pia[rmdof+i]-1; j<=m_pia[rmdof+i+1]-2; j++, n++)
    if (column*m_dof+i+1 == m_pja[j]) {
      idx = n;
      j = m_pia[rmdof+i+1]-1;
    }

  return m_pa[m_pia[rmdof+i]-1+idx];
}

Real
SymCompressedRowMatrix::get(Int row, Int column)
//******************************************************************************
//  Get value from matrix from specified position using absolute indexing
//!
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \return             matrix value
//!
//! \author  J. Bakosi
//******************************************************************************
{
  Int idx;

  for (Int n=0, j=m_pia[row]-1; j<=m_pia[row+1]-2; j++, n++)
    if (column+1 == m_pja[j]) {
      idx = n;
      j = m_pia[row+1]-1;
    }

  return m_pa[m_pia[row]-1+idx];
}

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
