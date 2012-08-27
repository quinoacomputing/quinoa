//******************************************************************************
/*!
  \file      src/Base/SparseMatrixCSR.h
  \brief     Symmetric compressed row sparse matrix class.
  \details   Derived sparse matrix class for symmetric compressed sparse row
             (CSR) storage format, with only the upper triangle stored,
             including the main diagonal.
  \author    J. Bakosi
  \date      Sat Aug 25 22:54:48 MDT 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
*/
//******************************************************************************

//! Spparse matrix base class
struct sparsemat {
      int size; //!< size: 1, matrix of (dof x size) x (dof x size)
      int rsize;//!< size: 1, real size of matrix (rsize = dof x size)
      int dof;  //!< size: 1, number of unknowns/node (degree of freedom)
      int nnz;  //!< size: 1, total number of nonzeros
      int *rnz; //!< size: size, number of nonzeros of each row
      int *ia;  //!< size: size*dof+1, row pointers
      int *ja;  //!< size: nnz, column indices
      double *a;//!< size: nnz, nonzero values
      double *r, *p, *z, *q, *d, *u;
                //!< sizes: rsize, auxiliary vectors for conjugate gradients
      int nlev; //!< size: 1, number of levels for level scheduling
      int *lev1, *lev2;	//!< sizes: rsize,nlev+1 linked lists for level scheduling
};

//! Initialize sparse matrix
void create_sparsemat( sparsemat *m, int DOF, int size, int *psup1, int *psup2 );

//! Free memory allocated for data structures of sparse matrix
void destroy_sparsemat( sparsemat *m );

void addmr( sparsemat *M, int row, int column, int i, double value );
void addma( sparsemat *M, int row, int column, double value );
void insmr( sparsemat *M, int row, int column, int i, double value );
void insma( sparsemat *M, int row, int column, double value );
double getmr( sparsemat *M, int row, int column, int i );
double getma( sparsemat *M, int row, int column );
void dpzero( double *ptr, int size, int nthreads );
void ipzero( int *ptr, int size, int nthreads );
void printmat_as_stored( sparsemat *M );
void printmat_as_structure( sparsemat *M );
void printmat_as_matrix( sparsemat *M );
void printmat_as_matlab( sparsemat *M );
void printmat2file_as_stored( sparsemat *M, char *filename, double *rhs );
