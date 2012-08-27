//******************************************************************************
/*!
  \file      src/Base/SparseMatrix.h
  \author    J. Bakosi
  \date      Sun 26 Aug 2012 09:02:02 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Sparse matrix declaration
  \details   Sparse matrix base class declaration
*/
//******************************************************************************
#ifndef SparseMatrix_h
#define SparseMatrix_h

namespace Quinoa {

//! Sparse matrix base class
class SparseMatrix {

  public:

    //! Constructor
             SparseMatrix(int size, int dof, int *psup1, int* psup2);
    //! Destructor
    virtual ~SparseMatrix();

  private:

    // Don't permit copy or assignment operators
    SparseMatrix(const SparseMatrix&);
    SparseMatrix& operator=(const SparseMatrix&);

    int m_size;   //!< Size of matrix: (dof x size) x (dof x size)
    int m_rsize;  //!< Width of matrix: dof x size
    int m_dof;    //!< Number of degrees of freedom
    int m_nnz;    //!< Total number of nonzeros
    int *m_rnz;   //!< Number of nonzeros of each row, vector size: size
    int *m_ia;    //!< Row pointers, vector size: size*dof+1
    int *m_ja;    //!< Column indices, vector size: nnz
    double *m_a;  //!< Nonzero values, vector size: nnz

//    double *m_r, *m_p, *m_z, *m_q, *m_d, *m_u;
//             //!< sizes: rsize, auxiliary vectors for conjugate gradients
//    int m_nlev; //!< size: 1, number of levels for level scheduling
//    int *m_lev1, *m_lev2;//!< sizes: rsize,nlev+1 linked lists for level scheduling
};

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

} // namespace Quinoa

#endif // SparseMatrix.h
