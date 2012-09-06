//******************************************************************************
/*!
  \file      src/Base/SparseMatrix.h
  \author    J. Bakosi
  \date      Thu Sep  6 14:52:10 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Sparse matrix declaration
  \details   Sparse matrix base class declaration
*/
//******************************************************************************
#ifndef SparseMatrix_h
#define SparseMatrix_h

#include <string>

using namespace std;

#include <QuinoaTypes.h>
#include <Memory.h>

namespace Quinoa {

//! Sparse matrix base class
class SparseMatrix {

  public:
    //! Constructor
    SparseMatrix(string name, Int size, Int dof) :
      m_name(name), m_dof(dof), m_size(size), m_rsize(size*dof) {}
    //! Destructor
    virtual ~SparseMatrix();

  protected:
    Memory *m_memory;  //!< Local copy of the memory store pointer
    string m_name;     //!< Name of the sparse matrix instance
    Int m_dof;         //!< Number of degrees of freedom
    Int m_nnz;         //!< Total number of nonzeros

  private:
    //! Don't permit copy operator
    SparseMatrix(const SparseMatrix&);
    //! Don't permit assigment operator
    SparseMatrix& operator=(const SparseMatrix&);

    Int m_size;   //!< Size of matrix: (dof x size) x (dof x size)
    Int m_rsize;  //!< Width of matrix: dof x size
};

// Free memory allocated for data structures of sparse matrix
// void destroy_sparsemat( sparsemat *m );
// 
// void addmr( sparsemat *M, int row, int column, int i, double value );
// void addma( sparsemat *M, int row, int column, double value );
// void insmr( sparsemat *M, int row, int column, int i, double value );
// void insma( sparsemat *M, int row, int column, double value );
// double getmr( sparsemat *M, int row, int column, int i );
// double getma( sparsemat *M, int row, int column );
// void dpzero( double *ptr, int size, int nthreads );
// void ipzero( int *ptr, int size, int nthreads );
// void printmat_as_stored( sparsemat *M );
// void printmat_as_structure( sparsemat *M );
// void printmat_as_matrix( sparsemat *M );
// void printmat_as_matlab( sparsemat *M );
// void printmat2file_as_stored( sparsemat *M, char *filename, double *rhs );

} // namespace Quinoa

#endif // SparseMatrix.h
