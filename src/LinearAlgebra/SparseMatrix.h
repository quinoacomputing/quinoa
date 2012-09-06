//******************************************************************************
/*!
  \file      src/Base/SparseMatrix.h
  \author    J. Bakosi
  \date      Thu Sep  6 15:37:11 2012
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

// void printmat_as_stored( sparsemat *M );
// void printmat_as_structure( sparsemat *M );
// void printmat_as_matrix( sparsemat *M );
// void printmat_as_matlab( sparsemat *M );
// void printmat2file_as_stored( sparsemat *M, char *filename, double *rhs );

} // namespace Quinoa

#endif // SparseMatrix.h
