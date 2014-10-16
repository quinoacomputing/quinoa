//******************************************************************************
/*!
  \file      src/LinearAlgebra/SparseMatrix.h
  \author    J. Bakosi
  \date      Mon Oct  7 10:18:12 2013
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Sparse matrix declaration
  \details   Sparse matrix base class declaration
*/
//******************************************************************************
#ifndef SparseMatrix_h
#define SparseMatrix_h

#include <Types.h>

namespace tk {

//! Sparse matrix base class
class SparseMatrix {

  protected:
    //! Constructor
    explicit SparseMatrix(int size, int dof) :
      m_size(size),
      m_rsize(size*dof),
      m_dof(dof) {}

    const int m_size;          //!< Size of matrix: (dof x size) x (dof x size)
    const int m_rsize;         //!< Width of matrix: dof x size
    const int m_dof;           //!< Number of degrees of freedom

    int m_nnz;                 //!< Total number of nonzeros

  private:
    //! Don't permit copy constructor
    SparseMatrix(const SparseMatrix&) = delete;
    //! Don't permit copy assigment
    SparseMatrix& operator=(const SparseMatrix&) = delete;
    //! Don't permit move constructor
    SparseMatrix(SparseMatrix&&) = delete;
    //! Don't permit move assigment
    SparseMatrix& operator=(SparseMatrix&&) = delete;
};

} // tk::

#endif // SparseMatrix_h
