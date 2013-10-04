//******************************************************************************
/*!
  \file      src/LinearAlgebra/SparseMatrix.h
  \author    J. Bakosi
  \date      Thu 03 Oct 2013 09:20:16 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Sparse matrix declaration
  \details   Sparse matrix base class declaration
*/
//******************************************************************************
#ifndef SparseMatrix_h
#define SparseMatrix_h

#include <QuinoaTypes.h>

namespace quinoa {

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

} // namespace quinoa

#endif // SparseMatrix_h
