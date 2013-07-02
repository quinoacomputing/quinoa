//******************************************************************************
/*!
  \file      src/LinearAlgebra/SparseMatrix.h
  \author    J. Bakosi
  \date      Tue Jul  2 15:36:29 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Sparse matrix declaration
  \details   Sparse matrix base class declaration
*/
//******************************************************************************
#ifndef SparseMatrix_h
#define SparseMatrix_h

#include <string>

#include <QuinoaTypes.h>
#include <Memory.h>

namespace Quinoa {

//! Sparse matrix base class
class SparseMatrix {

  protected:
    //! Constructor
    explicit SparseMatrix(Memory* const memory,
                          const std::string name,
                          const int size,
                          const int dof) :
      m_memory(memory),
      m_name(name),
      m_size(size),
      m_rsize(size*dof),
      m_dof(dof) {}

    //! Destructor
    virtual ~SparseMatrix() noexcept = 0;

    Memory* const m_memory;  //!< Local copy of the memory store pointer
    const std::string m_name;//!< Name of the sparse matrix instance
    const int m_size;        //!< Size of matrix: (dof x size) x (dof x size)
    const int m_rsize;       //!< Width of matrix: dof x size
    const int m_dof;         //!< Number of degrees of freedom

    int m_nnz;               //!< Total number of nonzeros

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

} // namespace Quinoa

#endif // SparseMatrix_h
