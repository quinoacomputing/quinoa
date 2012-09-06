//******************************************************************************
/*!
  \file      src/Base/SymCompressedRowMatrix.h
  \author    J. Bakosi
  \date      Wed 05 Sep 2012 08:57:13 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Symmetric compressed row sparse matrix declaration
  \details   Derived sparse matrix class for symmetric compressed sparse row
             (CSR) storage format, with only the upper triangle stored,
             including the main diagonal.
*/
//******************************************************************************
#ifndef SymCompressedRowMatrix_h
#define SymCompressedRowMatrix_h

#include <Memory.h>

namespace Quinoa {

//! Symmetric compressed row sparse matrix class
class SymCompressedRowMatrix : SparseMatrix {

  public:
    //! Constructor
    SymCompressedRowMatrix(Memory* memory,
                           string name,
                           Int size,
                           Int dof,
                           Int *psup1,
                           Int* psup2);
    //! Destructor
    ~SymCompressedRowMatrix() {};

  private:
    // Don't permit copy or assignment operators
    SymCompressedRowMatrix(const SymCompressedRowMatrix&);
    SymCompressedRowMatrix& operator=(const SymCompressedRowMatrix&);

    MemoryEntry *m_rnz;  //!< Number of nonzeros of each row, vector size: size
    MemoryEntry *m_ia;   //!< Row pointers, vector size: size*dof+1
    MemoryEntry *m_ja;   //!< Column indices, vector size: nnz
    MemoryEntry *m_a;   //!< Nonzero values, vector size: nnz
};

} // namespace Quinoa

#endif // SymCompressedRowMatrix.h
