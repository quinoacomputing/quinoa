//******************************************************************************
/*!
  \file      src/Base/SymCompressedRowMatrix.h
  \author    J. Bakosi
  \date      Wed 29 Aug 2012 07:25:07 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Symmetric compressed row sparse matrix declaration
  \details   Derived sparse matrix class for symmetric compressed sparse row
             (CSR) storage format, with only the upper triangle stored,
             including the main diagonal.
*/
//******************************************************************************
#ifndef SymCompressedRowMatrix_h
#define SymCompressedRowMatrix_h

namespace Quinoa {

//! Symmetric compressed row sparse matrix class
class SymCompressedRowMatrix : SparseMatrix {

  public:
    //! Constructor
             SymCompressedRowMatrix(int size, int dof, int *psup1, int* psup2);
    //! Destructor
    virtual ~SymCompressedRowMatrix();

  private:
    // Don't permit copy or assignment operators
    SymCompressedRowMatrix(const SymCompressedRowMatrix&);
    SymCompressedRowMatrix& operator=(const SymCompressedRowMatrix&);

    int *m_rnz;   //!< Number of nonzeros of each row, vector size: size
    int *m_ia;    //!< Row pointers, vector size: size*dof+1
    int *m_ja;    //!< Column indices, vector size: nnz
    double *m_a;  //!< Nonzero values, vector size: nnz
};

} // namespace Quinoa

#endif // SymCompressedRowMatrix.h
