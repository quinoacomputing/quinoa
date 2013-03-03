//******************************************************************************
/*!
  \file      src/LinearAlgebra/SymCompRowMatrix.h
  \author    J. Bakosi
  \date      Sun 03 Mar 2013 06:55:58 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Symmetric compressed row sparse matrix
  \details   Derived sparse matrix class for symmetric compressed sparse row
             (CSR) storage format, with only the upper triangle stored,
             including the main diagonal.
*/
//******************************************************************************
#ifndef SymCompRowMatrix_h
#define SymCompRowMatrix_h

#include <Memory.h>

namespace Quinoa {

//! Symmetric compressed row sparse matrix
class SymCompRowMatrix : private SparseMatrix {

  public:
    //! Constructor
    SymCompRowMatrix(Memory* memory,
                     string name,
                     int size,
                     int dof,
                     int *psup1,
                     int* psup2);
    //! Destructor
    ~SymCompRowMatrix();

    //! Add value to matrix in specified position using relative indexing
    void add(int row, int column, int i, real value);
    //! Add value to matrix in specified position using absolute indexing
    void add(int row, int column, real value);

    //! Insert value to matrix in specified position using relative indexing
    void ins(int row, int column, int i, real value);
    //! Insert value to matrix in specified position using absolute indexing
    void ins(int row, int column, real value);

    //! Get value from matrix from specified position using relative indexing
    real get(int row, int column, int i);
    //! Get value from matrix from specified position using absolute indexing
    real get(int row, int column);

    //! Print out matrix entries as stored
    void echoAsStored(ostream& ofs);

    //! Print out nonzero structure of matrix
    void echoNonzeroStructure(ostream& ofs);

    //! Print out matrix as a real matrix
    void echoAsMatrix(ostream& ofs);

    //! Print out matrix as a matlab matrix
    void echoAsMatlab(ostream& ofs);

  private:
    //! Don't permit copy constructor
    SymCompRowMatrix(const SymCompRowMatrix&) = delete;
    //! Don't permit copy assigment
    SymCompRowMatrix& operator=(const SymCompRowMatrix&) = delete;
    //! Don't permit move constructor
    SymCompRowMatrix(SymCompRowMatrix&&) = delete;
    //! Don't permit move assigment
    SymCompRowMatrix& operator=(SymCompRowMatrix&&) = delete;

    Data<int> m_ia;   //!< Row indices, vector size: size*dof+1
    Data<int> m_ja;   //!< Column indices, vector size: nnz
    Data<real> m_a;   //!< Nonzero matrix values, vector size: nnz
};

} // namespace Quinoa

#endif // SymCompRowMatrix_h
