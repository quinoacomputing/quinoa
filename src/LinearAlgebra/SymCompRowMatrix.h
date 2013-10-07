//******************************************************************************
/*!
  \file      src/LinearAlgebra/SymCompRowMatrix.h
  \author    J. Bakosi
  \date      Mon Oct  7 10:18:41 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Symmetric compressed row sparse matrix
  \details   Derived sparse matrix class for symmetric compressed sparse row
             (CSR) storage format, with only the upper triangle stored,
             including the main diagonal.
*/
//******************************************************************************
#ifndef SymCompRowMatrix_h
#define SymCompRowMatrix_h

#include <SparseMatrix.h>

namespace tk {

//! Symmetric compressed row sparse matrix
class SymCompRowMatrix : public SparseMatrix {

  public:
    //! Constructor
    explicit SymCompRowMatrix(int size,
                              int dof,
                              const int *psup1,
                              const int* psup2);

    //! Add value to matrix in specified position using relative indexing
    void add(int row, int column, int i, tk::real value);
    //! Add value to matrix in specified position using absolute indexing
    void add(int row, int column, tk::real value);

    //! Insert value to matrix in specified position using relative indexing
    void ins(int row, int column, int i, tk::real value);
    //! Insert value to matrix in specified position using absolute indexing
    void ins(int row, int column, tk::real value);

    //! Get value from matrix from specified position using relative indexing
    tk::real get(int row, int column, int i) const;
    //! Get value from matrix from specified position using absolute indexing
    tk::real get(int row, int column) const;

    //! Print out matrix entries as stored
    void echoAsStored(std::ostream& ofs) const;

    //! Print out nonzero structure of matrix
    void echoNonzeroStructure(std::ostream& ofs) const;

    //! Print out matrix as a real matrix
    void echoAsMatrix(std::ostream& ofs) const;

    //! Print out matrix as a matlab matrix
    void echoAsMatlab(std::ostream& ofs) const;

  private:
    //! Don't permit copy constructor
    SymCompRowMatrix(const SymCompRowMatrix&) = delete;
    //! Don't permit copy assigment
    SymCompRowMatrix& operator=(const SymCompRowMatrix&) = delete;
    //! Don't permit move constructor
    SymCompRowMatrix(SymCompRowMatrix&&) = delete;
    //! Don't permit move assigment
    SymCompRowMatrix& operator=(SymCompRowMatrix&&) = delete;

    std::unique_ptr<int[]> m_ia;    //!< Row indices, vector size: size*dof+1
    std::unique_ptr<int[]> m_ja;    //!< Column indices, vector size: nnz
    std::unique_ptr<tk::real[]> m_a;//!< Nonzero matrix values, vector size: nnz
};

} // tk::

#endif // SymCompRowMatrix_h
