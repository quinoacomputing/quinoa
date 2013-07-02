//******************************************************************************
/*!
  \file      src/LinearAlgebra/SymCompRowMatrix.h
  \author    J. Bakosi
  \date      Tue Jul  2 15:37:45 2013
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
class SymCompRowMatrix : public SparseMatrix {

  public:
    //! Constructor
    explicit SymCompRowMatrix(Memory* const memory,
                              const std::string name,
                              const int size,
                              const int dof,
                              const int *psup1,
                              const int* psup2);
    //! Destructor
    virtual ~SymCompRowMatrix() noexcept;

    //! Add value to matrix in specified position using relative indexing
    void add(int row, int column, int i, real value);
    //! Add value to matrix in specified position using absolute indexing
    void add(int row, int column, real value);

    //! Insert value to matrix in specified position using relative indexing
    void ins(int row, int column, int i, real value);
    //! Insert value to matrix in specified position using absolute indexing
    void ins(int row, int column, real value);

    //! Get value from matrix from specified position using relative indexing
    real get(int row, int column, int i) const;
    //! Get value from matrix from specified position using absolute indexing
    real get(int row, int column) const;

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

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

    Data<int> m_rnz;  //!< Nonzeros in each row (only in constructor)
    Data<int> m_ia;   //!< Row indices, vector size: size*dof+1
    Data<int> m_ja;   //!< Column indices, vector size: nnz
    Data<real> m_a;   //!< Nonzero matrix values, vector size: nnz
};

} // namespace Quinoa

#endif // SymCompRowMatrix_h
