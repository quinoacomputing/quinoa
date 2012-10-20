//******************************************************************************
/*!
  \file      src/LinearAlgebra/SymCompRowMatrix.C
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 04:17:46 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Symmetric compressed row sparse matrix
  \details   Derived sparse matrix class for symmetric compressed sparse row
             (CSR) storage format, with only the upper triangle stored,
             including the main diagonal.
*/
//******************************************************************************

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <SparseMatrix.h>
#include <SymCompRowMatrix.h>
#include <Macros.h>

using namespace Quinoa;

SymCompRowMatrix::SymCompRowMatrix(Memory* memory,
                                   string name,
                                   int size,
                                   int dof,
                                   int *psup1,
                                   int *psup2) : SparseMatrix(name, size, dof)
//******************************************************************************
//  Constructor
//! \details Creates a size x size compressed row sparse matrix with dof degrees
//!          of freedom, ie. the real size will be (size x dof) x (size x dof)
//!          and symmetric, storing only the upper triangle.
//! \param[in]  memory       Pointer to MemoryStore object
//! \param[in]  name         Name of the SparseMatrix instance
//! \param[in]  size         Size of matrix
//! \param[in]  dof          Number of degrees of freedom
//! \param[in]  psup1,psup2  Linked lists storing points surrounding points,
//!                          i.e. the graph of the nonzero structure
//! \author    J. Bakosi
//******************************************************************************
{
  // Store memory store pointer
  m_memory = memory;

  // Allocate array for storing the nonzeros in each row
  MemoryEntry* mrnz =
    memory->newZeroEntry(size, ValType::INT, VarType::SCALAR, name+"_rnz");
  // Get its raw pointer right away
  int* rnz = memory->getPtr<int>(mrnz);

  // Allocate array for row indices
  m_ia =
    memory->newZeroEntry(m_rsize+1, ValType::INT, VarType::SCALAR, name+"_ia");
  // Get and store its raw pointer right away
  m_pia = memory->getPtr<int>(m_ia);

  // Calculate number of nonzeros in each block row (rnz[]),
  // total number of nonzeros (nnz) and fill row indices (m_pia[])
  m_nnz = 0;
  m_pia[0] = 1;
  for (int i=0; i<size; i++) {
    // add up and store nonzeros of row i
    // (only upper triangular part, matrix is symmetric)
    rnz[i] = 1;
    for ( int j=psup2[i]+1; j<=psup2[i+1]; j++) {
      rnz[i]++;
    }

    // add up total number of nonzeros
    m_nnz += rnz[i]*dof;
    
    // fill up rowindex
    for (int k=0; k<dof; k++)
      m_pia[i*dof+k+1] = m_pia[i*dof+k] + rnz[i];
  }

  // Allocate array for column indices
  m_ja = memory->newZeroEntry(m_nnz, ValType::INT, VarType::SCALAR, name+"_ja");
  // Get and store its raw pointer right away
  m_pja = memory->getPtr<int>(m_ja);

  // Allocate array for nonzero matrix values
  m_a = memory->newZeroEntry(m_nnz, ValType::REAL, VarType::SCALAR, name+"_a");
  // Get and store its raw pointer right away
  m_pa = m_memory->getPtr<real>(m_a);

  // Fill column indices
  for (int i=0; i<size; i++) { // loop over all points
    for (int k=0; k<dof; k++) { // loop over all degrees of freedom in a point
      int itmp = i*dof+k;
      m_pja[m_pia[itmp]-1] = itmp+1;    // put in column index of main diagonal
      // loop over all points surrounding point i
      int j = psup2[i]+1;
      for (int n=1; j<=psup2[i+1]; j++) {
        int e = psup1[j]*dof+k+1;
        // put in column index of an off-diagonal
	m_pja[m_pia[itmp]-1+(n++)] = e;
      }
    }
  }

  // (bubble-)Sort column indices
  // loop over all points
  for (int i=0; i<size; i++) {
    // loop over all degrees of freedom in a point
    for (int k=0; k<dof; k++) {
      // loop over all points surrounding point i
      for (int j=psup2[i]+1; j<=psup2[i+1]; j++) {
        // sort column indices of row i
        for (int l=1; l<rnz[i]; l++) {
          for (int e=0; e<rnz[i]-l; e++) {
            if (m_pja[m_pia[i*dof+k]-1+e] > m_pja[m_pia[i*dof+k]+e]) {
              int itmp;
	      SWAP(m_pja[m_pia[i*dof+k]-1+e], m_pja[m_pia[i*dof+k]+e], itmp);
            }
          }
        }
      }
    }
  }

  // Free array for storing the nonzeros in each row
  m_memory->freeEntry(mrnz);
}

SymCompRowMatrix::~SymCompRowMatrix()
//******************************************************************************
//  Destructor
//! \details Free all memory allocated by SymCompRowMatrix in the
//!          constructor when leaving scope
//! \author  J. Bakosi
//******************************************************************************
{
  try {
    m_memory->freeEntry(m_ia);
    m_memory->freeEntry(m_ja);
    m_memory->freeEntry(m_a);
    // No exception leaves a destructor: if any of the above calls throws and
    // exception, e.g. m_ia points to an unallocated entry, a MemoryException
    // is thrown, caught inside here and we only emit a warning. This ensures
    // that terminate is not called and that we finish a potentially already
    // propagating exception.
  } catch (...) {
    cerr << "WARNING: Exception in SymCompRowMatrix::~SymCompRowMatrix" << endl;
  }
}

void
SymCompRowMatrix::add(int row, int column, int i, real value)
//******************************************************************************
//  Add value to matrix in specified position using relative indexing
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \param[in]  i       relative position in block
//! \param[in]  value   value to add
//! \author  J. Bakosi
//******************************************************************************
{
  int idx = 0;
  int rmdof = row*m_dof;

  for (int n=0, j=m_pia[rmdof+i]-1; j<=m_pia[rmdof+i+1]-2; j++, n++)
    if (column*m_dof+i+1 == m_pja[j]) {
      idx = n;
      j = m_pia[rmdof+i+1]-1;
    }

  #ifdef _OPENMP
  #pragma omp atomic
  #endif
  m_pa[m_pia[rmdof+i]-1+idx] += value;
}

void
SymCompRowMatrix::add(int row, int column, real value)
//******************************************************************************
//  Add value to matrix in specified position using absolute indexing
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \param[in]  value   value to add
//! \author  J. Bakosi
//******************************************************************************
{
  int idx = 0;

  for (int n=0, j=m_pia[row]-1; j<=m_pia[row+1]-2; j++, n++)
    if (column+1 == m_pja[j]) {
      idx = n;
      j = m_pia[row+1]-1;
    }

  #ifdef _OPENMP
  #pragma omp atomic
  #endif
  m_pa[m_pia[row]-1+idx] += value;
}

void
SymCompRowMatrix::ins(int row, int column, int i, real value)
//******************************************************************************
//  Insert value to matrix in specified position using relative indexing
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \param[in]  i       relative position in block
//! \param[in]  value   value to add
//! \author  J. Bakosi
//******************************************************************************
{
  int idx = 0;
  int rmdof = row*m_dof;
  
  for (int n=0, j=m_pia[rmdof+i]-1; j<=m_pia[rmdof+i+1]-2; j++, n++)
    if (column*m_dof+i+1 == m_pja[j]) {
      idx = n;
      j = m_pia[rmdof+i+1]-1;
    }

  m_pa[m_pia[rmdof+i]-1+idx] = value;
}

void
SymCompRowMatrix::ins(int row, int column, real value)
//******************************************************************************
//  Insert value to matrix in specified position using absolute indexing
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \param[in]  value   value to add
//! \author  J. Bakosi
//******************************************************************************
{
  int idx = 0;

  for (int n=0, j=m_pia[row]-1; j<=m_pia[row+1]-2; j++, n++)
    if (column+1 == m_pja[j]) {
      idx = n;
      j = m_pia[row+1]-1;
    }

  m_pa[m_pia[row]-1+idx] = value;
}

real
SymCompRowMatrix::get(int row, int column, int i)
//******************************************************************************
//  Get value from matrix from specified position using relative indexing
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \param[in]  i       relative position in block
//! \return             matrix value
//! \author  J. Bakosi
//******************************************************************************
{
  int idx = 0;
  int rmdof = row*m_dof;

  for (int n=0, j=m_pia[rmdof+i]-1; j<=m_pia[rmdof+i+1]-2; j++, n++)
    if (column*m_dof+i+1 == m_pja[j]) {
      idx = n;
      j = m_pia[rmdof+i+1]-1;
    }

  return m_pa[m_pia[rmdof+i]-1+idx];
}

real
SymCompRowMatrix::get(int row, int column)
//******************************************************************************
//  Get value from matrix from specified position using absolute indexing
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \return             matrix value
//! \author  J. Bakosi
//******************************************************************************
{
  int idx = 0;

  for (int n=0, j=m_pia[row]-1; j<=m_pia[row+1]-2; j++, n++)
    if (column+1 == m_pja[j]) {
      idx = n;
      j = m_pia[row+1]-1;
    }

  return m_pa[m_pia[row]-1+idx];
}

void
SymCompRowMatrix::echoAsStored(ostream& ofs)
//******************************************************************************
//  Print out matrix entries as stored
//! \param[in]  ofs  output stream
//! \author  J. Bakosi
//******************************************************************************
{
  ofs << "name  = " << m_name << endl;
  ofs << "size  = " << m_size << endl;
  ofs << "rsize = " << m_rsize << endl;
  ofs << "dof   = " << m_dof << endl;
  ofs << "nnz   = " << m_nnz << endl;

  int i;
  ofs << "ia[] = { ";
  for (i=0; i<m_rsize; i++) ofs << m_pia[i];
  ofs << m_pia[i] << endl;

  ofs << "ja[] = { ";
  for (i=0; i<m_nnz-1; i++) ofs << m_pja[i];
  ofs << m_pja[i] << endl;

  ofs << "a[] = { ";
  for (i=0; i<m_nnz-1; i++) ofs << m_pa[i];
  ofs << m_pa[i] << endl;
}

void
SymCompRowMatrix::echoNonzeroStructure(ostream& ofs)
//******************************************************************************
//  Print out nonzero structure of matrix
//! \param[in]  ofs  output stream
//! \author  J. Bakosi
//******************************************************************************
{
  for (int i=0; i<m_rsize; i++) {
    for (int j=1; j<m_pja[m_pia[i]-1]; j++)  // leading zeros
       ofs << ". ";

    for (int n=m_pia[i]-1; n<m_pia[i+1]-1; n++) {
      if (n>m_pia[i]-1) {  // zeros between nonzeros
        for (int j=m_pja[n-1]; j<m_pja[n]-1; j++)
           ofs << ". ";
      }
      ofs << "o ";  // nonzero
    }

    for (int j=m_pja[m_pia[i+1]-2]; j<m_rsize; j++)
       ofs << ". ";  // trailing zeros

    ofs << endl;
  }
}

void
SymCompRowMatrix::echoAsMatrix(ostream& ofs)
//******************************************************************************
//  Print out matrix as a real matrix
//! \param[in]  ofs  output stream
//! \author  J. Bakosi
//******************************************************************************
{
  for (int i=0; i<m_rsize; i++) {
    for (int j=1; j<m_pja[m_pia[i]-1]; j++)
       ofs << "0\t";

    for (int n=m_pia[i]-1; n<m_pia[i+1]-1; n++) {
      if (n>m_pia[i]-1) {
        for (int j=m_pja[n-1]; j<m_pja[n]-1; j++)
          ofs << "0\t";
      }
      ofs << m_pa[n] << "\t";
    }

    for (int j=m_pja[m_pia[i+1]-2]; j<m_rsize; j++)
      ofs << "0\t";

    ofs << endl;
  }
}

void
SymCompRowMatrix::echoAsMatlab(ostream& ofs)
//******************************************************************************
//  Print out matrix as a maltab matrix
//! \param[in]  ofs  output stream
//! \author  J. Bakosi
//******************************************************************************
{
  ofs << "A = [ ";
  for (int i=0; i<m_rsize; i++) {
    for (int j=1; j<m_pja[m_pia[i]-1]; j++)
       ofs << "0 ";

    for (int n=m_pia[i]-1; n<m_pia[i+1]-1; n++) {
      if (n>m_pia[i]-1)
        for (int j=m_pja[n-1]; j<m_pja[n]-1; j++)
          ofs << "0 ";
      ofs << m_pa[n] << " ";
    }

    for (int j=m_pja[m_pia[i+1]-2]; j<m_rsize; j++)
      ofs << "0 ";

    ofs << ";" << endl;
  }
  ofs << "]" << endl;
}
