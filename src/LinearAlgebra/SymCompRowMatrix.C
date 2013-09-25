//******************************************************************************
/*!
  \file      src/LinearAlgebra/SymCompRowMatrix.C
  \author    J. Bakosi
  \date      Wed Sep 25 14:31:05 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Symmetric compressed row sparse matrix
  \details   Derived sparse matrix class for symmetric compressed sparse row
             (CSR) storage format, with only the upper triangle stored,
             including the main diagonal.
*/
//******************************************************************************

#include <iostream>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Macro.h>
#include <SparseMatrix.h>
#include <SymCompRowMatrix.h>

using namespace quinoa;

SymCompRowMatrix::SymCompRowMatrix(const std::string& name,
                                   int size,
                                   int dof,
                                   const int *psup1,
                                   const int *psup2) :
  SparseMatrix(name, size, dof)
//******************************************************************************
//  Constructor
//! \details Creates a size x size compressed row sparse matrix with dof degrees
//!          of freedom, ie. the real size will be (size x dof) x (size x dof)
//!          and symmetric, storing only the upper triangle.
//! \param[in]  name         Name of the SparseMatrix instance
//! \param[in]  size         Size of matrix
//! \param[in]  dof          Number of degrees of freedom
//! \param[in]  psup1,psup2  Linked lists storing points surrounding points,
//!                          i.e. the graph of the nonzero structure
//! \author    J. Bakosi
//******************************************************************************
{
  // Allocate array for storing the nonzeros in each row
  // TODO: Converted to unique_ptr, but still need to zero!
  std::unique_ptr<int[]> rnz(new int [size]);

  // Allocate array for row indices
  // TODO: Converted to unique_ptr, but still need to zero!
  m_ia = std::unique_ptr<int[]>(new int [m_rsize+1]);

  // Calculate number of nonzeros in each block row (rnz[]),
  // total number of nonzeros (nnz) and fill row indices (m_ia[])
  m_nnz = 0;
  m_ia[0] = 1;
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
      m_ia[i*dof+k+1] = m_ia[i*dof+k] + rnz[i];
  }

  // Allocate array for column indices
  // TODO: Converted to unique_ptr, but still need to zero!
  m_ja = std::unique_ptr<int[]>(new int [m_nnz]);
  // Allocate array for nonzero matrix values
  // TODO: Converted to unique_ptr, but still need to zero!
  m_a = std::unique_ptr<real[]>(new real [m_nnz]);

  // Fill column indices
  for (int i=0; i<size; i++) { // loop over all points
    for (int k=0; k<dof; k++) { // loop over all degrees of freedom in a point
      int itmp = i*dof+k;
      m_ja[m_ia[itmp]-1] = itmp+1;    // put in column index of main diagonal
      // loop over all points surrounding point i
      int j = psup2[i]+1;
      for (int n=1; j<=psup2[i+1]; j++) {
        int e = psup1[j]*dof+k+1;
        // put in column index of an off-diagonal
	m_ja[m_ia[itmp]-1+(n++)] = e;
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
            if (m_ja[m_ia[i*dof+k]-1+e] > m_ja[m_ia[i*dof+k]+e]) {
              int itmp;
	      SWAP(m_ja[m_ia[i*dof+k]-1+e], m_ja[m_ia[i*dof+k]+e], itmp);
            }
          }
        }
      }
    }
  }

  // Free array for storing the nonzeros in each row
  rnz.reset();
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

  for (int n=0, j=m_ia[rmdof+i]-1; j<=m_ia[rmdof+i+1]-2; j++, n++)
    if (column*m_dof+i+1 == m_ja[j]) {
      idx = n;
      j = m_ia[rmdof+i+1]-1;
    }

  #ifdef _OPENMP
  #pragma omp atomic
  #endif
  m_a[m_ia[rmdof+i]-1+idx] += value;
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

  for (int n=0, j=m_ia[row]-1; j<=m_ia[row+1]-2; j++, n++)
    if (column+1 == m_ja[j]) {
      idx = n;
      j = m_ia[row+1]-1;
    }

  #ifdef _OPENMP
  #pragma omp atomic
  #endif
  m_a[m_ia[row]-1+idx] += value;
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
  
  for (int n=0, j=m_ia[rmdof+i]-1; j<=m_ia[rmdof+i+1]-2; j++, n++)
    if (column*m_dof+i+1 == m_ja[j]) {
      idx = n;
      j = m_ia[rmdof+i+1]-1;
    }

  m_a[m_ia[rmdof+i]-1+idx] = value;
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

  for (int n=0, j=m_ia[row]-1; j<=m_ia[row+1]-2; j++, n++)
    if (column+1 == m_ja[j]) {
      idx = n;
      j = m_ia[row+1]-1;
    }

  m_a[m_ia[row]-1+idx] = value;
}

real
SymCompRowMatrix::get(int row, int column, int i) const
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

  for (int n=0, j=m_ia[rmdof+i]-1; j<=m_ia[rmdof+i+1]-2; j++, n++)
    if (column*m_dof+i+1 == m_ja[j]) {
      idx = n;
      j = m_ia[rmdof+i+1]-1;
    }

  return m_a[m_ia[rmdof+i]-1+idx];
}

real
SymCompRowMatrix::get(int row, int column) const
//******************************************************************************
//  Get value from matrix from specified position using absolute indexing
//! \param[in]  row     block row
//! \param[in]  column  block column
//! \return             matrix value
//! \author  J. Bakosi
//******************************************************************************
{
  int idx = 0;

  for (int n=0, j=m_ia[row]-1; j<=m_ia[row+1]-2; j++, n++)
    if (column+1 == m_ja[j]) {
      idx = n;
      j = m_ia[row+1]-1;
    }

  return m_a[m_ia[row]-1+idx];
}

void
SymCompRowMatrix::echoAsStored(std::ostream& ofs) const
//******************************************************************************
//  Print out matrix entries as stored
//! \param[in]  ofs  output stream
//! \author  J. Bakosi
//******************************************************************************
{
  ofs << "name  = " << m_name << std::endl;
  ofs << "size  = " << m_size << std::endl;
  ofs << "rsize = " << m_rsize << std::endl;
  ofs << "dof   = " << m_dof << std::endl;
  ofs << "nnz   = " << m_nnz << std::endl;

  int i;
  ofs << "ia[] = { ";
  for (i=0; i<m_rsize; i++) ofs << m_ia[i];
  ofs << m_ia[i] << std::endl;

  ofs << "ja[] = { ";
  for (i=0; i<m_nnz-1; i++) ofs << m_ja[i];
  ofs << m_ja[i] << std::endl;

  ofs << "a[] = { ";
  for (i=0; i<m_nnz-1; i++) ofs << m_a[i];
  ofs << m_a[i] << std::endl;
}

void
SymCompRowMatrix::echoNonzeroStructure(std::ostream& ofs) const
//******************************************************************************
//  Print out nonzero structure of matrix
//! \param[in]  ofs  output stream
//! \author  J. Bakosi
//******************************************************************************
{
  for (int i=0; i<m_rsize; i++) {
    for (int j=1; j<m_ja[m_ia[i]-1]; j++)  // leading zeros
       ofs << ". ";

    for (int n=m_ia[i]-1; n<m_ia[i+1]-1; n++) {
      if (n>m_ia[i]-1) {  // zeros between nonzeros
        for (int j=m_ja[n-1]; j<m_ja[n]-1; j++)
           ofs << ". ";
      }
      ofs << "o ";  // nonzero
    }

    for (int j=m_ja[m_ia[i+1]-2]; j<m_rsize; j++)
       ofs << ". ";  // trailing zeros

    ofs << std::endl;
  }
}

void
SymCompRowMatrix::echoAsMatrix(std::ostream& ofs) const
//******************************************************************************
//  Print out matrix as a real matrix
//! \param[in]  ofs  output stream
//! \author  J. Bakosi
//******************************************************************************
{
  for (int i=0; i<m_rsize; i++) {
    for (int j=1; j<m_ja[m_ia[i]-1]; j++)
       ofs << "0\t";

    for (int n=m_ia[i]-1; n<m_ia[i+1]-1; n++) {
      if (n>m_ia[i]-1) {
        for (int j=m_ja[n-1]; j<m_ja[n]-1; j++)
          ofs << "0\t";
      }
      ofs << m_a[n] << "\t";
    }

    for (int j=m_ja[m_ia[i+1]-2]; j<m_rsize; j++)
      ofs << "0\t";

    ofs << std::endl;
  }
}

void
SymCompRowMatrix::echoAsMatlab(std::ostream& ofs) const
//******************************************************************************
//  Print out matrix as a maltab matrix
//! \param[in]  ofs  output stream
//! \author  J. Bakosi
//******************************************************************************
{
  ofs << "A = [ ";
  for (int i=0; i<m_rsize; i++) {
    for (int j=1; j<m_ja[m_ia[i]-1]; j++)
       ofs << "0 ";

    for (int n=m_ia[i]-1; n<m_ia[i+1]-1; n++) {
      if (n>m_ia[i]-1)
        for (int j=m_ja[n-1]; j<m_ja[n]-1; j++)
          ofs << "0 ";
      ofs << m_a[n] << " ";
    }

    for (int j=m_ja[m_ia[i+1]-2]; j<m_rsize; j++)
      ofs << "0 ";

    ofs << ";" << std::endl;
  }
  ofs << "]" << std::endl;
}
