/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_DenseContainer.hpp"
#include "Tpetra_RowMatrix.hpp"

//==============================================================================
int Ifpack2_DenseContainer::NumRows() const
{
  return(NumRows_);
}

//==============================================================================
int Ifpack2_DenseContainer::Initialize()
{
  
  IsInitialized_ = false;

  IFPACK2_CHK_ERR(LHS_.Reshape(NumRows_,NumVectors_));
  IFPACK2_CHK_ERR(RHS_.Reshape(NumRows_,NumVectors_));
  IFPACK2_CHK_ERR(ID_.Reshape(NumRows_,NumVectors_));
  IFPACK2_CHK_ERR(Matrix_.Reshape(NumRows_,NumRows_));

  // zero out matrix elements
  for (int i = 0 ; i < NumRows_ ; ++i)
    for (int j = 0 ; j < NumRows_ ; ++j)
      Matrix_(i,j) = 0.0;

  // zero out vector elements
  for (int i = 0 ; i < NumRows_ ; ++i)
    for (int j = 0 ; j < NumVectors_ ; ++j) {
      LHS_(i,j) = 0.0;
      RHS_(i,j) = 0.0;
    }

  // Set to -1 ID_'s
  for (int i = 0 ; i < NumRows_ ; ++i)
    ID_(i) = -1;  

  if (NumRows_ != 0) {
    IFPACK2_CHK_ERR(Solver_.SetMatrix(Matrix_));
    IFPACK2_CHK_ERR(Solver_.SetVectors(LHS_,RHS_));
  }

  IsInitialized_ = true;
  return(0);
  
}

//==============================================================================
double& Ifpack2_DenseContainer::LHS(const int i, const int Vector)
{
  return(LHS_.A()[Vector * NumRows_ + i]);
}
  
//==============================================================================
double& Ifpack2_DenseContainer::RHS(const int i, const int Vector)
{
  return(RHS_.A()[Vector * NumRows_ + i]);
}

//==============================================================================
int Ifpack2_DenseContainer::
SetMatrixElement(const int row, const int col, const double value)
{
  if (IsInitialized() == false)
    IFPACK2_CHK_ERR(Initialize());

  if ((row < 0) || (row >= NumRows())) {
    IFPACK2_CHK_ERR(-2); // not in range
  }

  if ((col < 0) || (col >= NumRows())) {
    IFPACK2_CHK_ERR(-2); // not in range
  }

  Matrix_(row, col) = value;

  return(0);

}

//==============================================================================
int Ifpack2_DenseContainer::ApplyInverse()
{

  if (!IsComputed()) {
    IFPACK2_CHK_ERR(-1);
  }
  
  if (NumRows_ != 0)
    IFPACK2_CHK_ERR(Solver_.Solve());

  ApplyInverseFlops_ += 2.0 * NumVectors_ * NumRows_ * NumRows_;
  return(0);
}

//==============================================================================
int& Ifpack2_DenseContainer::ID(const int i)
{
  return(ID_[i]);
}

//==============================================================================
// FIXME: optimize performances of this guy...
int Ifpack2_DenseContainer::Extract(const Tpetra_RowMatrix& Matrix_in)
{

  for (int j = 0 ; j < NumRows_ ; ++j) {
    // be sure that the user has set all the ID's
    if (ID(j) == -1)
      IFPACK2_CHK_ERR(-2);
    // be sure that all are local indices
    if (ID(j) > Matrix_in.NumMyRows())
      IFPACK2_CHK_ERR(-2);
  }

  // allocate storage to extract matrix rows.
  int Length = Matrix_in.MaxNumEntries();
  vector<double> Values;
  Values.resize(Length);
  vector<int> Indices;
  Indices.resize(Length);

  for (int j = 0 ; j < NumRows_ ; ++j) {

    int LRID = ID(j);

    int NumEntries;

    int ierr = 
      Matrix_in.ExtractMyRowCopy(LRID, Length, NumEntries, 
			      &Values[0], &Indices[0]);
    IFPACK2_CHK_ERR(ierr);

    for (int k = 0 ; k < NumEntries ; ++k) {

      int LCID = Indices[k];

      // skip off-processor elements
      if (LCID >= Matrix_in.NumMyRows()) 
	continue;

      // for local column IDs, look for each ID in the list
      // of columns hosted by this object
      // FIXME: use STL
      int jj = -1;
      for (int kk = 0 ; kk < NumRows_ ; ++kk)
	if (ID(kk) == LCID)
	  jj = kk;

      if (jj != -1)
	SetMatrixElement(j,jj,Values[k]);

    }
  }

  return(0);
}

//==============================================================================
int Ifpack2_DenseContainer::Compute(const Tpetra_RowMatrix& Matrix_in)
{
  IsComputed_ = false;
  if (IsInitialized() == false) {
    IFPACK2_CHK_ERR(Initialize());
  }

  if (KeepNonFactoredMatrix_)
    NonFactoredMatrix_ = Matrix_;

  // extract local rows and columns
  IFPACK2_CHK_ERR(Extract(Matrix_in));

  if (KeepNonFactoredMatrix_)
    NonFactoredMatrix_ = Matrix_;

  // factorize the matrix using LAPACK
  if (NumRows_ != 0)
    IFPACK2_CHK_ERR(Solver_.Factor());

  Label_ = "Ifpack2_DenseContainer";

  // not sure of count
  ComputeFlops_ += 4.0 * NumRows_ * NumRows_ * NumRows_ / 3;
  IsComputed_ = true;

  return(0);
}

//==============================================================================
int Ifpack2_DenseContainer::Apply()
{
  if (IsComputed() == false)
    IFPACK2_CHK_ERR(-3);

  if (KeepNonFactoredMatrix_) {
    IFPACK2_CHK_ERR(RHS_.Multiply('N','N', 1.0,NonFactoredMatrix_,LHS_,0.0));
  }
  else
    IFPACK2_CHK_ERR(RHS_.Multiply('N','N', 1.0,Matrix_,LHS_,0.0));

  ApplyFlops_ += 2 * NumRows_ * NumRows_;
  return(0);
}

//==============================================================================
ostream& Ifpack2_DenseContainer::Print(ostream & os) const
{
    os << "================================================================================" << endl;
  os << "Ifpack2_DenseContainer" << endl;
  os << "Number of rows          = " << NumRows() << endl;
  os << "Number of vectors       = " << NumVectors() << endl;
  os << "IsInitialized()         = " << IsInitialized() << endl;
  os << "IsComputed()            = " << IsComputed() << endl;
  os << "Flops in Compute()      = " << ComputeFlops() << endl; 
  os << "Flops in ApplyInverse() = " << ApplyInverseFlops() << endl; 
  os << "================================================================================" << endl;
  os << endl;

  return(os);
}
