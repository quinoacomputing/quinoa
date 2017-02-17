/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_SparsityFilter.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

//==============================================================================
Ifpack_SparsityFilter::Ifpack_SparsityFilter(const Teuchos::RefCountPtr<Epetra_RowMatrix>& Matrix,
					     int AllowedEntries, 
					     int AllowedBandwidth) :
  A_(Matrix),
  MaxNumEntries_(0),
  MaxNumEntriesA_(0),
  AllowedBandwidth_(AllowedBandwidth),
  AllowedEntries_(AllowedEntries),
  NumNonzeros_(0),
  NumRows_(0)
{
  // use this filter only on serial matrices
  if (A_->Comm().NumProc() != 1) {
    cerr << "Ifpack_SparsityFilter can be used with Comm().NumProc() == 1" << endl;
    cerr << "only. This class is a tool for Ifpack_AdditiveSchwarz," << endl;
    cerr << "and it is not meant to be used otherwise." << endl;
    exit(EXIT_FAILURE);
  }

  // only square serial matrices
  if ((A_->NumMyRows() != A_->NumMyCols()) ||
     (A_->NumMyRows() != A_->NumGlobalRows()))
    IFPACK_CHK_ERRV(-1);

  NumRows_ = A_->NumMyRows();
  MaxNumEntriesA_ = A_->MaxNumEntries();
  Indices_.resize(MaxNumEntriesA_);
  Values_.resize(MaxNumEntriesA_);

  // default value is to not consider bandwidth
  if (AllowedBandwidth_ == -1)
    AllowedBandwidth_ = NumRows_;
  
  // computes the number of nonzero elements per row in the 
  // dropped matrix. Stores this number in NumEntries_.
  // Also, computes the global number of nonzeros.
  vector<int>    Ind(MaxNumEntriesA_);
  vector<double> Val(MaxNumEntriesA_);

  NumEntries_.resize(NumRows_);
  for (int i = 0 ; i < NumRows_ ; ++i)
    NumEntries_[i] = MaxNumEntriesA_;

  for (int i = 0 ; i < A_->NumMyRows() ; ++i) {
    int Nnz;
    IFPACK_CHK_ERRV(ExtractMyRowCopy(i,MaxNumEntriesA_,Nnz,
				     &Val[0], &Ind[0]));

    NumEntries_[i] = Nnz;
    NumNonzeros_ += Nnz;
    if (Nnz > MaxNumEntries_)
      MaxNumEntries_ = Nnz;

  }
}

//==============================================================================
int Ifpack_SparsityFilter::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices) const
{
  if (Length < NumEntries_[MyRow])
    IFPACK_CHK_ERR(-1);

  int Nnz;
  IFPACK_CHK_ERR(A_->ExtractMyRowCopy(MyRow,MaxNumEntriesA_,Nnz,
				     &Values_[0],&Indices_[0]));

  double Threshold = 0.0;
    
  // this `if' is to define the cut-off value
  if (Nnz > AllowedEntries_) {
 
    vector<double> Values2(Nnz);
    int count = 0;
    for (int i = 0 ; i < Nnz ; ++i) {
      // skip diagonal entry (which is always inserted)
      if (Indices_[i] == MyRow)
	continue;
      // put absolute value
      Values2[count] = IFPACK_ABS(Values_[i]);
      count++;
    }

    for (int i = count ; i < Nnz ; ++i)
      Values2[i] = 0.0;

    // sort in descending order
    sort(Values2.rbegin(),Values2.rend());
    // get the cut-off value
    Threshold = Values2[AllowedEntries_ - 1];

  }

  // loop over all nonzero elements of row MyRow,
  // and drop elements below specified threshold.
  // Also, add value to zero diagonal
  NumEntries = 0;

  for (int i = 0 ; i < Nnz ; ++i) {

    if (IFPACK_ABS(Indices_[i] - MyRow) > AllowedBandwidth_)
      continue;

    if ((Indices_[i] != MyRow) && (IFPACK_ABS(Values_[i]) < Threshold))
      continue;

    Values[NumEntries] = Values_[i];
    Indices[NumEntries] = Indices_[i];

    NumEntries++;
    if (NumEntries > AllowedEntries_)
      break;
  }

  return(0);
}

//==============================================================================
int Ifpack_SparsityFilter::
ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  IFPACK_RETURN(A_->ExtractDiagonalCopy(Diagonal));
}

//==============================================================================
int Ifpack_SparsityFilter::
Multiply(bool TransA, const Epetra_MultiVector& X, 
	 Epetra_MultiVector& Y) const
{

  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors())
    IFPACK_CHK_ERR(-1);

  Y.PutScalar(0.0);

  vector<int> Indices(MaxNumEntries_);
  vector<double> Values(MaxNumEntries_);

  for (int i = 0 ; i < A_->NumMyRows() ; ++i) {

    int Nnz;
    ExtractMyRowCopy(i,MaxNumEntries_,Nnz,
		     &Values[0], &Indices[0]);
    if (!TransA) {
      // no transpose first
      for (int j = 0 ; j < NumVectors ; ++j) {
	for (int k = 0 ; k < Nnz ; ++k) {
	  Y[j][i] += Values[k] * X[j][Indices[k]];
	}
      }
    }
    else {
      // transpose here
      for (int j = 0 ; j < NumVectors ; ++j) {
	for (int k = 0 ; k < Nnz ; ++k) {
	  Y[j][Indices[k]] += Values[k] * X[j][i];
	}
      }
    }
  }

  return(0);
}

//==============================================================================
int Ifpack_SparsityFilter::
Solve(bool Upper, bool Trans, bool UnitDiagonal, 
      const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-98);
}

//==============================================================================
int Ifpack_SparsityFilter::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_RETURN(Multiply(UseTranspose(),X,Y));
}

//==============================================================================
int Ifpack_SparsityFilter::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-98); 
}
