/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_IC.h"
#include "Ifpack_IC_Utils.h"
#include "Ifpack_Condest.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Util.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
using Teuchos::RefCountPtr;
using Teuchos::rcp;

//==============================================================================
Ifpack_IC::Ifpack_IC(Epetra_RowMatrix* A) :
  A_(rcp(A,false)),
  Comm_(A->Comm()),
  UseTranspose_(false),
  Condest_(-1.0),
  Athresh_(0.0),
  Rthresh_(1.0),
  Droptol_(0.0),
  Lfil_(1.0),
  Aict_(0),
  Lict_(0),
  Ldiag_(0),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0),
  ApplyInverseTime_(0),
  Time_(Comm()),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0)
{
  Teuchos::ParameterList List;
  SetParameters(List);

}
//==============================================================================
Ifpack_IC::~Ifpack_IC()
{
  if (Lict_ != 0) {
    Ifpack_AIJMatrix * Lict = (Ifpack_AIJMatrix *) Lict_;
    delete [] Lict->ptr;
    delete [] Lict->col;
    delete [] Lict->val;
    delete Lict;
  }
  if (Aict_ != 0) {
    Ifpack_AIJMatrix * Aict = (Ifpack_AIJMatrix *) Aict_;
    delete Aict;
  }
  if (Ldiag_ != 0) delete [] Ldiag_;

  IsInitialized_ = false;
  IsComputed_ = false;
}

//==========================================================================
int Ifpack_IC::SetParameters(Teuchos::ParameterList& List)
{

  // Lfil_ = List.get("fact: level-of-fill",Lfil_); // Confusing parameter since Ifpack_IC can only do ICT not IC(k)
  Lfil_ = List.get("fact: ict level-of-fill", Lfil_); // Lfil_ is now the fill ratio as in ICT (1.0 means same #nonzeros as A)
  Athresh_ = List.get("fact: absolute threshold", Athresh_);
  Rthresh_ = List.get("fact: relative threshold", Rthresh_);
  Droptol_ = List.get("fact: drop tolerance", Droptol_);

  // set label
  sprintf(Label_, "IFPACK IC (fill=%f, drop=%f)",
          Lfil_, Droptol_);
  return(0);
}

//==========================================================================
int Ifpack_IC::Initialize()
{

  IsInitialized_ = false;

  // FIXME: construction of ILUK graph must be here

  ++NumInitialize_;
  IsInitialized_ = true;
  return(0);
}

//==========================================================================
int Ifpack_IC::ComputeSetup()
{
  // (re)allocate memory for ICT factors
  U_ = rcp(new Epetra_CrsMatrix(Copy, Matrix().RowMatrixRowMap(),
                                Matrix().RowMatrixRowMap(), 0));
  D_ = rcp(new Epetra_Vector(Matrix().RowMatrixRowMap()));

  if (U_.get() == 0 || D_.get() == 0)
    IFPACK_CHK_ERR(-5); // memory allocation error

  int ierr = 0;
  int i, j;
  int NumIn, NumU;
  bool DiagFound;
  int NumNonzeroDiags = 0;

  // Get Maximun Row length
  int MaxNumEntries = Matrix().MaxNumEntries();

  std::vector<int> InI(MaxNumEntries); // Allocate temp space
  std::vector<int> UI(MaxNumEntries);
  std::vector<double> InV(MaxNumEntries);
  std::vector<double> UV(MaxNumEntries);

  double *DV;
  ierr = D_->ExtractView(&DV); // Get view of diagonal

  // First we copy the user's matrix into diagonal vector and U, regardless of fill level

  int NumRows = Matrix().NumMyRows();

  for (i=0; i< NumRows; i++) {

    Matrix().ExtractMyRowCopy(i, MaxNumEntries, NumIn, &InV[0], &InI[0]); // Get Values and Indices

    // Split into L and U (we don't assume that indices are ordered).
    NumU = 0;
    DiagFound = false;

    for (j=0; j< NumIn; j++) {
      int k = InI[j];

      if (k==i) {
        DiagFound = true;
        // Store perturbed diagonal in Epetra_Vector D_
        DV[i] += Rthresh_ * InV[j] + EPETRA_SGN(InV[j]) * Athresh_;
      }

      else if (k < 0) return(-1); // Out of range
      else if (i<k && k<NumRows) {
        UI[NumU] = k;
        UV[NumU] = InV[j];
        NumU++;
      }
    }

    // Check in things for this row of L and U

    if (DiagFound) NumNonzeroDiags++;
    if (NumU) U_->InsertMyValues(i, NumU, &UV[0], &UI[0]);

  }

  U_->FillComplete(Matrix().OperatorDomainMap(),
                   Matrix().OperatorRangeMap());

  int ierr1 = 0;
  if (NumNonzeroDiags<U_->NumMyRows()) ierr1 = 1;
  Matrix().Comm().MaxAll(&ierr1, &ierr, 1);
  IFPACK_CHK_ERR(ierr);

  return(0);
}

//==========================================================================
int Ifpack_IC::Compute() {

  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());

  Time_.ResetStartTime();
  IsComputed_ = false;

  // copy matrix into L and U.
  IFPACK_CHK_ERR(ComputeSetup());

  int i;

  int m, n, nz, Nrhs, ldrhs, ldlhs;
  int * ptr=0, * ind;
  double * val, * rhs, * lhs;

  int ierr = Epetra_Util_ExtractHbData(U_.get(), 0, 0, m, n, nz, ptr, ind,
                                       val, Nrhs, rhs, ldrhs, lhs, ldlhs);
  if (ierr < 0)
    IFPACK_CHK_ERR(ierr);

  Ifpack_AIJMatrix * Aict;
  if (Aict_==0) {
    Aict = new Ifpack_AIJMatrix;
    Aict_ = (void *) Aict;
  }
  else Aict = (Ifpack_AIJMatrix *) Aict_;
  Ifpack_AIJMatrix * Lict;
  if (Lict_==0) {
    Lict = new Ifpack_AIJMatrix;
    Lict_ = (void *) Lict;
  }
  else {
    Lict = (Ifpack_AIJMatrix *) Lict_;
    Ifpack_AIJMatrix_dealloc( Lict );  // Delete storage, crout_ict will allocate it again.
  }
  if (Ldiag_ != 0) delete [] Ldiag_; // Delete storage, crout_ict will allocate it again.
  Aict->val = val;
  Aict->col = ind;
  Aict->ptr = ptr;
  double *DV;
  EPETRA_CHK_ERR(D_->ExtractView(&DV)); // Get view of diagonal

  // lfil is average number of nonzeros per row times fill ratio.
  // Currently crout_ict keeps a constant number of nonzeros per row.
  // TODO: Pass Lfil_ and make crout_ict keep variable #nonzeros per row.
  int lfil = (Lfil_ * nz)/n + .5;

  crout_ict(m, Aict, DV, Droptol_, lfil, Lict, &Ldiag_);

  // Get rid of unnecessary data
  delete [] ptr;

  // Create Epetra View of L from crout_ict
  U_ = rcp(new Epetra_CrsMatrix(View, A_->RowMatrixRowMap(), A_->RowMatrixRowMap(),0));
  D_ = rcp(new Epetra_Vector(View, A_->RowMatrixRowMap(), Ldiag_));

  ptr = Lict->ptr;
  ind = Lict->col;
  val = Lict->val;

  for (i=0; i< m; i++) {
    int NumEntries = ptr[i+1]-ptr[i];
    int * Indices = ind+ptr[i];
    double * Values = val+ptr[i];
    U_->InsertMyValues(i, NumEntries, Values, Indices);
  }

  U_->FillComplete(A_->OperatorDomainMap(), A_->OperatorRangeMap());
  D_->Reciprocal(*D_); // Put reciprocal of diagonal in this vector

#ifdef IFPACK_FLOPCOUNTERS
  double current_flops = 2 * nz; // Just an estimate
  double total_flops = 0;

  A_->Comm().SumAll(&current_flops, &total_flops, 1); // Get total madds across all PEs

  ComputeFlops_ += total_flops;
  // Now count the rest. NOTE: those flops are *global*
  ComputeFlops_ += (double) U_->NumGlobalNonzeros(); // Accounts for multiplier above
  ComputeFlops_ += (double) D_->GlobalLength(); // Accounts for reciprocal of diagonal
#endif
  ++NumCompute_;
  ComputeTime_ += Time_.ElapsedTime();


  IsComputed_ = true;

  return(0);

}

//=============================================================================
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_IC::ApplyInverse(const Epetra_MultiVector& X,
                            Epetra_MultiVector& Y) const
{

  if (!IsComputed())
    IFPACK_CHK_ERR(-3); // compute preconditioner first

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-2); // Return error: X and Y not the same size

  Time_.ResetStartTime();

  bool Upper = true;
  bool UnitDiagonal = true;

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  RefCountPtr< const Epetra_MultiVector > Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = rcp( new Epetra_MultiVector(X) );
  else
    Xcopy = rcp( &X, false );

  U_->Solve(Upper, true, UnitDiagonal, *Xcopy, Y);
  Y.Multiply(1.0, *D_, Y, 0.0); // y = D*y (D_ has inverse of diagonal)
  U_->Solve(Upper, false, UnitDiagonal, Y, Y); // Solve Uy = y

#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += 4.0 * U_->NumGlobalNonzeros();
  ApplyInverseFlops_ += D_->GlobalLength();
#endif

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_.ElapsedTime();

  return(0);

}

//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_IC::Apply(const Epetra_MultiVector& X,
                      Epetra_MultiVector& Y) const
{

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-2); // Return error: X and Y not the same size

  Epetra_MultiVector * X1 = (Epetra_MultiVector *) &X;
  Epetra_MultiVector * Y1 = (Epetra_MultiVector *) &Y;

  U_->Multiply(false, *X1, *Y1);
  Y1->Update(1.0, *X1, 1.0); // Y1 = Y1 + X1 (account for implicit unit diagonal)
  Y1->ReciprocalMultiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
  Epetra_MultiVector Y1temp(*Y1); // Need a temp copy of Y1
  U_->Multiply(true, Y1temp, *Y1);
  Y1->Update(1.0, Y1temp, 1.0); // (account for implicit unit diagonal)
  return(0);
}

//=============================================================================
double Ifpack_IC::Condest(const Ifpack_CondestType CT,
                          const int MaxIters, const double Tol,
                          Epetra_RowMatrix* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  if (Condest_ == -1.0)
    Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//=============================================================================
std::ostream&
Ifpack_IC::Print(std::ostream& os) const
{
  using std::endl;

  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_IC: " << Label() << endl << endl;
    os << "Level-of-fill      = " << LevelOfFill() << endl;
    os << "Absolute threshold = " << AbsoluteThreshold() << endl;
    os << "Relative threshold = " << RelativeThreshold() << endl;
    os << "Drop tolerance     = " << DropTolerance() << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows            = " << A_->NumGlobalRows64() << endl;
    if (IsComputed_) {
      os << "Number of nonzeros of H         = " << U_->NumGlobalNonzeros64() << endl;
      os << "nonzeros / rows                 = "
         << 1.0 * U_->NumGlobalNonzeros64() / U_->NumGlobalRows64() << endl;
    }
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize()
       << "  " << std::setw(15) << InitializeTime()
       << "               0.0            0.0" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute()
       << "  " << std::setw(15) << ComputeTime()
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops();
    if (ComputeTime() != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops() / ComputeTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse()
       << "  " << std::setw(15) << ApplyInverseTime()
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops();
    if (ApplyInverseTime() != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops() / ApplyInverseTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }


  return(os);
}
