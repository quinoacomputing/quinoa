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

#include "Ifpack_Euclid.h"
#if defined(HAVE_EUCLID) && defined(HAVE_MPI)

#include "Ifpack_Utils.h"
#include <algorithm>
#include "Epetra_MpiComm.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "getRow_dh.h"

using Teuchos::RCP;
using Teuchos::rcp;

Ifpack_Euclid::Ifpack_Euclid(Epetra_CrsMatrix* A):
  A_(rcp(A,false)),
  UseTranspose_(false),
  Condest_(-1),
  IsInitialized_(false),
  IsComputed_(false),
  Label_(),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  Time_(A_->Comm()),
  SetLevel_(1),
  SetBJ_(0),
  SetStats_(0),
  SetMem_(0),
  SetSparse_(0.0),
  SetRowScale_(0),
  SetILUT_(0.0)
{
  // Here we need to change the view of each row to have global indices. This is
  // because Euclid directly extracts a row view and expects global indices.
  for(int i = 0; i < A_->NumMyRows(); i++){
    int *indices;
    int len;
    A_->Graph().ExtractMyRowView(i, len, indices);
    for(int j = 0; j < len; j++){
      indices[j] = A_->GCID(indices[j]);
    }
  }
} //Constructor

//==============================================================================
void Ifpack_Euclid::Destroy(){
  // Destroy the euclid solver, only if it was setup
  if(IsComputed()){
    Euclid_dhDestroy(eu);
  }
  // Delete these euclid varaiables if they were created
  if(IsInitialized()){
    Parser_dhDestroy(parser_dh);
    parser_dh = NULL;
    TimeLog_dhDestroy(tlog_dh);
    tlog_dh = NULL;
    Mem_dhDestroy(mem_dh);
    mem_dh = NULL;
  }
  // Now that Euclid is done with the matrix, we change it back to having local indices.
  for(int i = 0; i < A_->NumMyRows(); i++){
    int *indices;
    int len;
    A_->Graph().ExtractMyRowView(i, len, indices);
    for(int j = 0; j < len; j++){
      indices[j] = A_->LCID(indices[j]);
    }
  }
} //Destroy()

//==============================================================================
int Ifpack_Euclid::Initialize(){
  //These are global variables in Euclid
  comm_dh = GetMpiComm();
  MPI_Comm_size(comm_dh, &np_dh);
  MPI_Comm_rank(comm_dh, &myid_dh);
  Time_.ResetStartTime();
  if(mem_dh == NULL){
    Mem_dhCreate(&mem_dh);
  }
  if (tlog_dh == NULL) {
    TimeLog_dhCreate(&tlog_dh);
  }

  if (parser_dh == NULL) {
    Parser_dhCreate(&parser_dh);
  }
  Parser_dhInit(parser_dh, 0, NULL);
  // Create the solver, this doesn't malloc anything yet, so it's only destroyed if Compute() is called.
  Euclid_dhCreate(&eu);
  IsInitialized_=true;
  NumInitialize_ = NumInitialize_ + 1;
  InitializeTime_ = InitializeTime_ + Time_.ElapsedTime();
  return 0;
} //Initialize()

//==============================================================================
int Ifpack_Euclid::SetParameters(Teuchos::ParameterList& list){
  List_ = list;
  SetLevel_ = list.get("SetLevel", (int)1);
  SetBJ_ = list.get("SetBJ", (int)0);
  SetStats_ = list.get("SetStats", (int)0);
  SetMem_ = list.get("SetMem", (int)0);
  SetSparse_ = list.get("SetSparse", (double)0.0);
  SetRowScale_ = list.get("SetRowScale", (int)0);
  SetILUT_ = list.get("SetILUT", (double)0.0);
  return 0;
} //SetParamters()

//==============================================================================
int Ifpack_Euclid::SetParameter(std::string name, int value){
  //Convert to lowercase (so it's case insensitive)
  locale loc;
  for(size_t i = 0; i < name.length(); i++){
    name[i] = (char) tolower(name[i], loc);
  }
  if(name.compare("setlevel") == 0){
    SetLevel_ = value;
  } else if(name.compare("setbj") == 0){
    SetBJ_ = value;
  } else if(name.compare("setstats") == 0){
    SetStats_ = value;
  } else if(name.compare("setmem") == 0){
    SetMem_ = value;
  } else if(name.compare("setrowscale") == 0){
    SetRowScale_ = value;
  } else {
    using std::cout;
    using std::endl;

    cout << "\nThe string " << name << " is not an available option." << endl;
    IFPACK_CHK_ERR(-1);
  }
  return 0;
} //SetParameter() (int)

//==============================================================================
int Ifpack_Euclid::SetParameter(std::string name, double value){
  //Convert to lowercase (so it's case insensitive)
  locale loc;
  for(size_t i; i < name.length(); i++){
    name[i] = (char) tolower(name[i], loc);
  }
  if(name.compare("setsparse") == 0){
    SetSparse_ = value;
  } else if(name.compare("setilut") == 0){
    SetILUT_ = value;
  } else {
    using std::cout;
    using std::endl;

    cout << "\nThe string " << name << " is not an available option." << endl;
    IFPACK_CHK_ERR(-1);
  }
  return 0;
} //SetParameter() (double)

//==============================================================================
int Ifpack_Euclid::Compute(){
  if(IsInitialized() == false){
    IFPACK_CHK_ERR(Initialize());
  }
  Time_.ResetStartTime();
  sprintf(Label_, "IFPACK_Euclid (level=%d, bj=%d, stats=%d, mem=%d, sparse=%f, rowscale=%d, ilut=%f)",
      SetLevel_, SetBJ_, SetStats_, SetMem_, SetSparse_, SetRowScale_, SetILUT_);
  // Set the parameters
  eu->level = SetLevel_;
  if(SetBJ_ != 0){
    strcpy("bj", eu->algo_par);
  }
  if(SetSparse_ != 0.0){
    eu->sparseTolA = SetSparse_;
  }
  if(SetRowScale_ != 0){
    eu->isScaled = true;
  }
  if(SetILUT_ != 0.0){
    eu->droptol = SetILUT_;
  }
  if(SetStats_ != 0 || SetMem_ != 0){
    eu->logging = true;
    Parser_dhInsert(parser_dh, "-eu_stats", "1");
  }
  // eu->A is the matrix as a void pointer, eu->m is local rows, eu->n is global rows
  eu->A = (void*) A_.get();
  eu->m = A_->NumMyRows();
  eu->n = A_->NumGlobalRows();
  Euclid_dhSetup(eu);
  IsComputed_ = true;
  NumCompute_ = NumCompute_ + 1;
  ComputeTime_ = ComputeTime_ + Time_.ElapsedTime();
  return 0;
} //Compute()

//==============================================================================
int Ifpack_Euclid::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
  if(IsComputed() == false){
    IFPACK_CHK_ERR(-1);
  }
  int NumVectors = X.NumVectors();
  if(NumVectors != Y.NumVectors()){
    IFPACK_CHK_ERR(-2);
  }
  Time_.ResetStartTime();
  // Loop through the vectors
  for(int vecNum = 0; vecNum < NumVectors; vecNum++){
    CallEuclid(X[vecNum], Y[vecNum]);
  }
  if(SetStats_ != 0){
    Euclid_dhPrintTestData(eu, stdout);
  }
  NumApplyInverse_ = NumApplyInverse_ + 1;
  ApplyInverseTime_ = ApplyInverseTime_ + Time_.ElapsedTime();
  return 0;
} //ApplyInverse()

//==============================================================================
int Ifpack_Euclid::CallEuclid(double *x, double *y) const{
  Euclid_dhApply(eu, x, y);
  return 0;
} //CallEuclid()

//==============================================================================
std::ostream& operator << (std::ostream& os, const Ifpack_Euclid& A){
  if (!A.Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_Euclid: " << A.Label () << endl << endl;
    os << "Using " << A.Comm().NumProc() << " processors." << endl;
    os << "Global number of rows            = " << A.Matrix().NumGlobalRows() << endl;
    os << "Global number of nonzeros        = " << A.Matrix().NumGlobalNonzeros() << endl;
    os << "Condition number estimate = " << A.Condest() << endl;
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << A.NumInitialize()
       << "  " << std::setw(15) << A.InitializeTime()
       << "              0.0              0.0" << endl;
    os << "Compute()       "   << std::setw(5) << A.NumCompute()
       << "  " << std::setw(15) << A.ComputeTime()
       << "  " << std::setw(15) << 1.0e-6 * A.ComputeFlops();
    if (A.ComputeTime() != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * A.ComputeFlops() / A.ComputeTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << A.NumApplyInverse()
       << "  " << std::setw(15) << A.ApplyInverseTime()
       << "  " << std::setw(15) << 1.0e-6 * A.ApplyInverseFlops();
    if (A.ApplyInverseTime() != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * A.ApplyInverseFlops() / A.ApplyInverseTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }
  return os;
} // <<

//==============================================================================
double Ifpack_Euclid::Condest(const Ifpack_CondestType CT,
                             const int MaxIters,
                             const double Tol,
                             Epetra_RowMatrix* Matrix_in){
  if (!IsComputed()) // cannot compute right now
    return(-1.0);
  return(Condest_);
} //Condest() - not implemented

#endif // HAVE_EUCLID && HAVE_MPI
