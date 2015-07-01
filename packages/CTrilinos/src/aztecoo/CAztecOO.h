#ifndef CAZTECOO_H
#define CAZTECOO_H

/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"


#ifdef HAVE_CTRILINOS_AZTECOO



/*! @file CAztecOO.h
 * @brief Wrappers for AztecOO */

/* True C header file! */

#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name AztecOO constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   AztecOO::AztecOO(Epetra_Operator * A, Epetra_MultiVector * X, Epetra_MultiVector * B)
*/
CT_AztecOO_ID_t AztecOO_Create_FromOperator ( 
  CT_Epetra_Operator_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   AztecOO::AztecOO(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B)
*/
CT_AztecOO_ID_t AztecOO_Create_FromRowMatrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   AztecOO::AztecOO(const Epetra_LinearProblem& LinearProblem)
*/
CT_AztecOO_ID_t AztecOO_Create_FromLinearProblem ( 
  CT_Epetra_LinearProblem_ID_t LinearProblemID );

/*! @brief Wrapper for 
   AztecOO::AztecOO()
*/
CT_AztecOO_ID_t AztecOO_Create (  );

/*! @brief Wrapper for 
   AztecOO::AztecOO(const AztecOO& Solver)
*/
CT_AztecOO_ID_t AztecOO_Duplicate ( CT_AztecOO_ID_t SolverID );

/*@}*/

/*! @name AztecOO destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual AztecOO::~AztecOO(void)
*/
void AztecOO_Destroy ( CT_AztecOO_ID_t * selfID );

/*@}*/

/*! @name AztecOO member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int AztecOO::SetProblem(const Epetra_LinearProblem& prob, bool call_SetPrecMatrix=false)
*/
int AztecOO_SetProblem ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_LinearProblem_ID_t probID, 
  boolean call_SetPrecMatrix );

/*! @brief Wrapper for 
   int AztecOO::SetUserOperator(Epetra_Operator * UserOperator)
*/
int AztecOO_SetUserOperator ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_Operator_ID_t UserOperatorID );

/*! @brief Wrapper for 
   int AztecOO::SetUserMatrix(Epetra_RowMatrix * UserMatrix, bool call_SetPrecMatrix=false)
*/
int AztecOO_SetUserMatrix ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t UserMatrixID, 
  boolean call_SetPrecMatrix );

/*! @brief Wrapper for 
   int AztecOO::SetLHS(Epetra_MultiVector * X)
*/
int AztecOO_SetLHS ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_MultiVector_ID_t XID );

/*! @brief Wrapper for 
   int AztecOO::SetRHS(Epetra_MultiVector * B)
*/
int AztecOO_SetRHS ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   int AztecOO::SetPrecMatrix(Epetra_RowMatrix * PrecMatrix)
*/
int AztecOO_SetPrecMatrix ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t PrecMatrixID );

/*! @brief Wrapper for 
   int AztecOO::SetPrecOperator(Epetra_Operator * PrecOperator)
*/
int AztecOO_SetPrecOperator ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_Operator_ID_t PrecOperatorID );

/*! @brief Wrapper for 
   int AztecOO::ConstructPreconditioner(double & condest)
*/
int AztecOO_ConstructPreconditioner ( 
  CT_AztecOO_ID_t selfID, double * condest );

/*! @brief Wrapper for 
   int AztecOO::DestroyPreconditioner()
*/
int AztecOO_DestroyPreconditioner ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   double AztecOO::Condest() const
*/
double AztecOO_Condest ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   int AztecOO::CheckInput() const
*/
int AztecOO_CheckInput ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_LinearProblem * AztecOO::GetProblem() const
*/
CT_Epetra_LinearProblem_ID_t AztecOO_GetProblem ( 
  CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_Operator * AztecOO::GetUserOperator() const
*/
CT_Epetra_Operator_ID_t AztecOO_GetUserOperator ( 
  CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_RowMatrix * AztecOO::GetUserMatrix() const
*/
CT_Epetra_RowMatrix_ID_t AztecOO_GetUserMatrix ( 
  CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_Operator * AztecOO::GetPrecOperator() const
*/
CT_Epetra_Operator_ID_t AztecOO_GetPrecOperator ( 
  CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_RowMatrix * AztecOO::GetPrecMatrix() const
*/
CT_Epetra_RowMatrix_ID_t AztecOO_GetPrecMatrix ( 
  CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_MultiVector * AztecOO::GetLHS() const
*/
CT_Epetra_MultiVector_ID_t AztecOO_GetLHS ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_MultiVector * AztecOO::GetRHS() const
*/
CT_Epetra_MultiVector_ID_t AztecOO_GetRHS ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   void AztecOO::PrintLinearSystem(const char* name)
*/
void AztecOO_PrintLinearSystem ( 
  CT_AztecOO_ID_t selfID, const char * name );

/*! @brief Wrapper for 
   int AztecOO::SetParameters(Teuchos::ParameterList& parameterlist, bool cerr_warning_if_unused=false)
*/
int AztecOO_SetParameters ( 
  CT_AztecOO_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t parameterlistID, 
  boolean cerr_warning_if_unused );

/*! @brief Wrapper for 
   int AztecOO::SetAztecDefaults()
*/
int AztecOO_SetAztecDefaults ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   int AztecOO::SetAztecOption(int option, int value)
*/
int AztecOO_SetAztecOption ( 
  CT_AztecOO_ID_t selfID, int option, int value );

/*! @brief Wrapper for 
   int AztecOO::GetAztecOption(int option)
*/
int AztecOO_GetAztecOption ( CT_AztecOO_ID_t selfID, int option );

/*! @brief Wrapper for 
   int AztecOO::SetAztecParam(int param, double value)
*/
int AztecOO_SetAztecParam ( 
  CT_AztecOO_ID_t selfID, int param, double value );

/*! @brief Wrapper for 
   const int* AztecOO::GetAllAztecOptions() const
*/
const int * AztecOO_GetAllAztecOptions ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   const double* AztecOO::GetAllAztecParams() const
*/
const double * AztecOO_GetAllAztecParams ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   int AztecOO::SetAllAztecOptions(const int * options)
*/
int AztecOO_SetAllAztecOptions ( 
  CT_AztecOO_ID_t selfID, const int * options );

/*! @brief Wrapper for 
   int AztecOO::SetAllAztecParams(const double * params)
*/
int AztecOO_SetAllAztecParams ( 
  CT_AztecOO_ID_t selfID, const double * params );

/*! @brief Wrapper for 
   int AztecOO::Iterate(int MaxIters, double Tolerance)
*/
int AztecOO_Iterate_Current ( 
  CT_AztecOO_ID_t selfID, int MaxIters, double Tolerance );

/*! @brief Wrapper for 
   int AztecOO::Iterate(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B, int MaxIters, double Tolerance)
*/
int AztecOO_Iterate ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t AID, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t BID, 
  int MaxIters, double Tolerance );

/*! @brief Wrapper for 
   int AztecOO::recursiveIterate(int MaxIters, double Tolerance)
*/
int AztecOO_recursiveIterate ( 
  CT_AztecOO_ID_t selfID, int MaxIters, double Tolerance );

/*! @brief Wrapper for 
   const double *AztecOO::GetAztecStatus() const
*/
const double * AztecOO_GetAztecStatus ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   int AztecOO::SetUseAdaptiveDefaultsTrue()
*/
int AztecOO_SetUseAdaptiveDefaultsTrue ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   int AztecOO::SetAdaptiveParams(int NumTrials, double * athresholds, double * rthresholds, double condestThreshold, double maxFill, int maxKspace)
*/
int AztecOO_SetAdaptiveParams ( 
  CT_AztecOO_ID_t selfID, int NumTrials, double * athresholds, 
  double * rthresholds, double condestThreshold, double maxFill, 
  int maxKspace );

/*! @brief Wrapper for 
   int AztecOO::AdaptiveIterate(int MaxIters, int MaxSolveAttempts, double Tolerance)
*/
int AztecOO_AdaptiveIterate ( 
  CT_AztecOO_ID_t selfID, int MaxIters, int MaxSolveAttempts, 
  double Tolerance );

/*! @brief Wrapper for 
   int AztecOO::NumIters() const
*/
int AztecOO_NumIters ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   double AztecOO::TrueResidual() const
*/
double AztecOO_TrueResidual ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   double AztecOO::ScaledResidual() const
*/
double AztecOO_ScaledResidual ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   double AztecOO::RecursiveResidual() const
*/
double AztecOO_RecursiveResidual ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   double AztecOO::SolveTime() const
*/
double AztecOO_SolveTime ( CT_AztecOO_ID_t selfID );

/*! @brief Wrapper for 
   int AztecOO::GetAllAztecStatus(double * status)
*/
int AztecOO_GetAllAztecStatus ( 
  CT_AztecOO_ID_t selfID, double * status );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* HAVE_CTRILINOS_AZTECOO */
#endif /* CAZTECOO_H */

