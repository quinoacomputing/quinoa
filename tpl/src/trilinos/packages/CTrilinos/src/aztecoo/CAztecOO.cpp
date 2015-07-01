
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


#include "CTrilinos_enums.h"
#include "CAztecOO.h"
#include "CAztecOO_Cpp.hpp"
#include "AztecOO.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type AztecOO */
Table<AztecOO>& tableOfAztecOOs()
{
    static Table<AztecOO> loc_tableOfAztecOOs(CT_AztecOO_ID);
    return loc_tableOfAztecOOs;
}


} // namespace


//
// Definitions from CAztecOO.h
//


extern "C" {


CT_AztecOO_ID_t AztecOO_Create_FromOperator ( 
  CT_Epetra_Operator_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_Operator> A = CEpetra::getOperator(AID);
    const Teuchos::RCP<Epetra_MultiVector> X = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> B = CEpetra::getMultiVector(BID);
    return CAztecOO::storeNewAztecOO(new AztecOO(A.getRawPtr(), X.getRawPtr(), 
        B.getRawPtr()));
}

CT_AztecOO_ID_t AztecOO_Create_FromRowMatrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_RowMatrix> A = CEpetra::getRowMatrix(AID);
    const Teuchos::RCP<Epetra_MultiVector> X = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> B = CEpetra::getMultiVector(BID);
    return CAztecOO::storeNewAztecOO(new AztecOO(A.getRawPtr(), X.getRawPtr(), 
        B.getRawPtr()));
}

CT_AztecOO_ID_t AztecOO_Create_FromLinearProblem ( 
  CT_Epetra_LinearProblem_ID_t LinearProblemID )
{
    const Teuchos::RCP<const Epetra_LinearProblem> LinearProblem = 
        CEpetra::getConstLinearProblem(LinearProblemID);
    return CAztecOO::storeNewAztecOO(new AztecOO(*LinearProblem));
}

CT_AztecOO_ID_t AztecOO_Create (  )
{
    return CAztecOO::storeNewAztecOO(new AztecOO());
}

CT_AztecOO_ID_t AztecOO_Duplicate ( CT_AztecOO_ID_t SolverID )
{
    const Teuchos::RCP<const AztecOO> Solver = CAztecOO::getConstAztecOO(
        SolverID);
    return CAztecOO::storeNewAztecOO(new AztecOO(*Solver));
}

void AztecOO_Destroy ( CT_AztecOO_ID_t * selfID )
{
    CAztecOO::removeAztecOO(selfID);
}

int AztecOO_SetProblem ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_LinearProblem_ID_t probID, 
  boolean call_SetPrecMatrix )
{
    const Teuchos::RCP<const Epetra_LinearProblem> prob = 
        CEpetra::getConstLinearProblem(probID);
    return CAztecOO::getAztecOO(selfID)->SetProblem(*prob, ((
        call_SetPrecMatrix) != FALSE ? true : false));
}

int AztecOO_SetUserOperator ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_Operator_ID_t UserOperatorID )
{
    const Teuchos::RCP<Epetra_Operator> UserOperator = CEpetra::getOperator(
        UserOperatorID);
    return CAztecOO::getAztecOO(selfID)->SetUserOperator(
        UserOperator.getRawPtr());
}

int AztecOO_SetUserMatrix ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t UserMatrixID, 
  boolean call_SetPrecMatrix )
{
    const Teuchos::RCP<Epetra_RowMatrix> UserMatrix = CEpetra::getRowMatrix(
        UserMatrixID);
    return CAztecOO::getAztecOO(selfID)->SetUserMatrix(UserMatrix.getRawPtr(), 
        ((call_SetPrecMatrix) != FALSE ? true : false));
}

int AztecOO_SetLHS ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_MultiVector_ID_t XID )
{
    const Teuchos::RCP<Epetra_MultiVector> X = CEpetra::getMultiVector(XID);
    return CAztecOO::getAztecOO(selfID)->SetLHS(X.getRawPtr());
}

int AztecOO_SetRHS ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_MultiVector> B = CEpetra::getMultiVector(BID);
    return CAztecOO::getAztecOO(selfID)->SetRHS(B.getRawPtr());
}

int AztecOO_SetPrecMatrix ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t PrecMatrixID )
{
    const Teuchos::RCP<Epetra_RowMatrix> PrecMatrix = CEpetra::getRowMatrix(
        PrecMatrixID);
    return CAztecOO::getAztecOO(selfID)->SetPrecMatrix(PrecMatrix.getRawPtr());
}

int AztecOO_SetPrecOperator ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_Operator_ID_t PrecOperatorID )
{
    const Teuchos::RCP<Epetra_Operator> PrecOperator = CEpetra::getOperator(
        PrecOperatorID);
    return CAztecOO::getAztecOO(selfID)->SetPrecOperator(
        PrecOperator.getRawPtr());
}

int AztecOO_ConstructPreconditioner ( 
  CT_AztecOO_ID_t selfID, double * condest )
{
    return CAztecOO::getAztecOO(selfID)->ConstructPreconditioner(*condest);
}

int AztecOO_DestroyPreconditioner ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getAztecOO(selfID)->DestroyPreconditioner();
}

double AztecOO_Condest ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getConstAztecOO(selfID)->Condest();
}

int AztecOO_CheckInput ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getConstAztecOO(selfID)->CheckInput();
}

CT_Epetra_LinearProblem_ID_t AztecOO_GetProblem ( 
  CT_AztecOO_ID_t selfID )
{
    return CEpetra::storeLinearProblem(CAztecOO::getConstAztecOO(
        selfID)->GetProblem());
}

CT_Epetra_Operator_ID_t AztecOO_GetUserOperator ( 
  CT_AztecOO_ID_t selfID )
{
    return CEpetra::storeOperator(CAztecOO::getConstAztecOO(
        selfID)->GetUserOperator());
}

CT_Epetra_RowMatrix_ID_t AztecOO_GetUserMatrix ( 
  CT_AztecOO_ID_t selfID )
{
    return CEpetra::storeRowMatrix(CAztecOO::getConstAztecOO(
        selfID)->GetUserMatrix());
}

CT_Epetra_Operator_ID_t AztecOO_GetPrecOperator ( 
  CT_AztecOO_ID_t selfID )
{
    return CEpetra::storeOperator(CAztecOO::getConstAztecOO(
        selfID)->GetPrecOperator());
}

CT_Epetra_RowMatrix_ID_t AztecOO_GetPrecMatrix ( 
  CT_AztecOO_ID_t selfID )
{
    return CEpetra::storeRowMatrix(CAztecOO::getConstAztecOO(
        selfID)->GetPrecMatrix());
}

CT_Epetra_MultiVector_ID_t AztecOO_GetLHS ( CT_AztecOO_ID_t selfID )
{
    return CEpetra::storeMultiVector(CAztecOO::getConstAztecOO(
        selfID)->GetLHS());
}

CT_Epetra_MultiVector_ID_t AztecOO_GetRHS ( CT_AztecOO_ID_t selfID )
{
    return CEpetra::storeMultiVector(CAztecOO::getConstAztecOO(
        selfID)->GetRHS());
}

void AztecOO_PrintLinearSystem ( 
  CT_AztecOO_ID_t selfID, const char * name )
{
    CAztecOO::getAztecOO(selfID)->PrintLinearSystem(name);
}

int AztecOO_SetParameters ( 
  CT_AztecOO_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t parameterlistID, 
  boolean cerr_warning_if_unused )
{
    const Teuchos::RCP<Teuchos::ParameterList> parameterlist = 
        CTeuchos::getParameterList(parameterlistID);
    return CAztecOO::getAztecOO(selfID)->SetParameters(*parameterlist, ((
        cerr_warning_if_unused) != FALSE ? true : false));
}

int AztecOO_SetAztecDefaults ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getAztecOO(selfID)->SetAztecDefaults();
}

int AztecOO_SetAztecOption ( 
  CT_AztecOO_ID_t selfID, int option, int value )
{
    return CAztecOO::getAztecOO(selfID)->SetAztecOption(option, value);
}

int AztecOO_GetAztecOption ( CT_AztecOO_ID_t selfID, int option )
{
    return CAztecOO::getAztecOO(selfID)->GetAztecOption(option);
}

int AztecOO_SetAztecParam ( 
  CT_AztecOO_ID_t selfID, int param, double value )
{
    return CAztecOO::getAztecOO(selfID)->SetAztecParam(param, value);
}

const int * AztecOO_GetAllAztecOptions ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getConstAztecOO(selfID)->GetAllAztecOptions();
}

const double * AztecOO_GetAllAztecParams ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getConstAztecOO(selfID)->GetAllAztecParams();
}

int AztecOO_SetAllAztecOptions ( 
  CT_AztecOO_ID_t selfID, const int * options )
{
    return CAztecOO::getAztecOO(selfID)->SetAllAztecOptions(options);
}

int AztecOO_SetAllAztecParams ( 
  CT_AztecOO_ID_t selfID, const double * params )
{
    return CAztecOO::getAztecOO(selfID)->SetAllAztecParams(params);
}

int AztecOO_Iterate_Current ( 
  CT_AztecOO_ID_t selfID, int MaxIters, double Tolerance )
{
    return CAztecOO::getAztecOO(selfID)->Iterate(MaxIters, Tolerance);
}

int AztecOO_Iterate ( 
  CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t AID, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t BID, 
  int MaxIters, double Tolerance )
{
    const Teuchos::RCP<Epetra_RowMatrix> A = CEpetra::getRowMatrix(AID);
    const Teuchos::RCP<Epetra_MultiVector> X = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> B = CEpetra::getMultiVector(BID);
    return CAztecOO::getAztecOO(selfID)->Iterate(A.getRawPtr(), X.getRawPtr(), 
        B.getRawPtr(), MaxIters, Tolerance);
}

int AztecOO_recursiveIterate ( 
  CT_AztecOO_ID_t selfID, int MaxIters, double Tolerance )
{
    return CAztecOO::getAztecOO(selfID)->recursiveIterate(MaxIters, Tolerance);
}

const double * AztecOO_GetAztecStatus ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getConstAztecOO(selfID)->GetAztecStatus();
}

int AztecOO_SetUseAdaptiveDefaultsTrue ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getAztecOO(selfID)->SetUseAdaptiveDefaultsTrue();
}

int AztecOO_SetAdaptiveParams ( 
  CT_AztecOO_ID_t selfID, int NumTrials, double * athresholds, 
  double * rthresholds, double condestThreshold, double maxFill, 
  int maxKspace )
{
    return CAztecOO::getAztecOO(selfID)->SetAdaptiveParams(NumTrials, 
        athresholds, rthresholds, condestThreshold, maxFill, maxKspace);
}

int AztecOO_AdaptiveIterate ( 
  CT_AztecOO_ID_t selfID, int MaxIters, int MaxSolveAttempts, 
  double Tolerance )
{
    return CAztecOO::getAztecOO(selfID)->AdaptiveIterate(MaxIters, 
        MaxSolveAttempts, Tolerance);
}

int AztecOO_NumIters ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getConstAztecOO(selfID)->NumIters();
}

double AztecOO_TrueResidual ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getConstAztecOO(selfID)->TrueResidual();
}

double AztecOO_ScaledResidual ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getConstAztecOO(selfID)->ScaledResidual();
}

double AztecOO_RecursiveResidual ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getConstAztecOO(selfID)->RecursiveResidual();
}

double AztecOO_SolveTime ( CT_AztecOO_ID_t selfID )
{
    return CAztecOO::getConstAztecOO(selfID)->SolveTime();
}

int AztecOO_GetAllAztecStatus ( 
  CT_AztecOO_ID_t selfID, double * status )
{
    return CAztecOO::getAztecOO(selfID)->GetAllAztecStatus(status);
}


} // extern "C"


//
// Definitions from CAztecOO_Cpp.hpp
//


/* get AztecOO from non-const table using CT_AztecOO_ID */
const Teuchos::RCP<AztecOO>
CAztecOO::getAztecOO( CT_AztecOO_ID_t id )
{
    return tableOfAztecOOs().get<AztecOO>(
        CTrilinos::abstractType<CT_AztecOO_ID_t>(id));
}

/* get AztecOO from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO>
CAztecOO::getAztecOO( CTrilinos_Universal_ID_t id )
{
    return tableOfAztecOOs().get<AztecOO>(id);
}

/* get const AztecOO from either the const or non-const table
 * using CT_AztecOO_ID */
const Teuchos::RCP<const AztecOO>
CAztecOO::getConstAztecOO( CT_AztecOO_ID_t id )
{
    return tableOfAztecOOs().getConst<AztecOO>(
        CTrilinos::abstractType<CT_AztecOO_ID_t>(id));
}

/* get const AztecOO from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO>
CAztecOO::getConstAztecOO( CTrilinos_Universal_ID_t id )
{
    return tableOfAztecOOs().getConst<AztecOO>(id);
}

/* store AztecOO (owned) in non-const table */
CT_AztecOO_ID_t
CAztecOO::storeNewAztecOO( AztecOO *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_ID_t>(
        tableOfAztecOOs().store<AztecOO>(pobj, true));
}

/* store AztecOO in non-const table */
CT_AztecOO_ID_t
CAztecOO::storeAztecOO( AztecOO *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_ID_t>(
        tableOfAztecOOs().store<AztecOO>(pobj, false));
}

/* store const AztecOO in const table */
CT_AztecOO_ID_t
CAztecOO::storeConstAztecOO( const AztecOO *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_ID_t>(
        tableOfAztecOOs().store<AztecOO>(pobj, false));
}

/* remove AztecOO from table using CT_AztecOO_ID */
void
CAztecOO::removeAztecOO( CT_AztecOO_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_ID_t>(*id);
    tableOfAztecOOs().remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_ID_t>(aid);
}

/* remove AztecOO from table using CTrilinos_Universal_ID_t */
void
CAztecOO::removeAztecOO( CTrilinos_Universal_ID_t *aid )
{
    tableOfAztecOOs().remove(aid);
}

/* purge AztecOO table */
void
CAztecOO::purgeAztecOO(  )
{
    tableOfAztecOOs().purge();
}



#endif /* HAVE_CTRILINOS_AZTECOO */


