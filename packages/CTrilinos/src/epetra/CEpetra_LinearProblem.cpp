
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

#include "CTrilinos_enums.h"
#include "CEpetra_LinearProblem.h"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "Epetra_LinearProblem.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_LinearProblem */
Table<Epetra_LinearProblem>& tableOfLinearProblems()
{
    static Table<Epetra_LinearProblem> loc_tableOfLinearProblems(CT_Epetra_LinearProblem_ID);
    return loc_tableOfLinearProblems;
}


} // namespace


//
// Definitions from CEpetra_LinearProblem.h
//


extern "C" {


CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_LinearProblem_Generalize ( 
  CT_Epetra_LinearProblem_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(id);
}

CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create (  )
{
    return CEpetra::storeNewLinearProblem(new Epetra_LinearProblem());
}

CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromMatrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_RowMatrix> A = CEpetra::getRowMatrix(AID);
    const Teuchos::RCP<Epetra_MultiVector> X = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> B = CEpetra::getMultiVector(BID);
    return CEpetra::storeNewLinearProblem(new Epetra_LinearProblem(
        A.getRawPtr(), X.getRawPtr(), B.getRawPtr()));
}

CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromOperator ( 
  CT_Epetra_Operator_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_Operator> A = CEpetra::getOperator(AID);
    const Teuchos::RCP<Epetra_MultiVector> X = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> B = CEpetra::getMultiVector(BID);
    return CEpetra::storeNewLinearProblem(new Epetra_LinearProblem(
        A.getRawPtr(), X.getRawPtr(), B.getRawPtr()));
}

CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Duplicate ( 
  CT_Epetra_LinearProblem_ID_t ProblemID )
{
    const Teuchos::RCP<const Epetra_LinearProblem> Problem = 
        CEpetra::getConstLinearProblem(ProblemID);
    return CEpetra::storeNewLinearProblem(new Epetra_LinearProblem(*Problem));
}

void Epetra_LinearProblem_Destroy ( 
  CT_Epetra_LinearProblem_ID_t * selfID )
{
    CEpetra::removeLinearProblem(selfID);
}

int Epetra_LinearProblem_CheckInput ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return CEpetra::getConstLinearProblem(selfID)->CheckInput();
}

void Epetra_LinearProblem_AssertSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    CEpetra::getLinearProblem(selfID)->AssertSymmetric();
}

void Epetra_LinearProblem_SetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_ProblemDifficultyLevel_E_t PDL )
{
    CEpetra::getLinearProblem(selfID)->SetPDL((ProblemDifficultyLevel) PDL);
}

void Epetra_LinearProblem_SetOperator_Matrix ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_RowMatrix_ID_t AID )
{
    const Teuchos::RCP<Epetra_RowMatrix> A = CEpetra::getRowMatrix(AID);
    CEpetra::getLinearProblem(selfID)->SetOperator(A.getRawPtr());
}

void Epetra_LinearProblem_SetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Operator_ID_t AID )
{
    const Teuchos::RCP<Epetra_Operator> A = CEpetra::getOperator(AID);
    CEpetra::getLinearProblem(selfID)->SetOperator(A.getRawPtr());
}

void Epetra_LinearProblem_SetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t XID )
{
    const Teuchos::RCP<Epetra_MultiVector> X = CEpetra::getMultiVector(XID);
    CEpetra::getLinearProblem(selfID)->SetLHS(X.getRawPtr());
}

void Epetra_LinearProblem_SetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_MultiVector> B = CEpetra::getMultiVector(BID);
    CEpetra::getLinearProblem(selfID)->SetRHS(B.getRawPtr());
}

int Epetra_LinearProblem_LeftScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID )
{
    const Teuchos::RCP<const Epetra_Vector> D = CEpetra::getConstVector(DID);
    return CEpetra::getLinearProblem(selfID)->LeftScale(*D);
}

int Epetra_LinearProblem_RightScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID )
{
    const Teuchos::RCP<const Epetra_Vector> D = CEpetra::getConstVector(DID);
    return CEpetra::getLinearProblem(selfID)->RightScale(*D);
}

CT_Epetra_Operator_ID_t Epetra_LinearProblem_GetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return CEpetra::storeOperator(CEpetra::getConstLinearProblem(
        selfID)->GetOperator());
}

CT_Epetra_RowMatrix_ID_t Epetra_LinearProblem_GetMatrix ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return CEpetra::storeRowMatrix(CEpetra::getConstLinearProblem(
        selfID)->GetMatrix());
}

CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return CEpetra::storeMultiVector(CEpetra::getConstLinearProblem(
        selfID)->GetLHS());
}

CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return CEpetra::storeMultiVector(CEpetra::getConstLinearProblem(
        selfID)->GetRHS());
}

CT_ProblemDifficultyLevel_E_t Epetra_LinearProblem_GetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return (CT_ProblemDifficultyLevel_E_t)( CEpetra::getConstLinearProblem(
        selfID)->GetPDL() );
}

boolean Epetra_LinearProblem_IsOperatorSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID )
{
    return ((CEpetra::getConstLinearProblem(
        selfID)->IsOperatorSymmetric()) ? TRUE : FALSE);
}


} // extern "C"


//
// Definitions from CEpetra_LinearProblem_Cpp.hpp
//


/* get Epetra_LinearProblem from non-const table using CT_Epetra_LinearProblem_ID */
const Teuchos::RCP<Epetra_LinearProblem>
CEpetra::getLinearProblem( CT_Epetra_LinearProblem_ID_t id )
{
    if (tableOfLinearProblems().isType(id.table))
        return tableOfLinearProblems().get<Epetra_LinearProblem>(
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_LinearProblem>(
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(id));
}

/* get Epetra_LinearProblem from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_LinearProblem>
CEpetra::getLinearProblem( CTrilinos_Universal_ID_t id )
{
    if (tableOfLinearProblems().isType(id.table))
        return tableOfLinearProblems().get<Epetra_LinearProblem>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_LinearProblem>(id);
}

/* get const Epetra_LinearProblem from either the const or non-const table
 * using CT_Epetra_LinearProblem_ID */
const Teuchos::RCP<const Epetra_LinearProblem>
CEpetra::getConstLinearProblem( CT_Epetra_LinearProblem_ID_t id )
{
    if (tableOfLinearProblems().isType(id.table))
        return tableOfLinearProblems().getConst<Epetra_LinearProblem>(
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_LinearProblem>(
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(id));
}

/* get const Epetra_LinearProblem from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_LinearProblem>
CEpetra::getConstLinearProblem( CTrilinos_Universal_ID_t id )
{
    if (tableOfLinearProblems().isType(id.table))
        return tableOfLinearProblems().getConst<Epetra_LinearProblem>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_LinearProblem>(id);
}

/* store Epetra_LinearProblem (owned) in non-const table */
CT_Epetra_LinearProblem_ID_t
CEpetra::storeNewLinearProblem( Epetra_LinearProblem *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(
        tableOfLinearProblems().store<Epetra_LinearProblem>(pobj, true));
}

/* store Epetra_LinearProblem in non-const table */
CT_Epetra_LinearProblem_ID_t
CEpetra::storeLinearProblem( Epetra_LinearProblem *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(
        tableOfLinearProblems().store<Epetra_LinearProblem>(pobj, false));
}

/* store const Epetra_LinearProblem in const table */
CT_Epetra_LinearProblem_ID_t
CEpetra::storeConstLinearProblem( const Epetra_LinearProblem *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(
        tableOfLinearProblems().store<Epetra_LinearProblem>(pobj, false));
}

/* remove Epetra_LinearProblem from table using CT_Epetra_LinearProblem_ID */
void
CEpetra::removeLinearProblem( CT_Epetra_LinearProblem_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(*id);
    if (tableOfLinearProblems().isType(aid.table))
        tableOfLinearProblems().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(aid);
}

/* remove Epetra_LinearProblem from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeLinearProblem( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfLinearProblems().isType(aid->table))
        tableOfLinearProblems().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_LinearProblem table */
void
CEpetra::purgeLinearProblem(  )
{
    tableOfLinearProblems().purge();
}

/* store Epetra_LinearProblem in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasLinearProblem( const Teuchos::RCP< Epetra_LinearProblem > & robj )
{
    return tableOfLinearProblems().alias(robj);
}

/* store const Epetra_LinearProblem in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstLinearProblem( const Teuchos::RCP< const Epetra_LinearProblem > & robj )
{
    return tableOfLinearProblems().alias(robj);
}



