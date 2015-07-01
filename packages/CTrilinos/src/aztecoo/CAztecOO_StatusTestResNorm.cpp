
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
#include "CAztecOO_StatusTestResNorm.h"
#include "CAztecOO_StatusTestResNorm_Cpp.hpp"
#include "AztecOO_StatusTestResNorm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type AztecOO_StatusTestResNorm */
Table<AztecOO_StatusTestResNorm>& tableOfStatusTestResNorms()
{
    static Table<AztecOO_StatusTestResNorm> loc_tableOfStatusTestResNorms(CT_AztecOO_StatusTestResNorm_ID);
    return loc_tableOfStatusTestResNorms;
}


} // namespace


//
// Definitions from CAztecOO_StatusTestResNorm.h
//


extern "C" {


CT_AztecOO_StatusTestResNorm_ID_t AztecOO_StatusTestResNorm_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(id);
}

CTrilinos_Universal_ID_t AztecOO_StatusTestResNorm_Generalize ( 
  CT_AztecOO_StatusTestResNorm_ID_t id )
{
    return CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(id);
}

CT_AztecOO_StatusTestResNorm_ID_t AztecOO_StatusTestResNorm_Create ( 
  CT_Epetra_Operator_ID_t OperatorID, CT_Epetra_Vector_ID_t LHSID, 
  CT_Epetra_Vector_ID_t RHSID, double Tolerance )
{
    const Teuchos::RCP<const Epetra_Operator> Operator = 
        CEpetra::getConstOperator(OperatorID);
    const Teuchos::RCP<const Epetra_Vector> LHS = CEpetra::getConstVector(
        LHSID);
    const Teuchos::RCP<const Epetra_Vector> RHS = CEpetra::getConstVector(
        RHSID);
    return CAztecOO::storeNewStatusTestResNorm(new AztecOO_StatusTestResNorm(
        *Operator, *LHS, *RHS, Tolerance));
}

void AztecOO_StatusTestResNorm_Destroy ( 
  CT_AztecOO_StatusTestResNorm_ID_t * selfID )
{
    CAztecOO::removeStatusTestResNorm(selfID);
}

int AztecOO_StatusTestResNorm_DefineResForm ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID, 
  CT_ResType_E_t TypeOfResidual, CT_NormType_E_t TypeOfNorm, 
  CT_Epetra_Vector_ID_t WeightsID )
{
    const Teuchos::RCP<Epetra_Vector> Weights = CEpetra::getVector(WeightsID);
    return CAztecOO::getStatusTestResNorm(selfID)->DefineResForm(
        (AztecOO_StatusTestResNorm::ResType) TypeOfResidual, 
        (AztecOO_StatusTestResNorm::NormType) TypeOfNorm, 
        Weights.getRawPtr());
}

int AztecOO_StatusTestResNorm_DefineScaleForm ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID, 
  CT_ScaleType_E_t TypeOfScaling, CT_NormType_E_t TypeOfNorm, 
  CT_Epetra_Vector_ID_t WeightsID, double ScaleValue )
{
    const Teuchos::RCP<Epetra_Vector> Weights = CEpetra::getVector(WeightsID);
    return CAztecOO::getStatusTestResNorm(selfID)->DefineScaleForm(
        (AztecOO_StatusTestResNorm::ScaleType) TypeOfScaling, 
        (AztecOO_StatusTestResNorm::NormType) TypeOfNorm, Weights.getRawPtr(), 
        ScaleValue);
}

int AztecOO_StatusTestResNorm_ResetTolerance ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID, double Tolerance )
{
    return CAztecOO::getStatusTestResNorm(selfID)->ResetTolerance(Tolerance);
}

int AztecOO_StatusTestResNorm_SetMaxNumExtraIterations ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID, 
  int maxNumExtraIterations )
{
    return CAztecOO::getStatusTestResNorm(selfID)->SetMaxNumExtraIterations(
        maxNumExtraIterations);
}

int AztecOO_StatusTestResNorm_GetMaxNumExtraIterations ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return CAztecOO::getStatusTestResNorm(selfID)->GetMaxNumExtraIterations();
}

boolean AztecOO_StatusTestResNorm_ResidualVectorRequired ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return ((CAztecOO::getConstStatusTestResNorm(
        selfID)->ResidualVectorRequired()) ? TRUE : FALSE);
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestResNorm_CheckStatus ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID, int CurrentIter, 
  CT_Epetra_MultiVector_ID_t CurrentResVectorID, 
  double CurrentResNormEst, boolean SolutionUpdated )
{
    const Teuchos::RCP<Epetra_MultiVector> CurrentResVector = 
        CEpetra::getMultiVector(CurrentResVectorID);
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getStatusTestResNorm(
        selfID)->CheckStatus(CurrentIter, CurrentResVector.getRawPtr(), 
        CurrentResNormEst, ((SolutionUpdated) != FALSE ? true : false)) );
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestResNorm_GetStatus ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getConstStatusTestResNorm(
        selfID)->GetStatus() );
}

void AztecOO_StatusTestResNorm_ResetStatus ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    CAztecOO::getStatusTestResNorm(selfID)->ResetStatus();
}

double AztecOO_StatusTestResNorm_GetTolerance ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return CAztecOO::getConstStatusTestResNorm(selfID)->GetTolerance();
}

double AztecOO_StatusTestResNorm_GetTestValue ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return CAztecOO::getConstStatusTestResNorm(selfID)->GetTestValue();
}

double AztecOO_StatusTestResNorm_GetResNormValue ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return CAztecOO::getConstStatusTestResNorm(selfID)->GetResNormValue();
}

double AztecOO_StatusTestResNorm_GetScaledNormValue ( 
  CT_AztecOO_StatusTestResNorm_ID_t selfID )
{
    return CAztecOO::getConstStatusTestResNorm(selfID)->GetScaledNormValue();
}


} // extern "C"


//
// Definitions from CAztecOO_StatusTestResNorm_Cpp.hpp
//


/* get AztecOO_StatusTestResNorm from non-const table using CT_AztecOO_StatusTestResNorm_ID */
const Teuchos::RCP<AztecOO_StatusTestResNorm>
CAztecOO::getStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t id )
{
    if (tableOfStatusTestResNorms().isType(id.table))
        return tableOfStatusTestResNorms().get<AztecOO_StatusTestResNorm>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<AztecOO_StatusTestResNorm>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(id));
}

/* get AztecOO_StatusTestResNorm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTestResNorm>
CAztecOO::getStatusTestResNorm( CTrilinos_Universal_ID_t id )
{
    if (tableOfStatusTestResNorms().isType(id.table))
        return tableOfStatusTestResNorms().get<AztecOO_StatusTestResNorm>(id);
    else
        return CTrilinos::TableRepos::get<AztecOO_StatusTestResNorm>(id);
}

/* get const AztecOO_StatusTestResNorm from either the const or non-const table
 * using CT_AztecOO_StatusTestResNorm_ID */
const Teuchos::RCP<const AztecOO_StatusTestResNorm>
CAztecOO::getConstStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t id )
{
    if (tableOfStatusTestResNorms().isType(id.table))
        return tableOfStatusTestResNorms().getConst<AztecOO_StatusTestResNorm>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<AztecOO_StatusTestResNorm>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(id));
}

/* get const AztecOO_StatusTestResNorm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTestResNorm>
CAztecOO::getConstStatusTestResNorm( CTrilinos_Universal_ID_t id )
{
    if (tableOfStatusTestResNorms().isType(id.table))
        return tableOfStatusTestResNorms().getConst<AztecOO_StatusTestResNorm>(id);
    else
        return CTrilinos::TableRepos::getConst<AztecOO_StatusTestResNorm>(id);
}

/* store AztecOO_StatusTestResNorm (owned) in non-const table */
CT_AztecOO_StatusTestResNorm_ID_t
CAztecOO::storeNewStatusTestResNorm( AztecOO_StatusTestResNorm *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(
        tableOfStatusTestResNorms().store<AztecOO_StatusTestResNorm>(pobj, true));
}

/* store AztecOO_StatusTestResNorm in non-const table */
CT_AztecOO_StatusTestResNorm_ID_t
CAztecOO::storeStatusTestResNorm( AztecOO_StatusTestResNorm *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(
        tableOfStatusTestResNorms().store<AztecOO_StatusTestResNorm>(pobj, false));
}

/* store const AztecOO_StatusTestResNorm in const table */
CT_AztecOO_StatusTestResNorm_ID_t
CAztecOO::storeConstStatusTestResNorm( const AztecOO_StatusTestResNorm *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(
        tableOfStatusTestResNorms().store<AztecOO_StatusTestResNorm>(pobj, false));
}

/* remove AztecOO_StatusTestResNorm from table using CT_AztecOO_StatusTestResNorm_ID */
void
CAztecOO::removeStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(*id);
    if (tableOfStatusTestResNorms().isType(aid.table))
        tableOfStatusTestResNorms().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(aid);
}

/* remove AztecOO_StatusTestResNorm from table using CTrilinos_Universal_ID_t */
void
CAztecOO::removeStatusTestResNorm( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfStatusTestResNorms().isType(aid->table))
        tableOfStatusTestResNorms().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge AztecOO_StatusTestResNorm table */
void
CAztecOO::purgeStatusTestResNorm(  )
{
    tableOfStatusTestResNorms().purge();
}

/* store AztecOO_StatusTestResNorm in non-const table */
CTrilinos_Universal_ID_t
CAztecOO::aliasStatusTestResNorm( const Teuchos::RCP< AztecOO_StatusTestResNorm > & robj )
{
    return tableOfStatusTestResNorms().alias(robj);
}

/* store const AztecOO_StatusTestResNorm in const table */
CTrilinos_Universal_ID_t
CAztecOO::aliasConstStatusTestResNorm( const Teuchos::RCP< const AztecOO_StatusTestResNorm > & robj )
{
    return tableOfStatusTestResNorms().alias(robj);
}



#endif /* HAVE_CTRILINOS_AZTECOO */


