
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
#include "CEpetra_Operator.h"
#include "CEpetra_Operator_Cpp.hpp"
#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Operator */
Table<Epetra_Operator>& tableOfOperators()
{
    static Table<Epetra_Operator> loc_tableOfOperators(CT_Epetra_Operator_ID);
    return loc_tableOfOperators;
}


} // namespace


//
// Definitions from CEpetra_Operator.h
//


extern "C" {


CT_Epetra_Operator_ID_t Epetra_Operator_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Operator_Generalize ( 
  CT_Epetra_Operator_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id);
}

void Epetra_Operator_Destroy ( CT_Epetra_Operator_ID_t * selfID )
{
    CEpetra::removeOperator(selfID);
}

int Epetra_Operator_SetUseTranspose ( 
  CT_Epetra_Operator_ID_t selfID, boolean UseTranspose )
{
    return CEpetra::getOperator(selfID)->SetUseTranspose(((UseTranspose) != 
        FALSE ? true : false));
}

int Epetra_Operator_Apply ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstOperator(selfID)->Apply(*X, *Y);
}

int Epetra_Operator_ApplyInverse ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstOperator(selfID)->ApplyInverse(*X, *Y);
}

double Epetra_Operator_NormInf ( CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::getConstOperator(selfID)->NormInf();
}

const char * Epetra_Operator_Label ( CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::getConstOperator(selfID)->Label();
}

boolean Epetra_Operator_UseTranspose ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return ((CEpetra::getConstOperator(
        selfID)->UseTranspose()) ? TRUE : FALSE);
}

boolean Epetra_Operator_HasNormInf ( CT_Epetra_Operator_ID_t selfID )
{
    return ((CEpetra::getConstOperator(selfID)->HasNormInf()) ? TRUE : FALSE);
}

CT_Epetra_Comm_ID_t Epetra_Operator_Comm ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::storeConstComm(&( CEpetra::getConstOperator(
        selfID)->Comm() ));
}

CT_Epetra_Map_ID_t Epetra_Operator_OperatorDomainMap ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstOperator(
        selfID)->OperatorDomainMap() ));
}

CT_Epetra_Map_ID_t Epetra_Operator_OperatorRangeMap ( 
  CT_Epetra_Operator_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstOperator(
        selfID)->OperatorRangeMap() ));
}


} // extern "C"


//
// Definitions from CEpetra_Operator_Cpp.hpp
//


/* get Epetra_Operator from non-const table using CT_Epetra_Operator_ID */
const Teuchos::RCP<Epetra_Operator>
CEpetra::getOperator( CT_Epetra_Operator_ID_t id )
{
    if (tableOfOperators().isType(id.table))
        return tableOfOperators().get<Epetra_Operator>(
        CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_Operator>(
        CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id));
}

/* get Epetra_Operator from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Operator>
CEpetra::getOperator( CTrilinos_Universal_ID_t id )
{
    if (tableOfOperators().isType(id.table))
        return tableOfOperators().get<Epetra_Operator>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_Operator>(id);
}

/* get const Epetra_Operator from either the const or non-const table
 * using CT_Epetra_Operator_ID */
const Teuchos::RCP<const Epetra_Operator>
CEpetra::getConstOperator( CT_Epetra_Operator_ID_t id )
{
    if (tableOfOperators().isType(id.table))
        return tableOfOperators().getConst<Epetra_Operator>(
        CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_Operator>(
        CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id));
}

/* get const Epetra_Operator from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Operator>
CEpetra::getConstOperator( CTrilinos_Universal_ID_t id )
{
    if (tableOfOperators().isType(id.table))
        return tableOfOperators().getConst<Epetra_Operator>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_Operator>(id);
}

/* store Epetra_Operator (owned) in non-const table */
CT_Epetra_Operator_ID_t
CEpetra::storeNewOperator( Epetra_Operator *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(
        tableOfOperators().store<Epetra_Operator>(pobj, true));
}

/* store Epetra_Operator in non-const table */
CT_Epetra_Operator_ID_t
CEpetra::storeOperator( Epetra_Operator *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(
        tableOfOperators().store<Epetra_Operator>(pobj, false));
}

/* store const Epetra_Operator in const table */
CT_Epetra_Operator_ID_t
CEpetra::storeConstOperator( const Epetra_Operator *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(
        tableOfOperators().store<Epetra_Operator>(pobj, false));
}

/* remove Epetra_Operator from table using CT_Epetra_Operator_ID */
void
CEpetra::removeOperator( CT_Epetra_Operator_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(*id);
    if (tableOfOperators().isType(aid.table))
        tableOfOperators().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(aid);
}

/* remove Epetra_Operator from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeOperator( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfOperators().isType(aid->table))
        tableOfOperators().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_Operator table */
void
CEpetra::purgeOperator(  )
{
    tableOfOperators().purge();
}

/* store Epetra_Operator in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasOperator( const Teuchos::RCP< Epetra_Operator > & robj )
{
    return tableOfOperators().alias(robj);
}

/* store const Epetra_Operator in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstOperator( const Teuchos::RCP< const Epetra_Operator > & robj )
{
    return tableOfOperators().alias(robj);
}



