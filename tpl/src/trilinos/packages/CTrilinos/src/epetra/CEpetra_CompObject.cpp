
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
#include "CEpetra_CompObject.h"
#include "CEpetra_CompObject_Cpp.hpp"
#include "Epetra_CompObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_Flops_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_CompObject */
Table<Epetra_CompObject>& tableOfCompObjects()
{
    static Table<Epetra_CompObject> loc_tableOfCompObjects(CT_Epetra_CompObject_ID);
    return loc_tableOfCompObjects;
}


} // namespace


//
// Definitions from CEpetra_CompObject.h
//


extern "C" {


CT_Epetra_CompObject_ID_t Epetra_CompObject_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_CompObject_Generalize ( 
  CT_Epetra_CompObject_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(id);
}

CT_Epetra_CompObject_ID_t Epetra_CompObject_Create (  )
{
    return CEpetra::storeNewCompObject(new Epetra_CompObject());
}

CT_Epetra_CompObject_ID_t Epetra_CompObject_Duplicate ( 
  CT_Epetra_CompObject_ID_t SourceID )
{
    const Teuchos::RCP<const Epetra_CompObject> Source = 
        CEpetra::getConstCompObject(SourceID);
    return CEpetra::storeNewCompObject(new Epetra_CompObject(*Source));
}

void Epetra_CompObject_Destroy ( CT_Epetra_CompObject_ID_t * selfID )
{
    CEpetra::removeCompObject(selfID);
}

void Epetra_CompObject_SetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_Flops_ID_t FlopCounter_inID )
{
    const Teuchos::RCP<const Epetra_Flops> FlopCounter_in = 
        CEpetra::getConstFlops(FlopCounter_inID);
    CEpetra::getCompObject(selfID)->SetFlopCounter(*FlopCounter_in);
}

void Epetra_CompObject_SetFlopCounter_Matching ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_CompObject_ID_t CompObjectID )
{
    const Teuchos::RCP<const Epetra_CompObject> CompObject = 
        CEpetra::getConstCompObject(CompObjectID);
    CEpetra::getCompObject(selfID)->SetFlopCounter(*CompObject);
}

void Epetra_CompObject_UnsetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID )
{
    CEpetra::getCompObject(selfID)->UnsetFlopCounter();
}

CT_Epetra_Flops_ID_t Epetra_CompObject_GetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID )
{
    return CEpetra::storeFlops(CEpetra::getConstCompObject(
        selfID)->GetFlopCounter());
}

void Epetra_CompObject_ResetFlops ( 
  CT_Epetra_CompObject_ID_t selfID )
{
    CEpetra::getConstCompObject(selfID)->ResetFlops();
}

double Epetra_CompObject_Flops ( CT_Epetra_CompObject_ID_t selfID )
{
    return CEpetra::getConstCompObject(selfID)->Flops();
}

void Epetra_CompObject_UpdateFlops_Int ( 
  CT_Epetra_CompObject_ID_t selfID, int Flops_in )
{
    CEpetra::getConstCompObject(selfID)->UpdateFlops(Flops_in);
}

void Epetra_CompObject_UpdateFlops_Long ( 
  CT_Epetra_CompObject_ID_t selfID, long int Flops_in )
{
    CEpetra::getConstCompObject(selfID)->UpdateFlops(Flops_in);
}

void Epetra_CompObject_UpdateFlops_Double ( 
  CT_Epetra_CompObject_ID_t selfID, double Flops_in )
{
    CEpetra::getConstCompObject(selfID)->UpdateFlops(Flops_in);
}

void Epetra_CompObject_UpdateFlops_Float ( 
  CT_Epetra_CompObject_ID_t selfID, float Flops_in )
{
    CEpetra::getConstCompObject(selfID)->UpdateFlops(Flops_in);
}

void Epetra_CompObject_Assign ( 
  CT_Epetra_CompObject_ID_t selfID, CT_Epetra_CompObject_ID_t srcID )
{
    Epetra_CompObject& self = *( CEpetra::getCompObject(selfID) );

    const Teuchos::RCP<const Epetra_CompObject> src = 
        CEpetra::getConstCompObject(srcID);
    self = *src;
}


} // extern "C"


//
// Definitions from CEpetra_CompObject_Cpp.hpp
//


/* get Epetra_CompObject from non-const table using CT_Epetra_CompObject_ID */
const Teuchos::RCP<Epetra_CompObject>
CEpetra::getCompObject( CT_Epetra_CompObject_ID_t id )
{
    if (tableOfCompObjects().isType(id.table))
        return tableOfCompObjects().get<Epetra_CompObject>(
        CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_CompObject>(
        CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(id));
}

/* get Epetra_CompObject from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_CompObject>
CEpetra::getCompObject( CTrilinos_Universal_ID_t id )
{
    if (tableOfCompObjects().isType(id.table))
        return tableOfCompObjects().get<Epetra_CompObject>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_CompObject>(id);
}

/* get const Epetra_CompObject from either the const or non-const table
 * using CT_Epetra_CompObject_ID */
const Teuchos::RCP<const Epetra_CompObject>
CEpetra::getConstCompObject( CT_Epetra_CompObject_ID_t id )
{
    if (tableOfCompObjects().isType(id.table))
        return tableOfCompObjects().getConst<Epetra_CompObject>(
        CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_CompObject>(
        CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(id));
}

/* get const Epetra_CompObject from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_CompObject>
CEpetra::getConstCompObject( CTrilinos_Universal_ID_t id )
{
    if (tableOfCompObjects().isType(id.table))
        return tableOfCompObjects().getConst<Epetra_CompObject>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_CompObject>(id);
}

/* store Epetra_CompObject (owned) in non-const table */
CT_Epetra_CompObject_ID_t
CEpetra::storeNewCompObject( Epetra_CompObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(
        tableOfCompObjects().store<Epetra_CompObject>(pobj, true));
}

/* store Epetra_CompObject in non-const table */
CT_Epetra_CompObject_ID_t
CEpetra::storeCompObject( Epetra_CompObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(
        tableOfCompObjects().store<Epetra_CompObject>(pobj, false));
}

/* store const Epetra_CompObject in const table */
CT_Epetra_CompObject_ID_t
CEpetra::storeConstCompObject( const Epetra_CompObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(
        tableOfCompObjects().store<Epetra_CompObject>(pobj, false));
}

/* remove Epetra_CompObject from table using CT_Epetra_CompObject_ID */
void
CEpetra::removeCompObject( CT_Epetra_CompObject_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(*id);
    if (tableOfCompObjects().isType(aid.table))
        tableOfCompObjects().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(aid);
}

/* remove Epetra_CompObject from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeCompObject( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfCompObjects().isType(aid->table))
        tableOfCompObjects().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_CompObject table */
void
CEpetra::purgeCompObject(  )
{
    tableOfCompObjects().purge();
}

/* store Epetra_CompObject in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasCompObject( const Teuchos::RCP< Epetra_CompObject > & robj )
{
    return tableOfCompObjects().alias(robj);
}

/* store const Epetra_CompObject in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstCompObject( const Teuchos::RCP< const Epetra_CompObject > & robj )
{
    return tableOfCompObjects().alias(robj);
}



