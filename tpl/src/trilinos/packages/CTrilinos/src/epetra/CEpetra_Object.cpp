
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
#include "CEpetra_Object.h"
#include "CEpetra_Object_Cpp.hpp"
#include "Epetra_Object.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Object */
Table<Epetra_Object>& tableOfObjects()
{
    static Table<Epetra_Object> loc_tableOfObjects(CT_Epetra_Object_ID);
    return loc_tableOfObjects;
}


} // namespace


//
// Definitions from CEpetra_Object.h
//


extern "C" {


CT_Epetra_Object_ID_t Epetra_Object_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Object_Generalize ( 
  CT_Epetra_Object_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Object_ID_t>(id);
}

CT_Epetra_Object_ID_t Epetra_Object_Create ( 
  int TracebackModeIn, boolean set_label )
{
    return CEpetra::storeNewObject(new Epetra_Object(TracebackModeIn, ((
        set_label) != FALSE ? true : false)));
}

CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( 
  const char * const Label, int TracebackModeIn )
{
    return CEpetra::storeNewObject(new Epetra_Object(Label, TracebackModeIn));
}

CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( 
  CT_Epetra_Object_ID_t ObjectID )
{
    const Teuchos::RCP<const Epetra_Object> Object = CEpetra::getConstObject(
        ObjectID);
    return CEpetra::storeNewObject(new Epetra_Object(*Object));
}

void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID )
{
    CEpetra::removeObject(selfID);
}

void Epetra_Object_SetLabel ( 
  CT_Epetra_Object_ID_t selfID, const char * const Label )
{
    CEpetra::getObject(selfID)->SetLabel(Label);
}

const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID )
{
    return CEpetra::getConstObject(selfID)->Label();
}

int Epetra_Object_ReportError ( 
  CT_Epetra_Object_ID_t selfID, const char Message[], int ErrorCode )
{
    return CEpetra::getConstObject(selfID)->ReportError(std::string(Message), 
        ErrorCode);
}

void Epetra_Object_SetTracebackMode ( int TracebackModeValue )
{
    Epetra_Object::SetTracebackMode(TracebackModeValue);
}

int Epetra_Object_GetTracebackMode (  )
{
    return Epetra_Object::GetTracebackMode();
}


} // extern "C"


//
// Definitions from CEpetra_Object_Cpp.hpp
//


/* get Epetra_Object from non-const table using CT_Epetra_Object_ID */
const Teuchos::RCP<Epetra_Object>
CEpetra::getObject( CT_Epetra_Object_ID_t id )
{
    if (tableOfObjects().isType(id.table))
        return tableOfObjects().get<Epetra_Object>(
        CTrilinos::abstractType<CT_Epetra_Object_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_Object>(
        CTrilinos::abstractType<CT_Epetra_Object_ID_t>(id));
}

/* get Epetra_Object from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Object>
CEpetra::getObject( CTrilinos_Universal_ID_t id )
{
    if (tableOfObjects().isType(id.table))
        return tableOfObjects().get<Epetra_Object>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_Object>(id);
}

/* get const Epetra_Object from either the const or non-const table
 * using CT_Epetra_Object_ID */
const Teuchos::RCP<const Epetra_Object>
CEpetra::getConstObject( CT_Epetra_Object_ID_t id )
{
    if (tableOfObjects().isType(id.table))
        return tableOfObjects().getConst<Epetra_Object>(
        CTrilinos::abstractType<CT_Epetra_Object_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_Object>(
        CTrilinos::abstractType<CT_Epetra_Object_ID_t>(id));
}

/* get const Epetra_Object from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Object>
CEpetra::getConstObject( CTrilinos_Universal_ID_t id )
{
    if (tableOfObjects().isType(id.table))
        return tableOfObjects().getConst<Epetra_Object>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_Object>(id);
}

/* store Epetra_Object (owned) in non-const table */
CT_Epetra_Object_ID_t
CEpetra::storeNewObject( Epetra_Object *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
        tableOfObjects().store<Epetra_Object>(pobj, true));
}

/* store Epetra_Object in non-const table */
CT_Epetra_Object_ID_t
CEpetra::storeObject( Epetra_Object *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
        tableOfObjects().store<Epetra_Object>(pobj, false));
}

/* store const Epetra_Object in const table */
CT_Epetra_Object_ID_t
CEpetra::storeConstObject( const Epetra_Object *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
        tableOfObjects().store<Epetra_Object>(pobj, false));
}

/* remove Epetra_Object from table using CT_Epetra_Object_ID */
void
CEpetra::removeObject( CT_Epetra_Object_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Object_ID_t>(*id);
    if (tableOfObjects().isType(aid.table))
        tableOfObjects().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Object_ID_t>(aid);
}

/* remove Epetra_Object from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeObject( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfObjects().isType(aid->table))
        tableOfObjects().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_Object table */
void
CEpetra::purgeObject(  )
{
    tableOfObjects().purge();
}

/* store Epetra_Object in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasObject( const Teuchos::RCP< Epetra_Object > & robj )
{
    return tableOfObjects().alias(robj);
}

/* store const Epetra_Object in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstObject( const Teuchos::RCP< const Epetra_Object > & robj )
{
    return tableOfObjects().alias(robj);
}



