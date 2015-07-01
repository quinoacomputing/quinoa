
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
#include "CEpetra_SrcDistObject.h"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "Epetra_SrcDistObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_SrcDistObject */
Table<Epetra_SrcDistObject>& tableOfSrcDistObjects()
{
    static Table<Epetra_SrcDistObject> loc_tableOfSrcDistObjects(CT_Epetra_SrcDistObject_ID);
    return loc_tableOfSrcDistObjects;
}


} // namespace


//
// Definitions from CEpetra_SrcDistObject.h
//


extern "C" {


CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_SrcDistObject_Generalize ( 
  CT_Epetra_SrcDistObject_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(id);
}

void Epetra_SrcDistObject_Destroy ( 
  CT_Epetra_SrcDistObject_ID_t * selfID )
{
    CEpetra::removeSrcDistObject(selfID);
}

CT_Epetra_BlockMap_ID_t Epetra_SrcDistObject_Map ( 
  CT_Epetra_SrcDistObject_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstSrcDistObject(
        selfID)->Map() ));
}


} // extern "C"


//
// Definitions from CEpetra_SrcDistObject_Cpp.hpp
//


/* get Epetra_SrcDistObject from non-const table using CT_Epetra_SrcDistObject_ID */
const Teuchos::RCP<Epetra_SrcDistObject>
CEpetra::getSrcDistObject( CT_Epetra_SrcDistObject_ID_t id )
{
    if (tableOfSrcDistObjects().isType(id.table))
        return tableOfSrcDistObjects().get<Epetra_SrcDistObject>(
        CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_SrcDistObject>(
        CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(id));
}

/* get Epetra_SrcDistObject from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_SrcDistObject>
CEpetra::getSrcDistObject( CTrilinos_Universal_ID_t id )
{
    if (tableOfSrcDistObjects().isType(id.table))
        return tableOfSrcDistObjects().get<Epetra_SrcDistObject>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_SrcDistObject>(id);
}

/* get const Epetra_SrcDistObject from either the const or non-const table
 * using CT_Epetra_SrcDistObject_ID */
const Teuchos::RCP<const Epetra_SrcDistObject>
CEpetra::getConstSrcDistObject( CT_Epetra_SrcDistObject_ID_t id )
{
    if (tableOfSrcDistObjects().isType(id.table))
        return tableOfSrcDistObjects().getConst<Epetra_SrcDistObject>(
        CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_SrcDistObject>(
        CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(id));
}

/* get const Epetra_SrcDistObject from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_SrcDistObject>
CEpetra::getConstSrcDistObject( CTrilinos_Universal_ID_t id )
{
    if (tableOfSrcDistObjects().isType(id.table))
        return tableOfSrcDistObjects().getConst<Epetra_SrcDistObject>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_SrcDistObject>(id);
}

/* store Epetra_SrcDistObject (owned) in non-const table */
CT_Epetra_SrcDistObject_ID_t
CEpetra::storeNewSrcDistObject( Epetra_SrcDistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(
        tableOfSrcDistObjects().store<Epetra_SrcDistObject>(pobj, true));
}

/* store Epetra_SrcDistObject in non-const table */
CT_Epetra_SrcDistObject_ID_t
CEpetra::storeSrcDistObject( Epetra_SrcDistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(
        tableOfSrcDistObjects().store<Epetra_SrcDistObject>(pobj, false));
}

/* store const Epetra_SrcDistObject in const table */
CT_Epetra_SrcDistObject_ID_t
CEpetra::storeConstSrcDistObject( const Epetra_SrcDistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(
        tableOfSrcDistObjects().store<Epetra_SrcDistObject>(pobj, false));
}

/* remove Epetra_SrcDistObject from table using CT_Epetra_SrcDistObject_ID */
void
CEpetra::removeSrcDistObject( CT_Epetra_SrcDistObject_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(*id);
    if (tableOfSrcDistObjects().isType(aid.table))
        tableOfSrcDistObjects().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(aid);
}

/* remove Epetra_SrcDistObject from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeSrcDistObject( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfSrcDistObjects().isType(aid->table))
        tableOfSrcDistObjects().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_SrcDistObject table */
void
CEpetra::purgeSrcDistObject(  )
{
    tableOfSrcDistObjects().purge();
}

/* store Epetra_SrcDistObject in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasSrcDistObject( const Teuchos::RCP< Epetra_SrcDistObject > & robj )
{
    return tableOfSrcDistObjects().alias(robj);
}

/* store const Epetra_SrcDistObject in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstSrcDistObject( const Teuchos::RCP< const Epetra_SrcDistObject > & robj )
{
    return tableOfSrcDistObjects().alias(robj);
}



