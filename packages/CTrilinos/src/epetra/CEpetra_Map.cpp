
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
#include "CEpetra_Map.h"
#include "CEpetra_Map_Cpp.hpp"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_Comm_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Map */
Table<Epetra_Map>& tableOfMaps()
{
    static Table<Epetra_Map> loc_tableOfMaps(CT_Epetra_Map_ID);
    return loc_tableOfMaps;
}


} // namespace


//
// Definitions from CEpetra_Map.h
//


extern "C" {


CT_Epetra_Map_ID_t Epetra_Map_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Map_Generalize ( 
  CT_Epetra_Map_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id);
}

CT_Epetra_Map_ID_t Epetra_Map_Create ( 
  int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = CEpetra::getConstComm(CommID);
    return CEpetra::storeNewMap(new Epetra_Map(NumGlobalElements, IndexBase, 
        *Comm));
}

CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = CEpetra::getConstComm(CommID);
    return CEpetra::storeNewMap(new Epetra_Map(NumGlobalElements, 
        NumMyElements, IndexBase, *Comm));
}

CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = CEpetra::getConstComm(CommID);
    return CEpetra::storeNewMap(new Epetra_Map(NumGlobalElements, 
        NumMyElements, MyGlobalElements, IndexBase, *Comm));
}

CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID )
{
    const Teuchos::RCP<const Epetra_Map> map = CEpetra::getConstMap(mapID);
    return CEpetra::storeNewMap(new Epetra_Map(*map));
}

void Epetra_Map_Destroy ( CT_Epetra_Map_ID_t * selfID )
{
    CEpetra::removeMap(selfID);
}

void Epetra_Map_Assign ( 
  CT_Epetra_Map_ID_t selfID, CT_Epetra_Map_ID_t mapID )
{
    Epetra_Map& self = *( CEpetra::getMap(selfID) );

    const Teuchos::RCP<const Epetra_Map> map = CEpetra::getConstMap(mapID);
    self = *map;
}


} // extern "C"


//
// Definitions from CEpetra_Map_Cpp.hpp
//


/* get Epetra_Map from non-const table using CT_Epetra_Map_ID */
const Teuchos::RCP<Epetra_Map>
CEpetra::getMap( CT_Epetra_Map_ID_t id )
{
    if (tableOfMaps().isType(id.table))
        return tableOfMaps().get<Epetra_Map>(
        CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_Map>(
        CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id));
}

/* get Epetra_Map from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Map>
CEpetra::getMap( CTrilinos_Universal_ID_t id )
{
    if (tableOfMaps().isType(id.table))
        return tableOfMaps().get<Epetra_Map>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_Map>(id);
}

/* get const Epetra_Map from either the const or non-const table
 * using CT_Epetra_Map_ID */
const Teuchos::RCP<const Epetra_Map>
CEpetra::getConstMap( CT_Epetra_Map_ID_t id )
{
    if (tableOfMaps().isType(id.table))
        return tableOfMaps().getConst<Epetra_Map>(
        CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_Map>(
        CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id));
}

/* get const Epetra_Map from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Map>
CEpetra::getConstMap( CTrilinos_Universal_ID_t id )
{
    if (tableOfMaps().isType(id.table))
        return tableOfMaps().getConst<Epetra_Map>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_Map>(id);
}

/* store Epetra_Map (owned) in non-const table */
CT_Epetra_Map_ID_t
CEpetra::storeNewMap( Epetra_Map *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        tableOfMaps().store<Epetra_Map>(pobj, true));
}

/* store Epetra_Map in non-const table */
CT_Epetra_Map_ID_t
CEpetra::storeMap( Epetra_Map *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        tableOfMaps().store<Epetra_Map>(pobj, false));
}

/* store const Epetra_Map in const table */
CT_Epetra_Map_ID_t
CEpetra::storeConstMap( const Epetra_Map *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        tableOfMaps().store<Epetra_Map>(pobj, false));
}

/* remove Epetra_Map from table using CT_Epetra_Map_ID */
void
CEpetra::removeMap( CT_Epetra_Map_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Map_ID_t>(*id);
    if (tableOfMaps().isType(aid.table))
        tableOfMaps().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Map_ID_t>(aid);
}

/* remove Epetra_Map from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeMap( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfMaps().isType(aid->table))
        tableOfMaps().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_Map table */
void
CEpetra::purgeMap(  )
{
    tableOfMaps().purge();
}

/* store Epetra_Map in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasMap( const Teuchos::RCP< Epetra_Map > & robj )
{
    return tableOfMaps().alias(robj);
}

/* store const Epetra_Map in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstMap( const Teuchos::RCP< const Epetra_Map > & robj )
{
    return tableOfMaps().alias(robj);
}



