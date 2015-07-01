
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
#include "CEpetra_Directory.h"
#include "CEpetra_Directory_Cpp.hpp"
#include "Epetra_Directory.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Directory */
Table<Epetra_Directory>& tableOfDirectorys()
{
    static Table<Epetra_Directory> loc_tableOfDirectorys(CT_Epetra_Directory_ID);
    return loc_tableOfDirectorys;
}


} // namespace


//
// Definitions from CEpetra_Directory.h
//


extern "C" {


CT_Epetra_Directory_ID_t Epetra_Directory_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Directory_Generalize ( 
  CT_Epetra_Directory_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(id);
}

void Epetra_Directory_Destroy ( CT_Epetra_Directory_ID_t * selfID )
{
    CEpetra::removeDirectory(selfID);
}

int Epetra_Directory_GetDirectoryEntries ( 
  CT_Epetra_Directory_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID, 
  const int NumEntries, const int * GlobalEntries, int * Procs, 
  int * LocalEntries, int * EntrySizes, 
  boolean high_rank_sharing_procs )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return CEpetra::getConstDirectory(selfID)->GetDirectoryEntries(*Map, 
        NumEntries, GlobalEntries, Procs, LocalEntries, EntrySizes, ((
        high_rank_sharing_procs) != FALSE ? true : false));
}

boolean Epetra_Directory_GIDsAllUniquelyOwned ( 
  CT_Epetra_Directory_ID_t selfID )
{
    return ((CEpetra::getConstDirectory(
        selfID)->GIDsAllUniquelyOwned()) ? TRUE : FALSE);
}


} // extern "C"


//
// Definitions from CEpetra_Directory_Cpp.hpp
//


/* get Epetra_Directory from non-const table using CT_Epetra_Directory_ID */
const Teuchos::RCP<Epetra_Directory>
CEpetra::getDirectory( CT_Epetra_Directory_ID_t id )
{
    if (tableOfDirectorys().isType(id.table))
        return tableOfDirectorys().get<Epetra_Directory>(
        CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_Directory>(
        CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(id));
}

/* get Epetra_Directory from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Directory>
CEpetra::getDirectory( CTrilinos_Universal_ID_t id )
{
    if (tableOfDirectorys().isType(id.table))
        return tableOfDirectorys().get<Epetra_Directory>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_Directory>(id);
}

/* get const Epetra_Directory from either the const or non-const table
 * using CT_Epetra_Directory_ID */
const Teuchos::RCP<const Epetra_Directory>
CEpetra::getConstDirectory( CT_Epetra_Directory_ID_t id )
{
    if (tableOfDirectorys().isType(id.table))
        return tableOfDirectorys().getConst<Epetra_Directory>(
        CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_Directory>(
        CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(id));
}

/* get const Epetra_Directory from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Directory>
CEpetra::getConstDirectory( CTrilinos_Universal_ID_t id )
{
    if (tableOfDirectorys().isType(id.table))
        return tableOfDirectorys().getConst<Epetra_Directory>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_Directory>(id);
}

/* store Epetra_Directory (owned) in non-const table */
CT_Epetra_Directory_ID_t
CEpetra::storeNewDirectory( Epetra_Directory *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(
        tableOfDirectorys().store<Epetra_Directory>(pobj, true));
}

/* store Epetra_Directory in non-const table */
CT_Epetra_Directory_ID_t
CEpetra::storeDirectory( Epetra_Directory *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(
        tableOfDirectorys().store<Epetra_Directory>(pobj, false));
}

/* store const Epetra_Directory in const table */
CT_Epetra_Directory_ID_t
CEpetra::storeConstDirectory( const Epetra_Directory *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(
        tableOfDirectorys().store<Epetra_Directory>(pobj, false));
}

/* remove Epetra_Directory from table using CT_Epetra_Directory_ID */
void
CEpetra::removeDirectory( CT_Epetra_Directory_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(*id);
    if (tableOfDirectorys().isType(aid.table))
        tableOfDirectorys().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(aid);
}

/* remove Epetra_Directory from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeDirectory( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfDirectorys().isType(aid->table))
        tableOfDirectorys().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_Directory table */
void
CEpetra::purgeDirectory(  )
{
    tableOfDirectorys().purge();
}

/* store Epetra_Directory in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasDirectory( const Teuchos::RCP< Epetra_Directory > & robj )
{
    return tableOfDirectorys().alias(robj);
}

/* store const Epetra_Directory in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstDirectory( const Teuchos::RCP< const Epetra_Directory > & robj )
{
    return tableOfDirectorys().alias(robj);
}



