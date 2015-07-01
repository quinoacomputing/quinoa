
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
#include "CEpetra_Distributor.h"
#include "CEpetra_Distributor_Cpp.hpp"
#include "Epetra_Distributor.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Distributor */
Table<Epetra_Distributor>& tableOfDistributors()
{
    static Table<Epetra_Distributor> loc_tableOfDistributors(CT_Epetra_Distributor_ID);
    return loc_tableOfDistributors;
}


} // namespace


//
// Definitions from CEpetra_Distributor.h
//


extern "C" {


CT_Epetra_Distributor_ID_t Epetra_Distributor_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Distributor_Generalize ( 
  CT_Epetra_Distributor_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id);
}

void Epetra_Distributor_Destroy ( 
  CT_Epetra_Distributor_ID_t * selfID )
{
    CEpetra::removeDistributor(selfID);
}

CT_Epetra_Distributor_ID_t Epetra_Distributor_Clone ( 
  CT_Epetra_Distributor_ID_t selfID )
{
    return CEpetra::storeDistributor(CEpetra::getDistributor(selfID)->Clone());
}

int Epetra_Distributor_CreateFromSends ( 
  CT_Epetra_Distributor_ID_t selfID, int NumExportIDs, 
  const int * ExportPIDs, boolean Deterministic, 
  int * NumRemoteIDs )
{
    return CEpetra::getDistributor(selfID)->CreateFromSends(NumExportIDs, 
        ExportPIDs, ((Deterministic) != FALSE ? true : false), *NumRemoteIDs);
}

int Epetra_Distributor_CreateFromRecvs ( 
  CT_Epetra_Distributor_ID_t selfID, int NumRemoteIDs, 
  const int * RemoteGIDs, const int * RemotePIDs, 
  boolean Deterministic, int * NumExportIDs, int ** ExportGIDs, 
  int ** ExportPIDs )
{
    return CEpetra::getDistributor(selfID)->CreateFromRecvs(NumRemoteIDs, 
        RemoteGIDs, RemotePIDs, ((Deterministic) != FALSE ? true : false), 
        *NumExportIDs, *ExportGIDs, *ExportPIDs);
}

int Epetra_Distributor_Do ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->Do(export_objs, obj_size, 
        *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReverse ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReverse(export_objs, obj_size, 
        *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoPosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoPosts(export_objs, obj_size, 
        *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoWaits ( CT_Epetra_Distributor_ID_t selfID )
{
    return CEpetra::getDistributor(selfID)->DoWaits();
}

int Epetra_Distributor_DoReversePosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReversePosts(export_objs, 
        obj_size, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReverseWaits ( 
  CT_Epetra_Distributor_ID_t selfID )
{
    return CEpetra::getDistributor(selfID)->DoReverseWaits();
}

int Epetra_Distributor_Do_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->Do(export_objs, obj_size, *sizes, 
        *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReverse_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReverse(export_objs, obj_size, 
        *sizes, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoPosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoPosts(export_objs, obj_size, 
        *sizes, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReversePosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReversePosts(export_objs, 
        obj_size, *sizes, *len_import_objs, *import_objs);
}


} // extern "C"


//
// Definitions from CEpetra_Distributor_Cpp.hpp
//


/* get Epetra_Distributor from non-const table using CT_Epetra_Distributor_ID */
const Teuchos::RCP<Epetra_Distributor>
CEpetra::getDistributor( CT_Epetra_Distributor_ID_t id )
{
    if (tableOfDistributors().isType(id.table))
        return tableOfDistributors().get<Epetra_Distributor>(
        CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_Distributor>(
        CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id));
}

/* get Epetra_Distributor from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Distributor>
CEpetra::getDistributor( CTrilinos_Universal_ID_t id )
{
    if (tableOfDistributors().isType(id.table))
        return tableOfDistributors().get<Epetra_Distributor>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_Distributor>(id);
}

/* get const Epetra_Distributor from either the const or non-const table
 * using CT_Epetra_Distributor_ID */
const Teuchos::RCP<const Epetra_Distributor>
CEpetra::getConstDistributor( CT_Epetra_Distributor_ID_t id )
{
    if (tableOfDistributors().isType(id.table))
        return tableOfDistributors().getConst<Epetra_Distributor>(
        CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_Distributor>(
        CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id));
}

/* get const Epetra_Distributor from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Distributor>
CEpetra::getConstDistributor( CTrilinos_Universal_ID_t id )
{
    if (tableOfDistributors().isType(id.table))
        return tableOfDistributors().getConst<Epetra_Distributor>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_Distributor>(id);
}

/* store Epetra_Distributor (owned) in non-const table */
CT_Epetra_Distributor_ID_t
CEpetra::storeNewDistributor( Epetra_Distributor *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(
        tableOfDistributors().store<Epetra_Distributor>(pobj, true));
}

/* store Epetra_Distributor in non-const table */
CT_Epetra_Distributor_ID_t
CEpetra::storeDistributor( Epetra_Distributor *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(
        tableOfDistributors().store<Epetra_Distributor>(pobj, false));
}

/* store const Epetra_Distributor in const table */
CT_Epetra_Distributor_ID_t
CEpetra::storeConstDistributor( const Epetra_Distributor *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(
        tableOfDistributors().store<Epetra_Distributor>(pobj, false));
}

/* remove Epetra_Distributor from table using CT_Epetra_Distributor_ID */
void
CEpetra::removeDistributor( CT_Epetra_Distributor_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(*id);
    if (tableOfDistributors().isType(aid.table))
        tableOfDistributors().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(aid);
}

/* remove Epetra_Distributor from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeDistributor( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfDistributors().isType(aid->table))
        tableOfDistributors().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_Distributor table */
void
CEpetra::purgeDistributor(  )
{
    tableOfDistributors().purge();
}

/* store Epetra_Distributor in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasDistributor( const Teuchos::RCP< Epetra_Distributor > & robj )
{
    return tableOfDistributors().alias(robj);
}

/* store const Epetra_Distributor in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstDistributor( const Teuchos::RCP< const Epetra_Distributor > & robj )
{
    return tableOfDistributors().alias(robj);
}



