
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
#include "CEpetra_Export.h"
#include "CEpetra_Export_Cpp.hpp"
#include "Epetra_Export.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Distributor_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Export */
Table<Epetra_Export>& tableOfExports()
{
    static Table<Epetra_Export> loc_tableOfExports(CT_Epetra_Export_ID);
    return loc_tableOfExports;
}


} // namespace


//
// Definitions from CEpetra_Export.h
//


extern "C" {


CT_Epetra_Export_ID_t Epetra_Export_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Export_Generalize ( 
  CT_Epetra_Export_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id);
}

CT_Epetra_Export_ID_t Epetra_Export_Create ( 
  CT_Epetra_BlockMap_ID_t SourceMapID, 
  CT_Epetra_BlockMap_ID_t TargetMapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> SourceMap = 
        CEpetra::getConstBlockMap(SourceMapID);
    const Teuchos::RCP<const Epetra_BlockMap> TargetMap = 
        CEpetra::getConstBlockMap(TargetMapID);
    return CEpetra::storeNewExport(new Epetra_Export(*SourceMap, *TargetMap));
}

CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( 
  CT_Epetra_Export_ID_t ExporterID )
{
    const Teuchos::RCP<const Epetra_Export> Exporter = CEpetra::getConstExport(
        ExporterID);
    return CEpetra::storeNewExport(new Epetra_Export(*Exporter));
}

void Epetra_Export_Destroy ( CT_Epetra_Export_ID_t * selfID )
{
    CEpetra::removeExport(selfID);
}

int Epetra_Export_NumSameIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->NumSameIDs();
}

int Epetra_Export_NumPermuteIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->NumPermuteIDs();
}

int * Epetra_Export_PermuteFromLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->PermuteFromLIDs();
}

int * Epetra_Export_PermuteToLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->PermuteToLIDs();
}

int Epetra_Export_NumRemoteIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->NumRemoteIDs();
}

int * Epetra_Export_RemoteLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->RemoteLIDs();
}

int Epetra_Export_NumExportIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->NumExportIDs();
}

int * Epetra_Export_ExportLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->ExportLIDs();
}

int * Epetra_Export_ExportPIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->ExportPIDs();
}

int Epetra_Export_NumSend ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->NumSend();
}

int Epetra_Export_NumRecv ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->NumRecv();
}

CT_Epetra_BlockMap_ID_t Epetra_Export_SourceMap ( 
  CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstExport(
        selfID)->SourceMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_Export_TargetMap ( 
  CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstExport(
        selfID)->TargetMap() ));
}

CT_Epetra_Distributor_ID_t Epetra_Export_Distributor ( 
  CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::storeDistributor(&( CEpetra::getConstExport(
        selfID)->Distributor() ));
}


} // extern "C"


//
// Definitions from CEpetra_Export_Cpp.hpp
//


/* get Epetra_Export from non-const table using CT_Epetra_Export_ID */
const Teuchos::RCP<Epetra_Export>
CEpetra::getExport( CT_Epetra_Export_ID_t id )
{
    if (tableOfExports().isType(id.table))
        return tableOfExports().get<Epetra_Export>(
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_Export>(
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id));
}

/* get Epetra_Export from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Export>
CEpetra::getExport( CTrilinos_Universal_ID_t id )
{
    if (tableOfExports().isType(id.table))
        return tableOfExports().get<Epetra_Export>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_Export>(id);
}

/* get const Epetra_Export from either the const or non-const table
 * using CT_Epetra_Export_ID */
const Teuchos::RCP<const Epetra_Export>
CEpetra::getConstExport( CT_Epetra_Export_ID_t id )
{
    if (tableOfExports().isType(id.table))
        return tableOfExports().getConst<Epetra_Export>(
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_Export>(
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id));
}

/* get const Epetra_Export from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Export>
CEpetra::getConstExport( CTrilinos_Universal_ID_t id )
{
    if (tableOfExports().isType(id.table))
        return tableOfExports().getConst<Epetra_Export>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_Export>(id);
}

/* store Epetra_Export (owned) in non-const table */
CT_Epetra_Export_ID_t
CEpetra::storeNewExport( Epetra_Export *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        tableOfExports().store<Epetra_Export>(pobj, true));
}

/* store Epetra_Export in non-const table */
CT_Epetra_Export_ID_t
CEpetra::storeExport( Epetra_Export *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        tableOfExports().store<Epetra_Export>(pobj, false));
}

/* store const Epetra_Export in const table */
CT_Epetra_Export_ID_t
CEpetra::storeConstExport( const Epetra_Export *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        tableOfExports().store<Epetra_Export>(pobj, false));
}

/* remove Epetra_Export from table using CT_Epetra_Export_ID */
void
CEpetra::removeExport( CT_Epetra_Export_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(*id);
    if (tableOfExports().isType(aid.table))
        tableOfExports().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Export_ID_t>(aid);
}

/* remove Epetra_Export from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeExport( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfExports().isType(aid->table))
        tableOfExports().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_Export table */
void
CEpetra::purgeExport(  )
{
    tableOfExports().purge();
}

/* store Epetra_Export in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasExport( const Teuchos::RCP< Epetra_Export > & robj )
{
    return tableOfExports().alias(robj);
}

/* store const Epetra_Export in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstExport( const Teuchos::RCP< const Epetra_Export > & robj )
{
    return tableOfExports().alias(robj);
}



