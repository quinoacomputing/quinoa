
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
#include "CEpetra_OffsetIndex.h"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "Epetra_OffsetIndex.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_OffsetIndex */
Table<Epetra_OffsetIndex>& tableOfOffsetIndexs()
{
    static Table<Epetra_OffsetIndex> loc_tableOfOffsetIndexs(CT_Epetra_OffsetIndex_ID);
    return loc_tableOfOffsetIndexs;
}


} // namespace


//
// Definitions from CEpetra_OffsetIndex.h
//


extern "C" {


CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_OffsetIndex_Generalize ( 
  CT_Epetra_OffsetIndex_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id);
}

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromImporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Import_ID_t ImporterID )
{
    const Teuchos::RCP<const Epetra_CrsGraph> SourceGraph = 
        CEpetra::getConstCrsGraph(SourceGraphID);
    const Teuchos::RCP<const Epetra_CrsGraph> TargetGraph = 
        CEpetra::getConstCrsGraph(TargetGraphID);
    const Teuchos::RCP<Epetra_Import> Importer = CEpetra::getImport(
        ImporterID);
    return CEpetra::storeNewOffsetIndex(new Epetra_OffsetIndex(*SourceGraph, 
        *TargetGraph, *Importer));
}

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromExporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Export_ID_t ExporterID )
{
    const Teuchos::RCP<const Epetra_CrsGraph> SourceGraph = 
        CEpetra::getConstCrsGraph(SourceGraphID);
    const Teuchos::RCP<const Epetra_CrsGraph> TargetGraph = 
        CEpetra::getConstCrsGraph(TargetGraphID);
    const Teuchos::RCP<Epetra_Export> Exporter = CEpetra::getExport(
        ExporterID);
    return CEpetra::storeNewOffsetIndex(new Epetra_OffsetIndex(*SourceGraph, 
        *TargetGraph, *Exporter));
}

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Duplicate ( 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<const Epetra_OffsetIndex> Indexor = 
        CEpetra::getConstOffsetIndex(IndexorID);
    return CEpetra::storeNewOffsetIndex(new Epetra_OffsetIndex(*Indexor));
}

void Epetra_OffsetIndex_Destroy ( 
  CT_Epetra_OffsetIndex_ID_t * selfID )
{
    CEpetra::removeOffsetIndex(selfID);
}

int ** Epetra_OffsetIndex_SameOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID )
{
    return CEpetra::getConstOffsetIndex(selfID)->SameOffsets();
}

int ** Epetra_OffsetIndex_PermuteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID )
{
    return CEpetra::getConstOffsetIndex(selfID)->PermuteOffsets();
}

int ** Epetra_OffsetIndex_RemoteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID )
{
    return CEpetra::getConstOffsetIndex(selfID)->RemoteOffsets();
}


} // extern "C"


//
// Definitions from CEpetra_OffsetIndex_Cpp.hpp
//


/* get Epetra_OffsetIndex from non-const table using CT_Epetra_OffsetIndex_ID */
const Teuchos::RCP<Epetra_OffsetIndex>
CEpetra::getOffsetIndex( CT_Epetra_OffsetIndex_ID_t id )
{
    if (tableOfOffsetIndexs().isType(id.table))
        return tableOfOffsetIndexs().get<Epetra_OffsetIndex>(
        CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_OffsetIndex>(
        CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id));
}

/* get Epetra_OffsetIndex from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_OffsetIndex>
CEpetra::getOffsetIndex( CTrilinos_Universal_ID_t id )
{
    if (tableOfOffsetIndexs().isType(id.table))
        return tableOfOffsetIndexs().get<Epetra_OffsetIndex>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_OffsetIndex>(id);
}

/* get const Epetra_OffsetIndex from either the const or non-const table
 * using CT_Epetra_OffsetIndex_ID */
const Teuchos::RCP<const Epetra_OffsetIndex>
CEpetra::getConstOffsetIndex( CT_Epetra_OffsetIndex_ID_t id )
{
    if (tableOfOffsetIndexs().isType(id.table))
        return tableOfOffsetIndexs().getConst<Epetra_OffsetIndex>(
        CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_OffsetIndex>(
        CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id));
}

/* get const Epetra_OffsetIndex from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_OffsetIndex>
CEpetra::getConstOffsetIndex( CTrilinos_Universal_ID_t id )
{
    if (tableOfOffsetIndexs().isType(id.table))
        return tableOfOffsetIndexs().getConst<Epetra_OffsetIndex>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_OffsetIndex>(id);
}

/* store Epetra_OffsetIndex (owned) in non-const table */
CT_Epetra_OffsetIndex_ID_t
CEpetra::storeNewOffsetIndex( Epetra_OffsetIndex *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
        tableOfOffsetIndexs().store<Epetra_OffsetIndex>(pobj, true));
}

/* store Epetra_OffsetIndex in non-const table */
CT_Epetra_OffsetIndex_ID_t
CEpetra::storeOffsetIndex( Epetra_OffsetIndex *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
        tableOfOffsetIndexs().store<Epetra_OffsetIndex>(pobj, false));
}

/* store const Epetra_OffsetIndex in const table */
CT_Epetra_OffsetIndex_ID_t
CEpetra::storeConstOffsetIndex( const Epetra_OffsetIndex *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
        tableOfOffsetIndexs().store<Epetra_OffsetIndex>(pobj, false));
}

/* remove Epetra_OffsetIndex from table using CT_Epetra_OffsetIndex_ID */
void
CEpetra::removeOffsetIndex( CT_Epetra_OffsetIndex_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(*id);
    if (tableOfOffsetIndexs().isType(aid.table))
        tableOfOffsetIndexs().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(aid);
}

/* remove Epetra_OffsetIndex from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeOffsetIndex( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfOffsetIndexs().isType(aid->table))
        tableOfOffsetIndexs().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_OffsetIndex table */
void
CEpetra::purgeOffsetIndex(  )
{
    tableOfOffsetIndexs().purge();
}

/* store Epetra_OffsetIndex in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasOffsetIndex( const Teuchos::RCP< Epetra_OffsetIndex > & robj )
{
    return tableOfOffsetIndexs().alias(robj);
}

/* store const Epetra_OffsetIndex in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstOffsetIndex( const Teuchos::RCP< const Epetra_OffsetIndex > & robj )
{
    return tableOfOffsetIndexs().alias(robj);
}



