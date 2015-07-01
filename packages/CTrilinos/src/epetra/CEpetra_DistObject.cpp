
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
#include "CEpetra_DistObject.h"
#include "CEpetra_DistObject_Cpp.hpp"
#include "Epetra_DistObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_DistObject */
Table<Epetra_DistObject>& tableOfDistObjects()
{
    static Table<Epetra_DistObject> loc_tableOfDistObjects(CT_Epetra_DistObject_ID);
    return loc_tableOfDistObjects;
}


} // namespace


//
// Definitions from CEpetra_DistObject.h
//


extern "C" {


CT_Epetra_DistObject_ID_t Epetra_DistObject_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_DistObject_Generalize ( 
  CT_Epetra_DistObject_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id);
}

void Epetra_DistObject_Destroy ( CT_Epetra_DistObject_ID_t * selfID )
{
    CEpetra::removeDistObject(selfID);
}

int Epetra_DistObject_Import ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<const Epetra_SrcDistObject> A = 
        CEpetra::getConstSrcDistObject(AID);
    const Teuchos::RCP<const Epetra_Import> Importer = CEpetra::getConstImport(
        ImporterID);
    const Teuchos::RCP<const Epetra_OffsetIndex> Indexor = 
        CEpetra::getConstOffsetIndex(IndexorID);
    return CEpetra::getDistObject(selfID)->Import(*A, *Importer, 
        (Epetra_CombineMode) CombineMode, Indexor.getRawPtr());
}

int Epetra_DistObject_Import_UsingExporter ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<const Epetra_SrcDistObject> A = 
        CEpetra::getConstSrcDistObject(AID);
    const Teuchos::RCP<const Epetra_Export> Exporter = CEpetra::getConstExport(
        ExporterID);
    const Teuchos::RCP<const Epetra_OffsetIndex> Indexor = 
        CEpetra::getConstOffsetIndex(IndexorID);
    return CEpetra::getDistObject(selfID)->Import(*A, *Exporter, 
        (Epetra_CombineMode) CombineMode, Indexor.getRawPtr());
}

int Epetra_DistObject_Export_UsingImporter ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<const Epetra_SrcDistObject> A = 
        CEpetra::getConstSrcDistObject(AID);
    const Teuchos::RCP<const Epetra_Import> Importer = CEpetra::getConstImport(
        ImporterID);
    const Teuchos::RCP<const Epetra_OffsetIndex> Indexor = 
        CEpetra::getConstOffsetIndex(IndexorID);
    return CEpetra::getDistObject(selfID)->Export(*A, *Importer, 
        (Epetra_CombineMode) CombineMode, Indexor.getRawPtr());
}

int Epetra_DistObject_Export ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<const Epetra_SrcDistObject> A = 
        CEpetra::getConstSrcDistObject(AID);
    const Teuchos::RCP<const Epetra_Export> Exporter = CEpetra::getConstExport(
        ExporterID);
    const Teuchos::RCP<const Epetra_OffsetIndex> Indexor = 
        CEpetra::getConstOffsetIndex(IndexorID);
    return CEpetra::getDistObject(selfID)->Export(*A, *Exporter, 
        (Epetra_CombineMode) CombineMode, Indexor.getRawPtr());
}

CT_Epetra_BlockMap_ID_t Epetra_DistObject_Map ( 
  CT_Epetra_DistObject_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstDistObject(
        selfID)->Map() ));
}

CT_Epetra_Comm_ID_t Epetra_DistObject_Comm ( 
  CT_Epetra_DistObject_ID_t selfID )
{
    return CEpetra::storeConstComm(&( CEpetra::getConstDistObject(
        selfID)->Comm() ));
}

boolean Epetra_DistObject_DistributedGlobal ( 
  CT_Epetra_DistObject_ID_t selfID )
{
    return ((CEpetra::getConstDistObject(
        selfID)->DistributedGlobal()) ? TRUE : FALSE);
}


} // extern "C"


//
// Definitions from CEpetra_DistObject_Cpp.hpp
//


/* get Epetra_DistObject from non-const table using CT_Epetra_DistObject_ID */
const Teuchos::RCP<Epetra_DistObject>
CEpetra::getDistObject( CT_Epetra_DistObject_ID_t id )
{
    if (tableOfDistObjects().isType(id.table))
        return tableOfDistObjects().get<Epetra_DistObject>(
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_DistObject>(
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id));
}

/* get Epetra_DistObject from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_DistObject>
CEpetra::getDistObject( CTrilinos_Universal_ID_t id )
{
    if (tableOfDistObjects().isType(id.table))
        return tableOfDistObjects().get<Epetra_DistObject>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_DistObject>(id);
}

/* get const Epetra_DistObject from either the const or non-const table
 * using CT_Epetra_DistObject_ID */
const Teuchos::RCP<const Epetra_DistObject>
CEpetra::getConstDistObject( CT_Epetra_DistObject_ID_t id )
{
    if (tableOfDistObjects().isType(id.table))
        return tableOfDistObjects().getConst<Epetra_DistObject>(
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_DistObject>(
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id));
}

/* get const Epetra_DistObject from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_DistObject>
CEpetra::getConstDistObject( CTrilinos_Universal_ID_t id )
{
    if (tableOfDistObjects().isType(id.table))
        return tableOfDistObjects().getConst<Epetra_DistObject>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_DistObject>(id);
}

/* store Epetra_DistObject (owned) in non-const table */
CT_Epetra_DistObject_ID_t
CEpetra::storeNewDistObject( Epetra_DistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
        tableOfDistObjects().store<Epetra_DistObject>(pobj, true));
}

/* store Epetra_DistObject in non-const table */
CT_Epetra_DistObject_ID_t
CEpetra::storeDistObject( Epetra_DistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
        tableOfDistObjects().store<Epetra_DistObject>(pobj, false));
}

/* store const Epetra_DistObject in const table */
CT_Epetra_DistObject_ID_t
CEpetra::storeConstDistObject( const Epetra_DistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
        tableOfDistObjects().store<Epetra_DistObject>(pobj, false));
}

/* remove Epetra_DistObject from table using CT_Epetra_DistObject_ID */
void
CEpetra::removeDistObject( CT_Epetra_DistObject_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(*id);
    if (tableOfDistObjects().isType(aid.table))
        tableOfDistObjects().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(aid);
}

/* remove Epetra_DistObject from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeDistObject( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfDistObjects().isType(aid->table))
        tableOfDistObjects().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_DistObject table */
void
CEpetra::purgeDistObject(  )
{
    tableOfDistObjects().purge();
}

/* store Epetra_DistObject in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasDistObject( const Teuchos::RCP< Epetra_DistObject > & robj )
{
    return tableOfDistObjects().alias(robj);
}

/* store const Epetra_DistObject in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstDistObject( const Teuchos::RCP< const Epetra_DistObject > & robj )
{
    return tableOfDistObjects().alias(robj);
}



