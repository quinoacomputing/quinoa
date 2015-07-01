
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
#include "CEpetra_SerialComm.h"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "Epetra_SerialComm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_SerialComm */
Table<Epetra_SerialComm>& tableOfSerialComms()
{
    static Table<Epetra_SerialComm> loc_tableOfSerialComms(CT_Epetra_SerialComm_ID);
    return loc_tableOfSerialComms;
}


} // namespace


//
// Definitions from CEpetra_SerialComm.h
//


extern "C" {


CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_SerialComm_Generalize ( 
  CT_Epetra_SerialComm_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id);
}

CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Create (  )
{
    return CEpetra::storeNewSerialComm(new Epetra_SerialComm());
}

CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Duplicate ( 
  CT_Epetra_SerialComm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_SerialComm> Comm = 
        CEpetra::getConstSerialComm(CommID);
    return CEpetra::storeNewSerialComm(new Epetra_SerialComm(*Comm));
}

void Epetra_SerialComm_Destroy ( CT_Epetra_SerialComm_ID_t * selfID )
{
    CEpetra::removeSerialComm(selfID);
}

CT_Epetra_Comm_ID_t Epetra_SerialComm_Clone ( 
  CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::storeComm(CEpetra::getConstSerialComm(selfID)->Clone());
}

void Epetra_SerialComm_Barrier ( CT_Epetra_SerialComm_ID_t selfID )
{
    CEpetra::getConstSerialComm(selfID)->Barrier();
}

int Epetra_SerialComm_Broadcast_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstSerialComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_SerialComm_Broadcast_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstSerialComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_SerialComm_Broadcast_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstSerialComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_SerialComm_Broadcast_Char ( 
  CT_Epetra_SerialComm_ID_t selfID, char * MyVals, int Count, 
  int Root )
{
    return CEpetra::getConstSerialComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_SerialComm_GatherAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  double * AllVals, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->GatherAll(MyVals, AllVals, 
        Count);
}

int Epetra_SerialComm_GatherAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->GatherAll(MyVals, AllVals, 
        Count);
}

int Epetra_SerialComm_GatherAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->GatherAll(MyVals, AllVals, 
        Count);
}

int Epetra_SerialComm_SumAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_SerialComm_SumAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialSums, 
  int * GlobalSums, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_SerialComm_SumAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_SerialComm_MaxAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_SerialComm_MaxAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMaxs, 
  int * GlobalMaxs, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_SerialComm_MaxAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_SerialComm_MinAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_SerialComm_MinAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMins, 
  int * GlobalMins, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_SerialComm_MinAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_SerialComm_ScanSum_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  double * ScanSums, int Count )
{
    return CEpetra::getConstSerialComm(selfID)->ScanSum(MyVals, ScanSums, 
        Count);
}

int Epetra_SerialComm_ScanSum_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->ScanSum(MyVals, ScanSums, 
        Count);
}

int Epetra_SerialComm_ScanSum_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count )
{
    return CEpetra::getConstSerialComm(selfID)->ScanSum(MyVals, ScanSums, 
        Count);
}

int Epetra_SerialComm_MyPID ( CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::getConstSerialComm(selfID)->MyPID();
}

int Epetra_SerialComm_NumProc ( CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::getConstSerialComm(selfID)->NumProc();
}

CT_Epetra_Distributor_ID_t Epetra_SerialComm_CreateDistributor ( 
  CT_Epetra_SerialComm_ID_t selfID )
{
    return CEpetra::storeDistributor(CEpetra::getConstSerialComm(
        selfID)->CreateDistributor());
}

CT_Epetra_Directory_ID_t Epetra_SerialComm_CreateDirectory ( 
  CT_Epetra_SerialComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return CEpetra::storeDirectory(CEpetra::getConstSerialComm(
        selfID)->CreateDirectory(*Map));
}

void Epetra_SerialComm_Assign ( 
  CT_Epetra_SerialComm_ID_t selfID, 
  CT_Epetra_SerialComm_ID_t CommID )
{
    Epetra_SerialComm& self = *( CEpetra::getSerialComm(selfID) );

    const Teuchos::RCP<const Epetra_SerialComm> Comm = 
        CEpetra::getConstSerialComm(CommID);
    self = *Comm;
}


} // extern "C"


//
// Definitions from CEpetra_SerialComm_Cpp.hpp
//


/* get Epetra_SerialComm from non-const table using CT_Epetra_SerialComm_ID */
const Teuchos::RCP<Epetra_SerialComm>
CEpetra::getSerialComm( CT_Epetra_SerialComm_ID_t id )
{
    if (tableOfSerialComms().isType(id.table))
        return tableOfSerialComms().get<Epetra_SerialComm>(
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_SerialComm>(
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id));
}

/* get Epetra_SerialComm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_SerialComm>
CEpetra::getSerialComm( CTrilinos_Universal_ID_t id )
{
    if (tableOfSerialComms().isType(id.table))
        return tableOfSerialComms().get<Epetra_SerialComm>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_SerialComm>(id);
}

/* get const Epetra_SerialComm from either the const or non-const table
 * using CT_Epetra_SerialComm_ID */
const Teuchos::RCP<const Epetra_SerialComm>
CEpetra::getConstSerialComm( CT_Epetra_SerialComm_ID_t id )
{
    if (tableOfSerialComms().isType(id.table))
        return tableOfSerialComms().getConst<Epetra_SerialComm>(
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_SerialComm>(
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id));
}

/* get const Epetra_SerialComm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_SerialComm>
CEpetra::getConstSerialComm( CTrilinos_Universal_ID_t id )
{
    if (tableOfSerialComms().isType(id.table))
        return tableOfSerialComms().getConst<Epetra_SerialComm>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_SerialComm>(id);
}

/* store Epetra_SerialComm (owned) in non-const table */
CT_Epetra_SerialComm_ID_t
CEpetra::storeNewSerialComm( Epetra_SerialComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
        tableOfSerialComms().store<Epetra_SerialComm>(pobj, true));
}

/* store Epetra_SerialComm in non-const table */
CT_Epetra_SerialComm_ID_t
CEpetra::storeSerialComm( Epetra_SerialComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
        tableOfSerialComms().store<Epetra_SerialComm>(pobj, false));
}

/* store const Epetra_SerialComm in const table */
CT_Epetra_SerialComm_ID_t
CEpetra::storeConstSerialComm( const Epetra_SerialComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
        tableOfSerialComms().store<Epetra_SerialComm>(pobj, false));
}

/* remove Epetra_SerialComm from table using CT_Epetra_SerialComm_ID */
void
CEpetra::removeSerialComm( CT_Epetra_SerialComm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(*id);
    if (tableOfSerialComms().isType(aid.table))
        tableOfSerialComms().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(aid);
}

/* remove Epetra_SerialComm from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeSerialComm( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfSerialComms().isType(aid->table))
        tableOfSerialComms().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_SerialComm table */
void
CEpetra::purgeSerialComm(  )
{
    tableOfSerialComms().purge();
}

/* store Epetra_SerialComm in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasSerialComm( const Teuchos::RCP< Epetra_SerialComm > & robj )
{
    return tableOfSerialComms().alias(robj);
}

/* store const Epetra_SerialComm in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstSerialComm( const Teuchos::RCP< const Epetra_SerialComm > & robj )
{
    return tableOfSerialComms().alias(robj);
}



