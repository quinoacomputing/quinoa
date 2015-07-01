
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
#include "CEpetra_Comm.h"
#include "CEpetra_Comm_Cpp.hpp"
#include "Epetra_Comm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Comm */
Table<Epetra_Comm>& tableOfComms()
{
    static Table<Epetra_Comm> loc_tableOfComms(CT_Epetra_Comm_ID);
    return loc_tableOfComms;
}


} // namespace


//
// Definitions from CEpetra_Comm.h
//


extern "C" {


CT_Epetra_Comm_ID_t Epetra_Comm_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Comm_Generalize ( 
  CT_Epetra_Comm_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(id);
}

void Epetra_Comm_Destroy ( CT_Epetra_Comm_ID_t * selfID )
{
    CEpetra::removeComm(selfID);
}

CT_Epetra_Comm_ID_t Epetra_Comm_Clone ( CT_Epetra_Comm_ID_t selfID )
{
    return CEpetra::storeComm(CEpetra::getConstComm(selfID)->Clone());
}

void Epetra_Comm_Barrier ( CT_Epetra_Comm_ID_t selfID )
{
    CEpetra::getConstComm(selfID)->Barrier();
}

int Epetra_Comm_Broadcast_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, int Count, int Root )
{
    return CEpetra::getConstComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_Comm_Broadcast_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int Count, int Root )
{
    return CEpetra::getConstComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_Comm_Broadcast_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, int Count, int Root )
{
    return CEpetra::getConstComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_Comm_Broadcast_Char ( 
  CT_Epetra_Comm_ID_t selfID, char * MyVals, int Count, int Root )
{
    return CEpetra::getConstComm(selfID)->Broadcast(MyVals, Count, Root);
}

int Epetra_Comm_GatherAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->GatherAll(MyVals, AllVals, Count);
}

int Epetra_Comm_GatherAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->GatherAll(MyVals, AllVals, Count);
}

int Epetra_Comm_GatherAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->GatherAll(MyVals, AllVals, Count);
}

int Epetra_Comm_SumAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count )
{
    return CEpetra::getConstComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_Comm_SumAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialSums, int * GlobalSums, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_Comm_SumAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialSums, long * GlobalSums, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->SumAll(PartialSums, GlobalSums, 
        Count);
}

int Epetra_Comm_MaxAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count )
{
    return CEpetra::getConstComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_Comm_MaxAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialMaxs, int * GlobalMaxs, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_Comm_MaxAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialMaxs, long * GlobalMaxs, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->MaxAll(PartialMaxs, GlobalMaxs, 
        Count);
}

int Epetra_Comm_MinAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count )
{
    return CEpetra::getConstComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_Comm_MinAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialMins, int * GlobalMins, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_Comm_MinAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialMins, long * GlobalMins, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->MinAll(PartialMins, GlobalMins, 
        Count);
}

int Epetra_Comm_ScanSum_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, double * ScanSums, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->ScanSum(MyVals, ScanSums, Count);
}

int Epetra_Comm_ScanSum_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->ScanSum(MyVals, ScanSums, Count);
}

int Epetra_Comm_ScanSum_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count )
{
    return CEpetra::getConstComm(selfID)->ScanSum(MyVals, ScanSums, Count);
}

int Epetra_Comm_MyPID ( CT_Epetra_Comm_ID_t selfID )
{
    return CEpetra::getConstComm(selfID)->MyPID();
}

int Epetra_Comm_NumProc ( CT_Epetra_Comm_ID_t selfID )
{
    return CEpetra::getConstComm(selfID)->NumProc();
}

CT_Epetra_Distributor_ID_t Epetra_Comm_CreateDistributor ( 
  CT_Epetra_Comm_ID_t selfID )
{
    return CEpetra::storeDistributor(CEpetra::getConstComm(
        selfID)->CreateDistributor());
}

CT_Epetra_Directory_ID_t Epetra_Comm_CreateDirectory ( 
  CT_Epetra_Comm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return CEpetra::storeDirectory(CEpetra::getConstComm(
        selfID)->CreateDirectory(*Map));
}


} // extern "C"


//
// Definitions from CEpetra_Comm_Cpp.hpp
//


/* get Epetra_Comm from non-const table using CT_Epetra_Comm_ID */
const Teuchos::RCP<Epetra_Comm>
CEpetra::getComm( CT_Epetra_Comm_ID_t id )
{
    if (tableOfComms().isType(id.table))
        return tableOfComms().get<Epetra_Comm>(
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_Comm>(
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(id));
}

/* get Epetra_Comm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Comm>
CEpetra::getComm( CTrilinos_Universal_ID_t id )
{
    if (tableOfComms().isType(id.table))
        return tableOfComms().get<Epetra_Comm>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_Comm>(id);
}

/* get const Epetra_Comm from either the const or non-const table
 * using CT_Epetra_Comm_ID */
const Teuchos::RCP<const Epetra_Comm>
CEpetra::getConstComm( CT_Epetra_Comm_ID_t id )
{
    if (tableOfComms().isType(id.table))
        return tableOfComms().getConst<Epetra_Comm>(
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_Comm>(
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(id));
}

/* get const Epetra_Comm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Comm>
CEpetra::getConstComm( CTrilinos_Universal_ID_t id )
{
    if (tableOfComms().isType(id.table))
        return tableOfComms().getConst<Epetra_Comm>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_Comm>(id);
}

/* store Epetra_Comm (owned) in non-const table */
CT_Epetra_Comm_ID_t
CEpetra::storeNewComm( Epetra_Comm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(
        tableOfComms().store<Epetra_Comm>(pobj, true));
}

/* store Epetra_Comm in non-const table */
CT_Epetra_Comm_ID_t
CEpetra::storeComm( Epetra_Comm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(
        tableOfComms().store<Epetra_Comm>(pobj, false));
}

/* store const Epetra_Comm in const table */
CT_Epetra_Comm_ID_t
CEpetra::storeConstComm( const Epetra_Comm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(
        tableOfComms().store<Epetra_Comm>(pobj, false));
}

/* remove Epetra_Comm from table using CT_Epetra_Comm_ID */
void
CEpetra::removeComm( CT_Epetra_Comm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(*id);
    if (tableOfComms().isType(aid.table))
        tableOfComms().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(aid);
}

/* remove Epetra_Comm from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeComm( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfComms().isType(aid->table))
        tableOfComms().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_Comm table */
void
CEpetra::purgeComm(  )
{
    tableOfComms().purge();
}

/* store Epetra_Comm in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasComm( const Teuchos::RCP< Epetra_Comm > & robj )
{
    return tableOfComms().alias(robj);
}

/* store const Epetra_Comm in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstComm( const Teuchos::RCP< const Epetra_Comm > & robj )
{
    return tableOfComms().alias(robj);
}



