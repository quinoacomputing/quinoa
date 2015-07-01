
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
#include "CEpetra_Time.h"
#include "CEpetra_Time_Cpp.hpp"
#include "Epetra_Time.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_Comm_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Time */
Table<Epetra_Time>& tableOfTimes()
{
    static Table<Epetra_Time> loc_tableOfTimes(CT_Epetra_Time_ID);
    return loc_tableOfTimes;
}


} // namespace


//
// Definitions from CEpetra_Time.h
//


extern "C" {


CT_Epetra_Time_ID_t Epetra_Time_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Time_Generalize ( 
  CT_Epetra_Time_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id);
}

CT_Epetra_Time_ID_t Epetra_Time_Create ( CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = CEpetra::getConstComm(CommID);
    return CEpetra::storeNewTime(new Epetra_Time(*Comm));
}

CT_Epetra_Time_ID_t Epetra_Time_Duplicate ( 
  CT_Epetra_Time_ID_t TimeID )
{
    const Teuchos::RCP<const Epetra_Time> Time = CEpetra::getConstTime(TimeID);
    return CEpetra::storeNewTime(new Epetra_Time(*Time));
}

void Epetra_Time_Destroy ( CT_Epetra_Time_ID_t * selfID )
{
    CEpetra::removeTime(selfID);
}

double Epetra_Time_WallTime ( CT_Epetra_Time_ID_t selfID )
{
    return CEpetra::getConstTime(selfID)->WallTime();
}

void Epetra_Time_ResetStartTime ( CT_Epetra_Time_ID_t selfID )
{
    CEpetra::getTime(selfID)->ResetStartTime();
}

double Epetra_Time_ElapsedTime ( CT_Epetra_Time_ID_t selfID )
{
    return CEpetra::getConstTime(selfID)->ElapsedTime();
}

void Epetra_Time_Assign ( 
  CT_Epetra_Time_ID_t selfID, CT_Epetra_Time_ID_t srcID )
{
    Epetra_Time& self = *( CEpetra::getTime(selfID) );

    const Teuchos::RCP<const Epetra_Time> src = CEpetra::getConstTime(srcID);
    self = *src;
}


} // extern "C"


//
// Definitions from CEpetra_Time_Cpp.hpp
//


/* get Epetra_Time from non-const table using CT_Epetra_Time_ID */
const Teuchos::RCP<Epetra_Time>
CEpetra::getTime( CT_Epetra_Time_ID_t id )
{
    if (tableOfTimes().isType(id.table))
        return tableOfTimes().get<Epetra_Time>(
        CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_Time>(
        CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id));
}

/* get Epetra_Time from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Time>
CEpetra::getTime( CTrilinos_Universal_ID_t id )
{
    if (tableOfTimes().isType(id.table))
        return tableOfTimes().get<Epetra_Time>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_Time>(id);
}

/* get const Epetra_Time from either the const or non-const table
 * using CT_Epetra_Time_ID */
const Teuchos::RCP<const Epetra_Time>
CEpetra::getConstTime( CT_Epetra_Time_ID_t id )
{
    if (tableOfTimes().isType(id.table))
        return tableOfTimes().getConst<Epetra_Time>(
        CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_Time>(
        CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id));
}

/* get const Epetra_Time from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Time>
CEpetra::getConstTime( CTrilinos_Universal_ID_t id )
{
    if (tableOfTimes().isType(id.table))
        return tableOfTimes().getConst<Epetra_Time>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_Time>(id);
}

/* store Epetra_Time (owned) in non-const table */
CT_Epetra_Time_ID_t
CEpetra::storeNewTime( Epetra_Time *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
        tableOfTimes().store<Epetra_Time>(pobj, true));
}

/* store Epetra_Time in non-const table */
CT_Epetra_Time_ID_t
CEpetra::storeTime( Epetra_Time *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
        tableOfTimes().store<Epetra_Time>(pobj, false));
}

/* store const Epetra_Time in const table */
CT_Epetra_Time_ID_t
CEpetra::storeConstTime( const Epetra_Time *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
        tableOfTimes().store<Epetra_Time>(pobj, false));
}

/* remove Epetra_Time from table using CT_Epetra_Time_ID */
void
CEpetra::removeTime( CT_Epetra_Time_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Time_ID_t>(*id);
    if (tableOfTimes().isType(aid.table))
        tableOfTimes().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Time_ID_t>(aid);
}

/* remove Epetra_Time from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeTime( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfTimes().isType(aid->table))
        tableOfTimes().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_Time table */
void
CEpetra::purgeTime(  )
{
    tableOfTimes().purge();
}

/* store Epetra_Time in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasTime( const Teuchos::RCP< Epetra_Time > & robj )
{
    return tableOfTimes().alias(robj);
}

/* store const Epetra_Time in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstTime( const Teuchos::RCP< const Epetra_Time > & robj )
{
    return tableOfTimes().alias(robj);
}



