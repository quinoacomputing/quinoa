
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


#ifdef HAVE_CTRILINOS_AZTECOO


#include "CTrilinos_enums.h"
#include "CAztecOO_StatusTestMaxIters.h"
#include "CAztecOO_StatusTestMaxIters_Cpp.hpp"
#include "AztecOO_StatusTestMaxIters.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type AztecOO_StatusTestMaxIters */
Table<AztecOO_StatusTestMaxIters>& tableOfStatusTestMaxIterss()
{
    static Table<AztecOO_StatusTestMaxIters> loc_tableOfStatusTestMaxIterss(CT_AztecOO_StatusTestMaxIters_ID);
    return loc_tableOfStatusTestMaxIterss;
}


} // namespace


//
// Definitions from CAztecOO_StatusTestMaxIters.h
//


extern "C" {


CT_AztecOO_StatusTestMaxIters_ID_t AztecOO_StatusTestMaxIters_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(id);
}

CTrilinos_Universal_ID_t AztecOO_StatusTestMaxIters_Generalize ( 
  CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    return CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id);
}

CT_AztecOO_StatusTestMaxIters_ID_t AztecOO_StatusTestMaxIters_Create ( 
  int MaxIters )
{
    return CAztecOO::storeNewStatusTestMaxIters(new AztecOO_StatusTestMaxIters(
        MaxIters));
}

void AztecOO_StatusTestMaxIters_Destroy ( 
  CT_AztecOO_StatusTestMaxIters_ID_t * selfID )
{
    CAztecOO::removeStatusTestMaxIters(selfID);
}

boolean AztecOO_StatusTestMaxIters_ResidualVectorRequired ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID )
{
    return ((CAztecOO::getConstStatusTestMaxIters(
        selfID)->ResidualVectorRequired()) ? TRUE : FALSE);
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestMaxIters_CheckStatus ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID, int CurrentIter, 
  CT_Epetra_MultiVector_ID_t CurrentResVectorID, 
  double CurrentResNormEst, boolean SolutionUpdated )
{
    const Teuchos::RCP<Epetra_MultiVector> CurrentResVector = 
        CEpetra::getMultiVector(CurrentResVectorID);
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getStatusTestMaxIters(
        selfID)->CheckStatus(CurrentIter, CurrentResVector.getRawPtr(), 
        CurrentResNormEst, ((SolutionUpdated) != FALSE ? true : false)) );
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTestMaxIters_GetStatus ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID )
{
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getConstStatusTestMaxIters(
        selfID)->GetStatus() );
}

int AztecOO_StatusTestMaxIters_GetMaxIters ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID )
{
    return CAztecOO::getConstStatusTestMaxIters(selfID)->GetMaxIters();
}

int AztecOO_StatusTestMaxIters_GetNumIters ( 
  CT_AztecOO_StatusTestMaxIters_ID_t selfID )
{
    return CAztecOO::getConstStatusTestMaxIters(selfID)->GetNumIters();
}


} // extern "C"


//
// Definitions from CAztecOO_StatusTestMaxIters_Cpp.hpp
//


/* get AztecOO_StatusTestMaxIters from non-const table using CT_AztecOO_StatusTestMaxIters_ID */
const Teuchos::RCP<AztecOO_StatusTestMaxIters>
CAztecOO::getStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    if (tableOfStatusTestMaxIterss().isType(id.table))
        return tableOfStatusTestMaxIterss().get<AztecOO_StatusTestMaxIters>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<AztecOO_StatusTestMaxIters>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id));
}

/* get AztecOO_StatusTestMaxIters from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTestMaxIters>
CAztecOO::getStatusTestMaxIters( CTrilinos_Universal_ID_t id )
{
    if (tableOfStatusTestMaxIterss().isType(id.table))
        return tableOfStatusTestMaxIterss().get<AztecOO_StatusTestMaxIters>(id);
    else
        return CTrilinos::TableRepos::get<AztecOO_StatusTestMaxIters>(id);
}

/* get const AztecOO_StatusTestMaxIters from either the const or non-const table
 * using CT_AztecOO_StatusTestMaxIters_ID */
const Teuchos::RCP<const AztecOO_StatusTestMaxIters>
CAztecOO::getConstStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    if (tableOfStatusTestMaxIterss().isType(id.table))
        return tableOfStatusTestMaxIterss().getConst<AztecOO_StatusTestMaxIters>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<AztecOO_StatusTestMaxIters>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id));
}

/* get const AztecOO_StatusTestMaxIters from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTestMaxIters>
CAztecOO::getConstStatusTestMaxIters( CTrilinos_Universal_ID_t id )
{
    if (tableOfStatusTestMaxIterss().isType(id.table))
        return tableOfStatusTestMaxIterss().getConst<AztecOO_StatusTestMaxIters>(id);
    else
        return CTrilinos::TableRepos::getConst<AztecOO_StatusTestMaxIters>(id);
}

/* store AztecOO_StatusTestMaxIters (owned) in non-const table */
CT_AztecOO_StatusTestMaxIters_ID_t
CAztecOO::storeNewStatusTestMaxIters( AztecOO_StatusTestMaxIters *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        tableOfStatusTestMaxIterss().store<AztecOO_StatusTestMaxIters>(pobj, true));
}

/* store AztecOO_StatusTestMaxIters in non-const table */
CT_AztecOO_StatusTestMaxIters_ID_t
CAztecOO::storeStatusTestMaxIters( AztecOO_StatusTestMaxIters *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        tableOfStatusTestMaxIterss().store<AztecOO_StatusTestMaxIters>(pobj, false));
}

/* store const AztecOO_StatusTestMaxIters in const table */
CT_AztecOO_StatusTestMaxIters_ID_t
CAztecOO::storeConstStatusTestMaxIters( const AztecOO_StatusTestMaxIters *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        tableOfStatusTestMaxIterss().store<AztecOO_StatusTestMaxIters>(pobj, false));
}

/* remove AztecOO_StatusTestMaxIters from table using CT_AztecOO_StatusTestMaxIters_ID */
void
CAztecOO::removeStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(*id);
    if (tableOfStatusTestMaxIterss().isType(aid.table))
        tableOfStatusTestMaxIterss().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(aid);
}

/* remove AztecOO_StatusTestMaxIters from table using CTrilinos_Universal_ID_t */
void
CAztecOO::removeStatusTestMaxIters( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfStatusTestMaxIterss().isType(aid->table))
        tableOfStatusTestMaxIterss().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge AztecOO_StatusTestMaxIters table */
void
CAztecOO::purgeStatusTestMaxIters(  )
{
    tableOfStatusTestMaxIterss().purge();
}

/* store AztecOO_StatusTestMaxIters in non-const table */
CTrilinos_Universal_ID_t
CAztecOO::aliasStatusTestMaxIters( const Teuchos::RCP< AztecOO_StatusTestMaxIters > & robj )
{
    return tableOfStatusTestMaxIterss().alias(robj);
}

/* store const AztecOO_StatusTestMaxIters in const table */
CTrilinos_Universal_ID_t
CAztecOO::aliasConstStatusTestMaxIters( const Teuchos::RCP< const AztecOO_StatusTestMaxIters > & robj )
{
    return tableOfStatusTestMaxIterss().alias(robj);
}



#endif /* HAVE_CTRILINOS_AZTECOO */


