
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
#include "CAztecOO_StatusTest.h"
#include "CAztecOO_StatusTest_Cpp.hpp"
#include "AztecOO_StatusTest.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type AztecOO_StatusTest */
Table<AztecOO_StatusTest>& tableOfStatusTests()
{
    static Table<AztecOO_StatusTest> loc_tableOfStatusTests(CT_AztecOO_StatusTest_ID);
    return loc_tableOfStatusTests;
}


} // namespace


//
// Definitions from CAztecOO_StatusTest.h
//


extern "C" {


CT_AztecOO_StatusTest_ID_t AztecOO_StatusTest_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(id);
}

CTrilinos_Universal_ID_t AztecOO_StatusTest_Generalize ( 
  CT_AztecOO_StatusTest_ID_t id )
{
    return CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id);
}

void AztecOO_StatusTest_Destroy ( 
  CT_AztecOO_StatusTest_ID_t * selfID )
{
    CAztecOO::removeStatusTest(selfID);
}

boolean AztecOO_StatusTest_ResidualVectorRequired ( 
  CT_AztecOO_StatusTest_ID_t selfID )
{
    return ((CAztecOO::getConstStatusTest(
        selfID)->ResidualVectorRequired()) ? TRUE : FALSE);
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTest_CheckStatus ( 
  CT_AztecOO_StatusTest_ID_t selfID, int CurrentIter, 
  CT_Epetra_MultiVector_ID_t CurrentResVectorID, 
  double CurrentResNormEst, boolean SolutionUpdated )
{
    const Teuchos::RCP<Epetra_MultiVector> CurrentResVector = 
        CEpetra::getMultiVector(CurrentResVectorID);
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getStatusTest(
        selfID)->CheckStatus(CurrentIter, CurrentResVector.getRawPtr(), 
        CurrentResNormEst, ((SolutionUpdated) != FALSE ? true : false)) );
}

CT_AztecOO_StatusType_E_t AztecOO_StatusTest_GetStatus ( 
  CT_AztecOO_StatusTest_ID_t selfID )
{
    return (CT_AztecOO_StatusType_E_t)( CAztecOO::getConstStatusTest(
        selfID)->GetStatus() );
}


} // extern "C"


//
// Definitions from CAztecOO_StatusTest_Cpp.hpp
//


/* get AztecOO_StatusTest from non-const table using CT_AztecOO_StatusTest_ID */
const Teuchos::RCP<AztecOO_StatusTest>
CAztecOO::getStatusTest( CT_AztecOO_StatusTest_ID_t id )
{
    if (tableOfStatusTests().isType(id.table))
        return tableOfStatusTests().get<AztecOO_StatusTest>(
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<AztecOO_StatusTest>(
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id));
}

/* get AztecOO_StatusTest from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTest>
CAztecOO::getStatusTest( CTrilinos_Universal_ID_t id )
{
    if (tableOfStatusTests().isType(id.table))
        return tableOfStatusTests().get<AztecOO_StatusTest>(id);
    else
        return CTrilinos::TableRepos::get<AztecOO_StatusTest>(id);
}

/* get const AztecOO_StatusTest from either the const or non-const table
 * using CT_AztecOO_StatusTest_ID */
const Teuchos::RCP<const AztecOO_StatusTest>
CAztecOO::getConstStatusTest( CT_AztecOO_StatusTest_ID_t id )
{
    if (tableOfStatusTests().isType(id.table))
        return tableOfStatusTests().getConst<AztecOO_StatusTest>(
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<AztecOO_StatusTest>(
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id));
}

/* get const AztecOO_StatusTest from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTest>
CAztecOO::getConstStatusTest( CTrilinos_Universal_ID_t id )
{
    if (tableOfStatusTests().isType(id.table))
        return tableOfStatusTests().getConst<AztecOO_StatusTest>(id);
    else
        return CTrilinos::TableRepos::getConst<AztecOO_StatusTest>(id);
}

/* store AztecOO_StatusTest (owned) in non-const table */
CT_AztecOO_StatusTest_ID_t
CAztecOO::storeNewStatusTest( AztecOO_StatusTest *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
        tableOfStatusTests().store<AztecOO_StatusTest>(pobj, true));
}

/* store AztecOO_StatusTest in non-const table */
CT_AztecOO_StatusTest_ID_t
CAztecOO::storeStatusTest( AztecOO_StatusTest *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
        tableOfStatusTests().store<AztecOO_StatusTest>(pobj, false));
}

/* store const AztecOO_StatusTest in const table */
CT_AztecOO_StatusTest_ID_t
CAztecOO::storeConstStatusTest( const AztecOO_StatusTest *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
        tableOfStatusTests().store<AztecOO_StatusTest>(pobj, false));
}

/* remove AztecOO_StatusTest from table using CT_AztecOO_StatusTest_ID */
void
CAztecOO::removeStatusTest( CT_AztecOO_StatusTest_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(*id);
    if (tableOfStatusTests().isType(aid.table))
        tableOfStatusTests().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(aid);
}

/* remove AztecOO_StatusTest from table using CTrilinos_Universal_ID_t */
void
CAztecOO::removeStatusTest( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfStatusTests().isType(aid->table))
        tableOfStatusTests().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge AztecOO_StatusTest table */
void
CAztecOO::purgeStatusTest(  )
{
    tableOfStatusTests().purge();
}

/* store AztecOO_StatusTest in non-const table */
CTrilinos_Universal_ID_t
CAztecOO::aliasStatusTest( const Teuchos::RCP< AztecOO_StatusTest > & robj )
{
    return tableOfStatusTests().alias(robj);
}

/* store const AztecOO_StatusTest in const table */
CTrilinos_Universal_ID_t
CAztecOO::aliasConstStatusTest( const Teuchos::RCP< const AztecOO_StatusTest > & robj )
{
    return tableOfStatusTests().alias(robj);
}



#endif /* HAVE_CTRILINOS_AZTECOO */


