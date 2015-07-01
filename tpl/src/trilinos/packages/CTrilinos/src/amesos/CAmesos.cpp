
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


#ifdef HAVE_CTRILINOS_AMESOS


#include "CTrilinos_enums.h"
#include "CAmesos.h"
#include "CAmesos_Cpp.hpp"
#include "Amesos.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"
#include "CAmesos_BaseSolver_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Amesos */
Table<Amesos>& tableOfAmesoss()
{
    static Table<Amesos> loc_tableOfAmesoss(CT_Amesos_ID);
    return loc_tableOfAmesoss;
}


} // namespace


//
// Definitions from CAmesos.h
//


extern "C" {


CT_Amesos_ID_t Amesos_Create (  )
{
    return CAmesos::storeNewAmesos(new Amesos());
}

void Amesos_Destroy ( CT_Amesos_ID_t * selfID )
{
    CAmesos::removeAmesos(selfID);
}

CT_Amesos_BaseSolver_ID_t Amesos_CreateSolver ( 
  CT_Amesos_ID_t selfID, const char * ClassType, 
  CT_Epetra_LinearProblem_ID_t LinearProblemID )
{
    const Teuchos::RCP<const Epetra_LinearProblem> LinearProblem = 
        CEpetra::getConstLinearProblem(LinearProblemID);
    return CAmesos::storeBaseSolver(CAmesos::getAmesos(selfID)->Create(
        ClassType, *LinearProblem));
}

boolean Amesos_Query ( 
  CT_Amesos_ID_t selfID, const char * ClassType )
{
    return ((CAmesos::getAmesos(selfID)->Query(ClassType)) ? TRUE : FALSE);
}

CT_Teuchos_ParameterList_ID_t Amesos_GetValidParameters (  )
{
    return CTeuchos::storeParameterList(new Teuchos::ParameterList(
        Amesos::GetValidParameters()));
}


} // extern "C"


//
// Definitions from CAmesos_Cpp.hpp
//


/* get Amesos from non-const table using CT_Amesos_ID */
const Teuchos::RCP<Amesos>
CAmesos::getAmesos( CT_Amesos_ID_t id )
{
    return tableOfAmesoss().get<Amesos>(
        CTrilinos::abstractType<CT_Amesos_ID_t>(id));
}

/* get Amesos from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Amesos>
CAmesos::getAmesos( CTrilinos_Universal_ID_t id )
{
    return tableOfAmesoss().get<Amesos>(id);
}

/* get const Amesos from either the const or non-const table
 * using CT_Amesos_ID */
const Teuchos::RCP<const Amesos>
CAmesos::getConstAmesos( CT_Amesos_ID_t id )
{
    return tableOfAmesoss().getConst<Amesos>(
        CTrilinos::abstractType<CT_Amesos_ID_t>(id));
}

/* get const Amesos from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Amesos>
CAmesos::getConstAmesos( CTrilinos_Universal_ID_t id )
{
    return tableOfAmesoss().getConst<Amesos>(id);
}

/* store Amesos (owned) in non-const table */
CT_Amesos_ID_t
CAmesos::storeNewAmesos( Amesos *pobj )
{
    return CTrilinos::concreteType<CT_Amesos_ID_t>(
        tableOfAmesoss().store<Amesos>(pobj, true));
}

/* store Amesos in non-const table */
CT_Amesos_ID_t
CAmesos::storeAmesos( Amesos *pobj )
{
    return CTrilinos::concreteType<CT_Amesos_ID_t>(
        tableOfAmesoss().store<Amesos>(pobj, false));
}

/* store const Amesos in const table */
CT_Amesos_ID_t
CAmesos::storeConstAmesos( const Amesos *pobj )
{
    return CTrilinos::concreteType<CT_Amesos_ID_t>(
        tableOfAmesoss().store<Amesos>(pobj, false));
}

/* remove Amesos from table using CT_Amesos_ID */
void
CAmesos::removeAmesos( CT_Amesos_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Amesos_ID_t>(*id);
    tableOfAmesoss().remove(&aid);
    *id = CTrilinos::concreteType<CT_Amesos_ID_t>(aid);
}

/* remove Amesos from table using CTrilinos_Universal_ID_t */
void
CAmesos::removeAmesos( CTrilinos_Universal_ID_t *aid )
{
    tableOfAmesoss().remove(aid);
}

/* purge Amesos table */
void
CAmesos::purgeAmesos(  )
{
    tableOfAmesoss().purge();
}



#endif /* HAVE_CTRILINOS_AMESOS */


