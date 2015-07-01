
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
#include "CAmesos_BaseSolver.h"
#include "CAmesos_BaseSolver_Cpp.hpp"
#include "Amesos_BaseSolver.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Amesos_BaseSolver */
Table<Amesos_BaseSolver>& tableOfBaseSolvers()
{
    static Table<Amesos_BaseSolver> loc_tableOfBaseSolvers(CT_Amesos_BaseSolver_ID);
    return loc_tableOfBaseSolvers;
}


} // namespace


//
// Definitions from CAmesos_BaseSolver.h
//


extern "C" {


CT_Amesos_BaseSolver_ID_t Amesos_BaseSolver_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(id);
}

CTrilinos_Universal_ID_t Amesos_BaseSolver_Generalize ( 
  CT_Amesos_BaseSolver_ID_t id )
{
    return CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id);
}

void Amesos_BaseSolver_Destroy ( CT_Amesos_BaseSolver_ID_t * selfID )
{
    CAmesos::removeBaseSolver(selfID);
}

int Amesos_BaseSolver_SymbolicFactorization ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getBaseSolver(selfID)->SymbolicFactorization();
}

int Amesos_BaseSolver_NumericFactorization ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getBaseSolver(selfID)->NumericFactorization();
}

int Amesos_BaseSolver_Solve ( CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getBaseSolver(selfID)->Solve();
}

int Amesos_BaseSolver_SetUseTranspose ( 
  CT_Amesos_BaseSolver_ID_t selfID, boolean UseTranspose )
{
    return CAmesos::getBaseSolver(selfID)->SetUseTranspose(((UseTranspose) != 
        FALSE ? true : false));
}

boolean Amesos_BaseSolver_UseTranspose ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return ((CAmesos::getConstBaseSolver(
        selfID)->UseTranspose()) ? TRUE : FALSE);
}

int Amesos_BaseSolver_SetParameters ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t ParameterListID )
{
    const Teuchos::RCP<Teuchos::ParameterList> ParameterList = 
        CTeuchos::getParameterList(ParameterListID);
    return CAmesos::getBaseSolver(selfID)->SetParameters(*ParameterList);
}

CT_Epetra_LinearProblem_ID_t Amesos_BaseSolver_GetProblem ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CEpetra::storeConstLinearProblem(CAmesos::getConstBaseSolver(
        selfID)->GetProblem());
}

boolean Amesos_BaseSolver_MatrixShapeOK ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return ((CAmesos::getConstBaseSolver(
        selfID)->MatrixShapeOK()) ? TRUE : FALSE);
}

CT_Epetra_Comm_ID_t Amesos_BaseSolver_Comm ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CEpetra::storeConstComm(&( CAmesos::getConstBaseSolver(
        selfID)->Comm() ));
}

int Amesos_BaseSolver_NumSymbolicFact ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getConstBaseSolver(selfID)->NumSymbolicFact();
}

int Amesos_BaseSolver_NumNumericFact ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getConstBaseSolver(selfID)->NumNumericFact();
}

int Amesos_BaseSolver_NumSolve ( CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getConstBaseSolver(selfID)->NumSolve();
}

void Amesos_BaseSolver_PrintStatus ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    CAmesos::getConstBaseSolver(selfID)->PrintStatus();
}

void Amesos_BaseSolver_PrintTiming ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    CAmesos::getConstBaseSolver(selfID)->PrintTiming();
}

void Amesos_BaseSolver_setParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t paramListID )
{
    const Teuchos::RCP<Teuchos::ParameterList> paramList = 
        CTeuchos::getParameterList(paramListID);
    CAmesos::getBaseSolver(selfID)->setParameterList(paramList);
}

CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_getNonconstParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CTeuchos::storeParameterList(CAmesos::getBaseSolver(
        selfID)->getNonconstParameterList().getRawPtr());
}

CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_unsetParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CTeuchos::storeParameterList(CAmesos::getBaseSolver(
        selfID)->unsetParameterList().getRawPtr());
}

void Amesos_BaseSolver_GetTiming ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t TimingParameterListID )
{
    const Teuchos::RCP<Teuchos::ParameterList> TimingParameterList = 
        CTeuchos::getParameterList(TimingParameterListID);
    CAmesos::getConstBaseSolver(selfID)->GetTiming(*TimingParameterList);
}


} // extern "C"


//
// Definitions from CAmesos_BaseSolver_Cpp.hpp
//


/* get Amesos_BaseSolver from non-const table using CT_Amesos_BaseSolver_ID */
const Teuchos::RCP<Amesos_BaseSolver>
CAmesos::getBaseSolver( CT_Amesos_BaseSolver_ID_t id )
{
    if (tableOfBaseSolvers().isType(id.table))
        return tableOfBaseSolvers().get<Amesos_BaseSolver>(
        CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Amesos_BaseSolver>(
        CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id));
}

/* get Amesos_BaseSolver from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Amesos_BaseSolver>
CAmesos::getBaseSolver( CTrilinos_Universal_ID_t id )
{
    if (tableOfBaseSolvers().isType(id.table))
        return tableOfBaseSolvers().get<Amesos_BaseSolver>(id);
    else
        return CTrilinos::TableRepos::get<Amesos_BaseSolver>(id);
}

/* get const Amesos_BaseSolver from either the const or non-const table
 * using CT_Amesos_BaseSolver_ID */
const Teuchos::RCP<const Amesos_BaseSolver>
CAmesos::getConstBaseSolver( CT_Amesos_BaseSolver_ID_t id )
{
    if (tableOfBaseSolvers().isType(id.table))
        return tableOfBaseSolvers().getConst<Amesos_BaseSolver>(
        CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Amesos_BaseSolver>(
        CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id));
}

/* get const Amesos_BaseSolver from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Amesos_BaseSolver>
CAmesos::getConstBaseSolver( CTrilinos_Universal_ID_t id )
{
    if (tableOfBaseSolvers().isType(id.table))
        return tableOfBaseSolvers().getConst<Amesos_BaseSolver>(id);
    else
        return CTrilinos::TableRepos::getConst<Amesos_BaseSolver>(id);
}

/* store Amesos_BaseSolver (owned) in non-const table */
CT_Amesos_BaseSolver_ID_t
CAmesos::storeNewBaseSolver( Amesos_BaseSolver *pobj )
{
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(
        tableOfBaseSolvers().store<Amesos_BaseSolver>(pobj, true));
}

/* store Amesos_BaseSolver in non-const table */
CT_Amesos_BaseSolver_ID_t
CAmesos::storeBaseSolver( Amesos_BaseSolver *pobj )
{
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(
        tableOfBaseSolvers().store<Amesos_BaseSolver>(pobj, false));
}

/* store const Amesos_BaseSolver in const table */
CT_Amesos_BaseSolver_ID_t
CAmesos::storeConstBaseSolver( const Amesos_BaseSolver *pobj )
{
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(
        tableOfBaseSolvers().store<Amesos_BaseSolver>(pobj, false));
}

/* remove Amesos_BaseSolver from table using CT_Amesos_BaseSolver_ID */
void
CAmesos::removeBaseSolver( CT_Amesos_BaseSolver_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(*id);
    if (tableOfBaseSolvers().isType(aid.table))
        tableOfBaseSolvers().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(aid);
}

/* remove Amesos_BaseSolver from table using CTrilinos_Universal_ID_t */
void
CAmesos::removeBaseSolver( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfBaseSolvers().isType(aid->table))
        tableOfBaseSolvers().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Amesos_BaseSolver table */
void
CAmesos::purgeBaseSolver(  )
{
    tableOfBaseSolvers().purge();
}

/* store Amesos_BaseSolver in non-const table */
CTrilinos_Universal_ID_t
CAmesos::aliasBaseSolver( const Teuchos::RCP< Amesos_BaseSolver > & robj )
{
    return tableOfBaseSolvers().alias(robj);
}

/* store const Amesos_BaseSolver in const table */
CTrilinos_Universal_ID_t
CAmesos::aliasConstBaseSolver( const Teuchos::RCP< const Amesos_BaseSolver > & robj )
{
    return tableOfBaseSolvers().alias(robj);
}



#endif /* HAVE_CTRILINOS_AMESOS */


