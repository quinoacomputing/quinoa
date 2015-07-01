
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
#include "CTeuchos_ParameterList.h"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CTeuchos_ParameterEntry_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Teuchos::ParameterList */
Table<Teuchos::ParameterList>& tableOfParameterLists()
{
    static Table<Teuchos::ParameterList> loc_tableOfParameterLists(CT_Teuchos_ParameterList_ID);
    return loc_tableOfParameterLists;
}


} // namespace


//
// Definitions from CTeuchos_ParameterList.h
//


extern "C" {


CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(id);
}

CTrilinos_Universal_ID_t Teuchos_ParameterList_Generalize ( 
  CT_Teuchos_ParameterList_ID_t id )
{
    return CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id);
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create (  )
{
    return CTeuchos::storeNewParameterList(new Teuchos::ParameterList());
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_WithName ( 
  const char name[] )
{
    return CTeuchos::storeNewParameterList(new Teuchos::ParameterList(
        std::string(name)));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_FromSource ( 
  CT_Teuchos_ParameterList_ID_t sourceID )
{
    const Teuchos::RCP<const Teuchos::ParameterList> source = 
        CTeuchos::getConstParameterList(sourceID);
    return CTeuchos::storeNewParameterList(new Teuchos::ParameterList(
        *source));
}

void Teuchos_ParameterList_Destroy ( 
  CT_Teuchos_ParameterList_ID_t * selfID )
{
    CTeuchos::removeParameterList(selfID);
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setName ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->setName(std::string(name)) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID )
{
    const Teuchos::RCP<const Teuchos::ParameterList> source = 
        CTeuchos::getConstParameterList(sourceID);
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->setParameters(*source) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParametersNotAlreadySet ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID )
{
    const Teuchos::RCP<const Teuchos::ParameterList> source = 
        CTeuchos::getConstParameterList(sourceID);
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->setParametersNotAlreadySet(*source) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_disableRecursiveValidation ( 
  CT_Teuchos_ParameterList_ID_t selfID )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->disableRecursiveValidation() ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  double value, char const docString[] )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->set<double>(std::string(name), value, std::string(docString), 
         Teuchos::null ) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  int value, char const docString[] )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->set<int>(std::string(name), value, std::string(docString), 
         Teuchos::null ) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_str ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  char value[], char const docString[] )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->set(std::string(name), value, std::string(docString), 
         Teuchos::null ) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  CT_Teuchos_ParameterList_ID_t valueID, char const docString[] )
{
    const Teuchos::RCP<const Teuchos::ParameterList> value = 
        CTeuchos::getConstParameterList(valueID);
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->set(std::string(name), *value, std::string(docString)) ));
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  CT_Teuchos_ParameterEntry_ID_t entryID )
{
    const Teuchos::RCP<const Teuchos::ParameterEntry> entry = 
        CTeuchos::getConstParameterEntry(entryID);
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->setEntry(std::string(name), *entry) ));
}

double Teuchos_ParameterList_get_def_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double def_value )
{
    return CTeuchos::getParameterList(selfID)->get<double>(std::string(name), 
        def_value);
}

int Teuchos_ParameterList_get_def_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int def_value )
{
    return CTeuchos::getParameterList(selfID)->get<int>(std::string(name), 
        def_value);
}

const char * Teuchos_ParameterList_get_def_char ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  char def_value[] )
{
    return CTeuchos::getParameterList(selfID)->get(std::string(name), 
        def_value).c_str();
}

const char * Teuchos_ParameterList_get_def_const_char ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  const char def_value[] )
{
    return CTeuchos::getParameterList(selfID)->get(std::string(name), 
        def_value).c_str();
}

double Teuchos_ParameterList_get_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getParameterList(selfID)->get<double>(std::string(name));
}

int Teuchos_ParameterList_get_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getParameterList(selfID)->get<int>(std::string(name));
}

double Teuchos_ParameterList_get_const_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getConstParameterList(selfID)->get<double>(std::string(
        name));
}

int Teuchos_ParameterList_get_const_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getConstParameterList(selfID)->get<int>(std::string(
        name));
}

double * Teuchos_ParameterList_getPtr_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getParameterList(selfID)->getPtr<double>(std::string(
        name));
}

int * Teuchos_ParameterList_getPtr_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getParameterList(selfID)->getPtr<int>(std::string(name));
}

const double * Teuchos_ParameterList_getPtr_const_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getConstParameterList(selfID)->getPtr<double>(std::string(
        name));
}

const int * Teuchos_ParameterList_getPtr_const_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::getConstParameterList(selfID)->getPtr<int>(std::string(
        name));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::storeParameterEntry(&( CTeuchos::getParameterList(
        selfID)->getEntry(std::string(name)) ));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::storeConstParameterEntry(
        &( CTeuchos::getConstParameterList(selfID)->getEntry(std::string(
        name)) ));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::storeParameterEntry(CTeuchos::getParameterList(
        selfID)->getEntryPtr(std::string(name)));
}

CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::storeConstParameterEntry(CTeuchos::getConstParameterList(
        selfID)->getEntryPtr(std::string(name)));
}

boolean Teuchos_ParameterList_remove ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  boolean throwIfNotExists )
{
    return ((CTeuchos::getParameterList(selfID)->remove(std::string(name), ((
        throwIfNotExists) != FALSE ? true : false))) ? TRUE : FALSE);
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  boolean mustAlreadyExist, const char docString[] )
{
    return CTeuchos::storeParameterList(&( CTeuchos::getParameterList(
        selfID)->sublist(std::string(name), ((mustAlreadyExist) != 
        FALSE ? true : false), std::string(docString)) ));
    // RAB: This is a segfault waiting to happen!  If the parent PL is
    // destored, this will become a dangling reference and any access will
    // result in undefined behavior!
}

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist_existing ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return CTeuchos::storeConstParameterList(&( CTeuchos::getConstParameterList(
        selfID)->sublist(std::string(name)) ));
    // RAB: This is a segfault waiting to happen (see above)!
}

const char * Teuchos_ParameterList_name_it ( 
  CT_Teuchos_ParameterList_ID_t selfID )
{
    return CTeuchos::getConstParameterList(selfID)->name().c_str();
}

boolean Teuchos_ParameterList_isParameter ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return ((CTeuchos::getConstParameterList(selfID)->isParameter(std::string(
        name))) ? TRUE : FALSE);
}

boolean Teuchos_ParameterList_isSublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return ((CTeuchos::getConstParameterList(selfID)->isSublist(std::string(
        name))) ? TRUE : FALSE);
}

boolean Teuchos_ParameterList_isType_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return ((CTeuchos::getConstParameterList(selfID)->isType<double>(
        std::string(name))) ? TRUE : FALSE);
}

boolean Teuchos_ParameterList_isType_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] )
{
    return ((CTeuchos::getConstParameterList(selfID)->isType<int>(std::string(
        name))) ? TRUE : FALSE);
}

boolean Teuchos_ParameterList_isType_type_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double * ptr )
{
    return ((CTeuchos::getConstParameterList(selfID)->isType<double>(
        std::string(name), ptr)) ? TRUE : FALSE);
}

boolean Teuchos_ParameterList_isType_type_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int * ptr )
{
    return ((CTeuchos::getConstParameterList(selfID)->isType<int>(std::string(
        name), ptr)) ? TRUE : FALSE);
}

const char * Teuchos_ParameterList_currentParametersString ( 
  CT_Teuchos_ParameterList_ID_t selfID )
{
    return CTeuchos::getConstParameterList(
        selfID)->currentParametersString().c_str();
}

void Teuchos_ParameterList_validateParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth, 
  const CT_EValidateUsed_E_t validateUsed, 
  const CT_EValidateDefaults_E_t validateDefaults )
{
    const Teuchos::RCP<const Teuchos::ParameterList> validParamList = 
        CTeuchos::getConstParameterList(validParamListID);
    CTeuchos::getConstParameterList(selfID)->validateParameters(
        *validParamList, depth, (Teuchos::EValidateUsed) validateUsed, 
        (Teuchos::EValidateDefaults) validateDefaults);
}

void Teuchos_ParameterList_validateParametersAndSetDefaults ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth )
{
    const Teuchos::RCP<const Teuchos::ParameterList> validParamList = 
        CTeuchos::getConstParameterList(validParamListID);
    CTeuchos::getParameterList(selfID)->validateParametersAndSetDefaults(
        *validParamList, depth);
}

void Teuchos_ParameterList_Assign ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID )
{
    Teuchos::ParameterList& self = *( CTeuchos::getParameterList(selfID) );

    const Teuchos::RCP<const Teuchos::ParameterList> source = 
        CTeuchos::getConstParameterList(sourceID);
    self = *source;
}


} // extern "C"


//
// Definitions from CTeuchos_ParameterList_Cpp.hpp
//


/* get Teuchos::ParameterList from non-const table using CT_Teuchos_ParameterList_ID */
const Teuchos::RCP<Teuchos::ParameterList>
CTeuchos::getParameterList( CT_Teuchos_ParameterList_ID_t id )
{
    if (tableOfParameterLists().isType(id.table))
        return tableOfParameterLists().get<Teuchos::ParameterList>(
        CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Teuchos::ParameterList>(
        CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id));
}

/* get Teuchos::ParameterList from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Teuchos::ParameterList>
CTeuchos::getParameterList( CTrilinos_Universal_ID_t id )
{
    if (tableOfParameterLists().isType(id.table))
        return tableOfParameterLists().get<Teuchos::ParameterList>(id);
    else
        return CTrilinos::TableRepos::get<Teuchos::ParameterList>(id);
}

/* get const Teuchos::ParameterList from either the const or non-const table
 * using CT_Teuchos_ParameterList_ID */
const Teuchos::RCP<const Teuchos::ParameterList>
CTeuchos::getConstParameterList( CT_Teuchos_ParameterList_ID_t id )
{
    if (tableOfParameterLists().isType(id.table))
        return tableOfParameterLists().getConst<Teuchos::ParameterList>(
        CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Teuchos::ParameterList>(
        CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id));
}

/* get const Teuchos::ParameterList from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Teuchos::ParameterList>
CTeuchos::getConstParameterList( CTrilinos_Universal_ID_t id )
{
    if (tableOfParameterLists().isType(id.table))
        return tableOfParameterLists().getConst<Teuchos::ParameterList>(id);
    else
        return CTrilinos::TableRepos::getConst<Teuchos::ParameterList>(id);
}

/* store Teuchos::ParameterList (owned) in non-const table */
CT_Teuchos_ParameterList_ID_t
CTeuchos::storeNewParameterList( Teuchos::ParameterList *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
        tableOfParameterLists().store<Teuchos::ParameterList>(pobj, true));
}

/* store Teuchos::ParameterList in non-const table */
CT_Teuchos_ParameterList_ID_t
CTeuchos::storeParameterList( Teuchos::ParameterList *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
        tableOfParameterLists().store<Teuchos::ParameterList>(pobj, false));
}

/* store const Teuchos::ParameterList in const table */
CT_Teuchos_ParameterList_ID_t
CTeuchos::storeConstParameterList( const Teuchos::ParameterList *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
        tableOfParameterLists().store<Teuchos::ParameterList>(pobj, false));
}

/* remove Teuchos::ParameterList from table using CT_Teuchos_ParameterList_ID */
void
CTeuchos::removeParameterList( CT_Teuchos_ParameterList_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(*id);
    if (tableOfParameterLists().isType(aid.table))
        tableOfParameterLists().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(aid);
}

/* remove Teuchos::ParameterList from table using CTrilinos_Universal_ID_t */
void
CTeuchos::removeParameterList( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfParameterLists().isType(aid->table))
        tableOfParameterLists().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Teuchos::ParameterList table */
void
CTeuchos::purgeParameterList(  )
{
    tableOfParameterLists().purge();
}

/* store Teuchos::ParameterList in non-const table */
CTrilinos_Universal_ID_t
CTeuchos::aliasParameterList( const Teuchos::RCP< Teuchos::ParameterList > & robj )
{
    return tableOfParameterLists().alias(robj);
}

/* store const Teuchos::ParameterList in const table */
CTrilinos_Universal_ID_t
CTeuchos::aliasConstParameterList( const Teuchos::RCP< const Teuchos::ParameterList > & robj )
{
    return tableOfParameterLists().alias(robj);
}



