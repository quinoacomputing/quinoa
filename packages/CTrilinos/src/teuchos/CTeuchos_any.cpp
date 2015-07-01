
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
#include "CTeuchos_any.h"
#include "CTeuchos_any_Cpp.hpp"
#include "Teuchos_any.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Teuchos::any */
Table<Teuchos::any>& tableOfanys()
{
    static Table<Teuchos::any> loc_tableOfanys(CT_Teuchos_any_ID);
    return loc_tableOfanys;
}


} // namespace


//
// Definitions from CTeuchos_any.h
//


extern "C" {


CT_Teuchos_any_ID_t Teuchos_any_Create (  )
{
    return CTeuchos::storeNewany(new Teuchos::any());
}

CT_Teuchos_any_ID_t Teuchos_any_Create_double ( double value )
{
    return CTeuchos::storeNewany(new Teuchos::any(value));
}

CT_Teuchos_any_ID_t Teuchos_any_Create_int ( int value )
{
    return CTeuchos::storeNewany(new Teuchos::any(value));
}

CT_Teuchos_any_ID_t Teuchos_any_Duplicate ( 
  CT_Teuchos_any_ID_t otherID )
{
    const Teuchos::RCP<const Teuchos::any> other = CTeuchos::getConstany(
        otherID);
    return CTeuchos::storeNewany(new Teuchos::any(*other));
}

void Teuchos_any_Destroy ( CT_Teuchos_any_ID_t * selfID )
{
    CTeuchos::removeany(selfID);
}

CT_Teuchos_any_ID_t Teuchos_any_swap ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID )
{
    const Teuchos::RCP<Teuchos::any> rhs = CTeuchos::getany(rhsID);
    return CTeuchos::storeany(&( CTeuchos::getany(selfID)->swap(*rhs) ));
}

boolean Teuchos_any_empty ( CT_Teuchos_any_ID_t selfID )
{
    return ((CTeuchos::getConstany(selfID)->empty()) ? TRUE : FALSE);
}

const char * Teuchos_any_typeName ( CT_Teuchos_any_ID_t selfID )
{
    return CTeuchos::getConstany(selfID)->typeName().c_str();
}

boolean Teuchos_any_same ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t otherID )
{
    const Teuchos::RCP<const Teuchos::any> other = CTeuchos::getConstany(
        otherID);
    return ((CTeuchos::getConstany(selfID)->same(*other)) ? TRUE : FALSE);
}

void Teuchos_any_Assign ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID )
{
    Teuchos::any& self = *( CTeuchos::getany(selfID) );

    const Teuchos::RCP<const Teuchos::any> rhs = CTeuchos::getConstany(rhsID);
    self = *rhs;
}


} // extern "C"


//
// Definitions from CTeuchos_any_Cpp.hpp
//


/* get Teuchos::any from non-const table using CT_Teuchos_any_ID */
const Teuchos::RCP<Teuchos::any>
CTeuchos::getany( CT_Teuchos_any_ID_t id )
{
    return tableOfanys().get<Teuchos::any>(
        CTrilinos::abstractType<CT_Teuchos_any_ID_t>(id));
}

/* get Teuchos::any from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Teuchos::any>
CTeuchos::getany( CTrilinos_Universal_ID_t id )
{
    return tableOfanys().get<Teuchos::any>(id);
}

/* get const Teuchos::any from either the const or non-const table
 * using CT_Teuchos_any_ID */
const Teuchos::RCP<const Teuchos::any>
CTeuchos::getConstany( CT_Teuchos_any_ID_t id )
{
    return tableOfanys().getConst<Teuchos::any>(
        CTrilinos::abstractType<CT_Teuchos_any_ID_t>(id));
}

/* get const Teuchos::any from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Teuchos::any>
CTeuchos::getConstany( CTrilinos_Universal_ID_t id )
{
    return tableOfanys().getConst<Teuchos::any>(id);
}

/* store Teuchos::any (owned) in non-const table */
CT_Teuchos_any_ID_t
CTeuchos::storeNewany( Teuchos::any *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_any_ID_t>(
        tableOfanys().store<Teuchos::any>(pobj, true));
}

/* store Teuchos::any in non-const table */
CT_Teuchos_any_ID_t
CTeuchos::storeany( Teuchos::any *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_any_ID_t>(
        tableOfanys().store<Teuchos::any>(pobj, false));
}

/* store const Teuchos::any in const table */
CT_Teuchos_any_ID_t
CTeuchos::storeConstany( const Teuchos::any *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_any_ID_t>(
        tableOfanys().store<Teuchos::any>(pobj, false));
}

/* remove Teuchos::any from table using CT_Teuchos_any_ID */
void
CTeuchos::removeany( CT_Teuchos_any_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Teuchos_any_ID_t>(*id);
    tableOfanys().remove(&aid);
    *id = CTrilinos::concreteType<CT_Teuchos_any_ID_t>(aid);
}

/* remove Teuchos::any from table using CTrilinos_Universal_ID_t */
void
CTeuchos::removeany( CTrilinos_Universal_ID_t *aid )
{
    tableOfanys().remove(aid);
}

/* purge Teuchos::any table */
void
CTeuchos::purgeany(  )
{
    tableOfanys().purge();
}



