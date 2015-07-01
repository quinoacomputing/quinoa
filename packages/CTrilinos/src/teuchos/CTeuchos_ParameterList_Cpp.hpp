#ifndef CTEUCHOS_PARAMETERLIST_CPP_HPP
#define CTEUCHOS_PARAMETERLIST_CPP_HPP

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
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"


namespace CTeuchos {


using Teuchos::RCP;


/*! get Teuchos::ParameterList from non-const table using CT_Teuchos_ParameterList_ID */
const RCP<Teuchos::ParameterList>
getParameterList( CT_Teuchos_ParameterList_ID_t id );

/*! get Teuchos::ParameterList from non-const table using CTrilinos_Universal_ID_t */
const RCP<Teuchos::ParameterList>
getParameterList( CTrilinos_Universal_ID_t id );

/*! get const Teuchos::ParameterList from either the const or non-const table
 * using CT_Teuchos_ParameterList_ID */
const RCP<const Teuchos::ParameterList>
getConstParameterList( CT_Teuchos_ParameterList_ID_t id );

/*! get const Teuchos::ParameterList from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Teuchos::ParameterList>
getConstParameterList( CTrilinos_Universal_ID_t id );

/*! store Teuchos::ParameterList (owned) in non-const table */
CT_Teuchos_ParameterList_ID_t
storeNewParameterList( Teuchos::ParameterList *pobj );

/*! store Teuchos::ParameterList in non-const table */
CT_Teuchos_ParameterList_ID_t
storeParameterList( Teuchos::ParameterList *pobj );

/*! store const Teuchos::ParameterList in const table */
CT_Teuchos_ParameterList_ID_t
storeConstParameterList( const Teuchos::ParameterList *pobj );

/* remove Teuchos::ParameterList from table using CT_Teuchos_ParameterList_ID */
void
removeParameterList( CT_Teuchos_ParameterList_ID_t *id );

/* remove Teuchos::ParameterList from table using CTrilinos_Universal_ID_t */
void
removeParameterList( CTrilinos_Universal_ID_t *aid );

/* purge Teuchos::ParameterList table */
void
purgeParameterList(  );

/* store Teuchos::ParameterList in non-const table */
CTrilinos_Universal_ID_t
aliasParameterList( const Teuchos::RCP< Teuchos::ParameterList > & robj );

/* store const Teuchos::ParameterList in const table */
CTrilinos_Universal_ID_t
aliasConstParameterList( const Teuchos::RCP< const Teuchos::ParameterList > & robj );

} // namespace CTeuchos


#endif // CTEUCHOS_PARAMETERLIST_CPP_HPP


