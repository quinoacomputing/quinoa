#ifndef CTEUCHOS_PARAMETERENTRY_CPP_HPP
#define CTEUCHOS_PARAMETERENTRY_CPP_HPP

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
#include "Teuchos_ParameterEntry.hpp"


namespace CTeuchos {


using Teuchos::RCP;


/*! get Teuchos::ParameterEntry from non-const table using CT_Teuchos_ParameterEntry_ID */
const RCP<Teuchos::ParameterEntry>
getParameterEntry( CT_Teuchos_ParameterEntry_ID_t id );

/*! get Teuchos::ParameterEntry from non-const table using CTrilinos_Universal_ID_t */
const RCP<Teuchos::ParameterEntry>
getParameterEntry( CTrilinos_Universal_ID_t id );

/*! get const Teuchos::ParameterEntry from either the const or non-const table
 * using CT_Teuchos_ParameterEntry_ID */
const RCP<const Teuchos::ParameterEntry>
getConstParameterEntry( CT_Teuchos_ParameterEntry_ID_t id );

/*! get const Teuchos::ParameterEntry from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Teuchos::ParameterEntry>
getConstParameterEntry( CTrilinos_Universal_ID_t id );

/*! store Teuchos::ParameterEntry (owned) in non-const table */
CT_Teuchos_ParameterEntry_ID_t
storeNewParameterEntry( Teuchos::ParameterEntry *pobj );

/*! store Teuchos::ParameterEntry in non-const table */
CT_Teuchos_ParameterEntry_ID_t
storeParameterEntry( Teuchos::ParameterEntry *pobj );

/*! store const Teuchos::ParameterEntry in const table */
CT_Teuchos_ParameterEntry_ID_t
storeConstParameterEntry( const Teuchos::ParameterEntry *pobj );

/* remove Teuchos::ParameterEntry from table using CT_Teuchos_ParameterEntry_ID */
void
removeParameterEntry( CT_Teuchos_ParameterEntry_ID_t *id );

/* remove Teuchos::ParameterEntry from table using CTrilinos_Universal_ID_t */
void
removeParameterEntry( CTrilinos_Universal_ID_t *aid );

/* purge Teuchos::ParameterEntry table */
void
purgeParameterEntry(  );

} // namespace CTeuchos


#endif // CTEUCHOS_PARAMETERENTRY_CPP_HPP


