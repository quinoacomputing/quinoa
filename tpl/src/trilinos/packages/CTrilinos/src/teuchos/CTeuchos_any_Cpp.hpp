#ifndef CTEUCHOS_ANY_CPP_HPP
#define CTEUCHOS_ANY_CPP_HPP

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
#include "Teuchos_any.hpp"


namespace CTeuchos {


using Teuchos::RCP;


/*! get Teuchos::any from non-const table using CT_Teuchos_any_ID */
const RCP<Teuchos::any>
getany( CT_Teuchos_any_ID_t id );

/*! get Teuchos::any from non-const table using CTrilinos_Universal_ID_t */
const RCP<Teuchos::any>
getany( CTrilinos_Universal_ID_t id );

/*! get const Teuchos::any from either the const or non-const table
 * using CT_Teuchos_any_ID */
const RCP<const Teuchos::any>
getConstany( CT_Teuchos_any_ID_t id );

/*! get const Teuchos::any from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Teuchos::any>
getConstany( CTrilinos_Universal_ID_t id );

/*! store Teuchos::any (owned) in non-const table */
CT_Teuchos_any_ID_t
storeNewany( Teuchos::any *pobj );

/*! store Teuchos::any in non-const table */
CT_Teuchos_any_ID_t
storeany( Teuchos::any *pobj );

/*! store const Teuchos::any in const table */
CT_Teuchos_any_ID_t
storeConstany( const Teuchos::any *pobj );

/* remove Teuchos::any from table using CT_Teuchos_any_ID */
void
removeany( CT_Teuchos_any_ID_t *id );

/* remove Teuchos::any from table using CTrilinos_Universal_ID_t */
void
removeany( CTrilinos_Universal_ID_t *aid );

/* purge Teuchos::any table */
void
purgeany(  );

} // namespace CTeuchos


#endif // CTEUCHOS_ANY_CPP_HPP


