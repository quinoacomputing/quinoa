#ifndef CAZTECOO_STATUSTESTRESNORM_CPP_HPP
#define CAZTECOO_STATUSTESTRESNORM_CPP_HPP

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
#include "Teuchos_RCP.hpp"
#include "AztecOO_StatusTestResNorm.h"


namespace CAztecOO {


using Teuchos::RCP;


/*! get AztecOO_StatusTestResNorm from non-const table using CT_AztecOO_StatusTestResNorm_ID */
const RCP<AztecOO_StatusTestResNorm>
getStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t id );

/*! get AztecOO_StatusTestResNorm from non-const table using CTrilinos_Universal_ID_t */
const RCP<AztecOO_StatusTestResNorm>
getStatusTestResNorm( CTrilinos_Universal_ID_t id );

/*! get const AztecOO_StatusTestResNorm from either the const or non-const table
 * using CT_AztecOO_StatusTestResNorm_ID */
const RCP<const AztecOO_StatusTestResNorm>
getConstStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t id );

/*! get const AztecOO_StatusTestResNorm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const AztecOO_StatusTestResNorm>
getConstStatusTestResNorm( CTrilinos_Universal_ID_t id );

/*! store AztecOO_StatusTestResNorm (owned) in non-const table */
CT_AztecOO_StatusTestResNorm_ID_t
storeNewStatusTestResNorm( AztecOO_StatusTestResNorm *pobj );

/*! store AztecOO_StatusTestResNorm in non-const table */
CT_AztecOO_StatusTestResNorm_ID_t
storeStatusTestResNorm( AztecOO_StatusTestResNorm *pobj );

/*! store const AztecOO_StatusTestResNorm in const table */
CT_AztecOO_StatusTestResNorm_ID_t
storeConstStatusTestResNorm( const AztecOO_StatusTestResNorm *pobj );

/* remove AztecOO_StatusTestResNorm from table using CT_AztecOO_StatusTestResNorm_ID */
void
removeStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t *id );

/* remove AztecOO_StatusTestResNorm from table using CTrilinos_Universal_ID_t */
void
removeStatusTestResNorm( CTrilinos_Universal_ID_t *aid );

/* purge AztecOO_StatusTestResNorm table */
void
purgeStatusTestResNorm(  );

/* store AztecOO_StatusTestResNorm in non-const table */
CTrilinos_Universal_ID_t
aliasStatusTestResNorm( const Teuchos::RCP< AztecOO_StatusTestResNorm > & robj );

/* store const AztecOO_StatusTestResNorm in const table */
CTrilinos_Universal_ID_t
aliasConstStatusTestResNorm( const Teuchos::RCP< const AztecOO_StatusTestResNorm > & robj );

} // namespace CAztecOO


#endif /* HAVE_CTRILINOS_AZTECOO */
#endif // CAZTECOO_STATUSTESTRESNORM_CPP_HPP


