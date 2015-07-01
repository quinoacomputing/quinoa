#ifndef CEPETRA_COMPOBJECT_CPP_HPP
#define CEPETRA_COMPOBJECT_CPP_HPP

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
#include "Epetra_CompObject.h"


namespace CEpetra {


using Teuchos::RCP;


/*! get Epetra_CompObject from non-const table using CT_Epetra_CompObject_ID */
const RCP<Epetra_CompObject>
getCompObject( CT_Epetra_CompObject_ID_t id );

/*! get Epetra_CompObject from non-const table using CTrilinos_Universal_ID_t */
const RCP<Epetra_CompObject>
getCompObject( CTrilinos_Universal_ID_t id );

/*! get const Epetra_CompObject from either the const or non-const table
 * using CT_Epetra_CompObject_ID */
const RCP<const Epetra_CompObject>
getConstCompObject( CT_Epetra_CompObject_ID_t id );

/*! get const Epetra_CompObject from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Epetra_CompObject>
getConstCompObject( CTrilinos_Universal_ID_t id );

/*! store Epetra_CompObject (owned) in non-const table */
CT_Epetra_CompObject_ID_t
storeNewCompObject( Epetra_CompObject *pobj );

/*! store Epetra_CompObject in non-const table */
CT_Epetra_CompObject_ID_t
storeCompObject( Epetra_CompObject *pobj );

/*! store const Epetra_CompObject in const table */
CT_Epetra_CompObject_ID_t
storeConstCompObject( const Epetra_CompObject *pobj );

/* remove Epetra_CompObject from table using CT_Epetra_CompObject_ID */
void
removeCompObject( CT_Epetra_CompObject_ID_t *id );

/* remove Epetra_CompObject from table using CTrilinos_Universal_ID_t */
void
removeCompObject( CTrilinos_Universal_ID_t *aid );

/* purge Epetra_CompObject table */
void
purgeCompObject(  );

/* store Epetra_CompObject in non-const table */
CTrilinos_Universal_ID_t
aliasCompObject( const Teuchos::RCP< Epetra_CompObject > & robj );

/* store const Epetra_CompObject in const table */
CTrilinos_Universal_ID_t
aliasConstCompObject( const Teuchos::RCP< const Epetra_CompObject > & robj );

} // namespace CEpetra


#endif // CEPETRA_COMPOBJECT_CPP_HPP


