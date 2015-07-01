#ifndef CEPETRA_MAP_CPP_HPP
#define CEPETRA_MAP_CPP_HPP

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
#include "Epetra_Map.h"


namespace CEpetra {


using Teuchos::RCP;


/*! get Epetra_Map from non-const table using CT_Epetra_Map_ID */
const RCP<Epetra_Map>
getMap( CT_Epetra_Map_ID_t id );

/*! get Epetra_Map from non-const table using CTrilinos_Universal_ID_t */
const RCP<Epetra_Map>
getMap( CTrilinos_Universal_ID_t id );

/*! get const Epetra_Map from either the const or non-const table
 * using CT_Epetra_Map_ID */
const RCP<const Epetra_Map>
getConstMap( CT_Epetra_Map_ID_t id );

/*! get const Epetra_Map from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Epetra_Map>
getConstMap( CTrilinos_Universal_ID_t id );

/*! store Epetra_Map (owned) in non-const table */
CT_Epetra_Map_ID_t
storeNewMap( Epetra_Map *pobj );

/*! store Epetra_Map in non-const table */
CT_Epetra_Map_ID_t
storeMap( Epetra_Map *pobj );

/*! store const Epetra_Map in const table */
CT_Epetra_Map_ID_t
storeConstMap( const Epetra_Map *pobj );

/* remove Epetra_Map from table using CT_Epetra_Map_ID */
void
removeMap( CT_Epetra_Map_ID_t *id );

/* remove Epetra_Map from table using CTrilinos_Universal_ID_t */
void
removeMap( CTrilinos_Universal_ID_t *aid );

/* purge Epetra_Map table */
void
purgeMap(  );

/* store Epetra_Map in non-const table */
CTrilinos_Universal_ID_t
aliasMap( const Teuchos::RCP< Epetra_Map > & robj );

/* store const Epetra_Map in const table */
CTrilinos_Universal_ID_t
aliasConstMap( const Teuchos::RCP< const Epetra_Map > & robj );

} // namespace CEpetra


#endif // CEPETRA_MAP_CPP_HPP


