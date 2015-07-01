#ifndef CEPETRA_TIME_CPP_HPP
#define CEPETRA_TIME_CPP_HPP

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
#include "Epetra_Time.h"


namespace CEpetra {


using Teuchos::RCP;


/*! get Epetra_Time from non-const table using CT_Epetra_Time_ID */
const RCP<Epetra_Time>
getTime( CT_Epetra_Time_ID_t id );

/*! get Epetra_Time from non-const table using CTrilinos_Universal_ID_t */
const RCP<Epetra_Time>
getTime( CTrilinos_Universal_ID_t id );

/*! get const Epetra_Time from either the const or non-const table
 * using CT_Epetra_Time_ID */
const RCP<const Epetra_Time>
getConstTime( CT_Epetra_Time_ID_t id );

/*! get const Epetra_Time from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Epetra_Time>
getConstTime( CTrilinos_Universal_ID_t id );

/*! store Epetra_Time (owned) in non-const table */
CT_Epetra_Time_ID_t
storeNewTime( Epetra_Time *pobj );

/*! store Epetra_Time in non-const table */
CT_Epetra_Time_ID_t
storeTime( Epetra_Time *pobj );

/*! store const Epetra_Time in const table */
CT_Epetra_Time_ID_t
storeConstTime( const Epetra_Time *pobj );

/* remove Epetra_Time from table using CT_Epetra_Time_ID */
void
removeTime( CT_Epetra_Time_ID_t *id );

/* remove Epetra_Time from table using CTrilinos_Universal_ID_t */
void
removeTime( CTrilinos_Universal_ID_t *aid );

/* purge Epetra_Time table */
void
purgeTime(  );

/* store Epetra_Time in non-const table */
CTrilinos_Universal_ID_t
aliasTime( const Teuchos::RCP< Epetra_Time > & robj );

/* store const Epetra_Time in const table */
CTrilinos_Universal_ID_t
aliasConstTime( const Teuchos::RCP< const Epetra_Time > & robj );

} // namespace CEpetra


#endif // CEPETRA_TIME_CPP_HPP


