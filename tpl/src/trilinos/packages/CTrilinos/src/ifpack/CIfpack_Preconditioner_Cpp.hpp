#ifndef CIFPACK_PRECONDITIONER_CPP_HPP
#define CIFPACK_PRECONDITIONER_CPP_HPP

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


#ifdef HAVE_CTRILINOS_IFPACK


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Ifpack_Preconditioner.h"


namespace CIfpack {


using Teuchos::RCP;


/*! get Ifpack_Preconditioner from non-const table using CT_Ifpack_Preconditioner_ID */
const RCP<Ifpack_Preconditioner>
getPreconditioner( CT_Ifpack_Preconditioner_ID_t id );

/*! get Ifpack_Preconditioner from non-const table using CTrilinos_Universal_ID_t */
const RCP<Ifpack_Preconditioner>
getPreconditioner( CTrilinos_Universal_ID_t id );

/*! get const Ifpack_Preconditioner from either the const or non-const table
 * using CT_Ifpack_Preconditioner_ID */
const RCP<const Ifpack_Preconditioner>
getConstPreconditioner( CT_Ifpack_Preconditioner_ID_t id );

/*! get const Ifpack_Preconditioner from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Ifpack_Preconditioner>
getConstPreconditioner( CTrilinos_Universal_ID_t id );

/*! store Ifpack_Preconditioner (owned) in non-const table */
CT_Ifpack_Preconditioner_ID_t
storeNewPreconditioner( Ifpack_Preconditioner *pobj );

/*! store Ifpack_Preconditioner in non-const table */
CT_Ifpack_Preconditioner_ID_t
storePreconditioner( Ifpack_Preconditioner *pobj );

/*! store const Ifpack_Preconditioner in const table */
CT_Ifpack_Preconditioner_ID_t
storeConstPreconditioner( const Ifpack_Preconditioner *pobj );

/* remove Ifpack_Preconditioner from table using CT_Ifpack_Preconditioner_ID */
void
removePreconditioner( CT_Ifpack_Preconditioner_ID_t *id );

/* remove Ifpack_Preconditioner from table using CTrilinos_Universal_ID_t */
void
removePreconditioner( CTrilinos_Universal_ID_t *aid );

/* purge Ifpack_Preconditioner table */
void
purgePreconditioner(  );

/* store Ifpack_Preconditioner in non-const table */
CTrilinos_Universal_ID_t
aliasPreconditioner( const Teuchos::RCP< Ifpack_Preconditioner > & robj );

/* store const Ifpack_Preconditioner in const table */
CTrilinos_Universal_ID_t
aliasConstPreconditioner( const Teuchos::RCP< const Ifpack_Preconditioner > & robj );

} // namespace CIfpack


#endif /* HAVE_CTRILINOS_IFPACK */
#endif // CIFPACK_PRECONDITIONER_CPP_HPP


