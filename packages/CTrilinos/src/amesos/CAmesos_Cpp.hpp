#ifndef CAMESOS_CPP_HPP
#define CAMESOS_CPP_HPP

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


#ifdef HAVE_CTRILINOS_AMESOS


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Amesos.h"


namespace CAmesos {


using Teuchos::RCP;


/*! get Amesos from non-const table using CT_Amesos_ID */
const RCP<Amesos>
getAmesos( CT_Amesos_ID_t id );

/*! get Amesos from non-const table using CTrilinos_Universal_ID_t */
const RCP<Amesos>
getAmesos( CTrilinos_Universal_ID_t id );

/*! get const Amesos from either the const or non-const table
 * using CT_Amesos_ID */
const RCP<const Amesos>
getConstAmesos( CT_Amesos_ID_t id );

/*! get const Amesos from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Amesos>
getConstAmesos( CTrilinos_Universal_ID_t id );

/*! store Amesos (owned) in non-const table */
CT_Amesos_ID_t
storeNewAmesos( Amesos *pobj );

/*! store Amesos in non-const table */
CT_Amesos_ID_t
storeAmesos( Amesos *pobj );

/*! store const Amesos in const table */
CT_Amesos_ID_t
storeConstAmesos( const Amesos *pobj );

/* remove Amesos from table using CT_Amesos_ID */
void
removeAmesos( CT_Amesos_ID_t *id );

/* remove Amesos from table using CTrilinos_Universal_ID_t */
void
removeAmesos( CTrilinos_Universal_ID_t *aid );

/* purge Amesos table */
void
purgeAmesos(  );

} // namespace CAmesos


#endif /* HAVE_CTRILINOS_AMESOS */
#endif // CAMESOS_CPP_HPP


