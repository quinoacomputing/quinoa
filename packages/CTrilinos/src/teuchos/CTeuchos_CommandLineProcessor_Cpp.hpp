#ifndef CTEUCHOS_COMMANDLINEPROCESSOR_CPP_HPP
#define CTEUCHOS_COMMANDLINEPROCESSOR_CPP_HPP

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
#include "Teuchos_CommandLineProcessor.hpp"


namespace CTeuchos {


using Teuchos::RCP;


/*! get Teuchos::CommandLineProcessor from non-const table using CT_Teuchos_CommandLineProcessor_ID */
const RCP<Teuchos::CommandLineProcessor>
getCommandLineProcessor( CT_Teuchos_CommandLineProcessor_ID_t id );

/*! get Teuchos::CommandLineProcessor from non-const table using CTrilinos_Universal_ID_t */
const RCP<Teuchos::CommandLineProcessor>
getCommandLineProcessor( CTrilinos_Universal_ID_t id );

/*! get const Teuchos::CommandLineProcessor from either the const or non-const table
 * using CT_Teuchos_CommandLineProcessor_ID */
const RCP<const Teuchos::CommandLineProcessor>
getConstCommandLineProcessor( CT_Teuchos_CommandLineProcessor_ID_t id );

/*! get const Teuchos::CommandLineProcessor from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Teuchos::CommandLineProcessor>
getConstCommandLineProcessor( CTrilinos_Universal_ID_t id );

/*! store Teuchos::CommandLineProcessor (owned) in non-const table */
CT_Teuchos_CommandLineProcessor_ID_t
storeNewCommandLineProcessor( Teuchos::CommandLineProcessor *pobj );

/*! store Teuchos::CommandLineProcessor in non-const table */
CT_Teuchos_CommandLineProcessor_ID_t
storeCommandLineProcessor( Teuchos::CommandLineProcessor *pobj );

/*! store const Teuchos::CommandLineProcessor in const table */
CT_Teuchos_CommandLineProcessor_ID_t
storeConstCommandLineProcessor( const Teuchos::CommandLineProcessor *pobj );

/* remove Teuchos::CommandLineProcessor from table using CT_Teuchos_CommandLineProcessor_ID */
void
removeCommandLineProcessor( CT_Teuchos_CommandLineProcessor_ID_t *id );

/* remove Teuchos::CommandLineProcessor from table using CTrilinos_Universal_ID_t */
void
removeCommandLineProcessor( CTrilinos_Universal_ID_t *aid );

/* purge Teuchos::CommandLineProcessor table */
void
purgeCommandLineProcessor(  );

} // namespace CTeuchos


#endif // CTEUCHOS_COMMANDLINEPROCESSOR_CPP_HPP


