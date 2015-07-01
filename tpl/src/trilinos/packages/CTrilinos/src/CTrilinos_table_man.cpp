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

#include "CTrilinos_table_man.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_TableRepos.hpp"


extern "C" {


/*! Copies the RCP from one table into a second table. The new ID
 *  will be returned from the function. Both the old and the new
 *  IDs will need to be removed from the tables in order to destroy
 *  the object. */
CTrilinos_Universal_ID_t CT_Alias(CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t new_table)
{
    return CTrilinos::TableRepos::alias(aid, new_table, true);
}

/*! Removes the RCP from one table and puts it in another. *aid will
 *  hold the new struct value afterward. Only the new RCP will need
 *  to be removed in order to destroy the object. */
void CT_Migrate(CTrilinos_Universal_ID_t *aid, CTrilinos_Table_ID_t new_table)
{
    CTrilinos_Universal_ID_t newid = CTrilinos::TableRepos::alias(*aid, new_table, false);
    *aid = newid;
}

/*! Checks to see if the underlying object referenced by a table
 *  entry is dynamic_cast'able to a given type (can be used to
 *  distinguish, e.g., an Epetra_SerialComm from an Epetra_MpiComm). */
boolean CT_TypeCheck(CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t type)
{
    bool ret = CTrilinos::TableRepos::typeCheck(aid, type);
    return (ret ? TRUE : FALSE);
}


} // extern "C"
