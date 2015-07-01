
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
#include "CIfpack.h"
#include "CIfpack_Cpp.hpp"
#include "Ifpack.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"
#include "CTrilinos_difficult.hpp"
#include "CIfpack_Preconditioner_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Ifpack */
Table<Ifpack>& tableOfIfpacks()
{
    static Table<Ifpack> loc_tableOfIfpacks(CT_Ifpack_ID);
    return loc_tableOfIfpacks;
}


} // namespace


//
// Definitions from CIfpack.h
//


extern "C" {


CT_Ifpack_ID_t Ifpack_Create (  )
{
    return CIfpack::storeNewIfpack(new Ifpack());
}

void Ifpack_Destroy ( CT_Ifpack_ID_t * selfID )
{
    CIfpack::removeIfpack(selfID);
}

CT_Ifpack_Preconditioner_ID_t Ifpack_CreatePreconditioner_UsingName ( 
  CT_Ifpack_ID_t selfID, const char PrecType[], 
  CT_Epetra_RowMatrix_ID_t MatrixID, const int overlap )
{
    const Teuchos::RCP<Epetra_RowMatrix> Matrix = CEpetra::getRowMatrix(
        MatrixID);
    return CIfpack::storePreconditioner(CIfpack::getIfpack(selfID)->Create(
        std::string(PrecType), Matrix.getRawPtr(), overlap));
}

int Ifpack_SetParameters ( 
  CT_Ifpack_ID_t selfID, int argc, char * argv[], 
  CT_Teuchos_ParameterList_ID_t ListID, char * PrecType[], 
  int * Overlap )
{
    const Teuchos::RCP<Teuchos::ParameterList> List = 
        CTeuchos::getParameterList(ListID);
    std::string *tmp_PrecType = NULL;
    CTrilinos::pass_string_in(PrecType, tmp_PrecType);
    int ret = CIfpack::getIfpack(selfID)->SetParameters(argc, argv, *List, 
        *tmp_PrecType, *Overlap);
    CTrilinos::pass_string_out(tmp_PrecType, PrecType);
    delete tmp_PrecType;

    return ret;
}

const char * Ifpack_toString ( const CT_EPrecType_E_t precType )
{
    return Ifpack::toString(CTrilinos::convert_to_difficult_enum(precType));
}

CT_Ifpack_Preconditioner_ID_t Ifpack_CreatePreconditioner_UsingType ( 
  CT_EPrecType_E_t PrecType, CT_Epetra_RowMatrix_ID_t MatrixID, 
  const int overlap )
{
    const Teuchos::RCP<Epetra_RowMatrix> Matrix = CEpetra::getRowMatrix(
        MatrixID);
    return CIfpack::storePreconditioner(Ifpack::Create(
        CTrilinos::convert_to_difficult_enum(PrecType), Matrix.getRawPtr(), 
        overlap));
}


} // extern "C"


//
// Definitions from CIfpack_Cpp.hpp
//


/* get Ifpack from non-const table using CT_Ifpack_ID */
const Teuchos::RCP<Ifpack>
CIfpack::getIfpack( CT_Ifpack_ID_t id )
{
    return tableOfIfpacks().get<Ifpack>(
        CTrilinos::abstractType<CT_Ifpack_ID_t>(id));
}

/* get Ifpack from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Ifpack>
CIfpack::getIfpack( CTrilinos_Universal_ID_t id )
{
    return tableOfIfpacks().get<Ifpack>(id);
}

/* get const Ifpack from either the const or non-const table
 * using CT_Ifpack_ID */
const Teuchos::RCP<const Ifpack>
CIfpack::getConstIfpack( CT_Ifpack_ID_t id )
{
    return tableOfIfpacks().getConst<Ifpack>(
        CTrilinos::abstractType<CT_Ifpack_ID_t>(id));
}

/* get const Ifpack from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Ifpack>
CIfpack::getConstIfpack( CTrilinos_Universal_ID_t id )
{
    return tableOfIfpacks().getConst<Ifpack>(id);
}

/* store Ifpack (owned) in non-const table */
CT_Ifpack_ID_t
CIfpack::storeNewIfpack( Ifpack *pobj )
{
    return CTrilinos::concreteType<CT_Ifpack_ID_t>(
        tableOfIfpacks().store<Ifpack>(pobj, true));
}

/* store Ifpack in non-const table */
CT_Ifpack_ID_t
CIfpack::storeIfpack( Ifpack *pobj )
{
    return CTrilinos::concreteType<CT_Ifpack_ID_t>(
        tableOfIfpacks().store<Ifpack>(pobj, false));
}

/* store const Ifpack in const table */
CT_Ifpack_ID_t
CIfpack::storeConstIfpack( const Ifpack *pobj )
{
    return CTrilinos::concreteType<CT_Ifpack_ID_t>(
        tableOfIfpacks().store<Ifpack>(pobj, false));
}

/* remove Ifpack from table using CT_Ifpack_ID */
void
CIfpack::removeIfpack( CT_Ifpack_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Ifpack_ID_t>(*id);
    tableOfIfpacks().remove(&aid);
    *id = CTrilinos::concreteType<CT_Ifpack_ID_t>(aid);
}

/* remove Ifpack from table using CTrilinos_Universal_ID_t */
void
CIfpack::removeIfpack( CTrilinos_Universal_ID_t *aid )
{
    tableOfIfpacks().remove(aid);
}

/* purge Ifpack table */
void
CIfpack::purgeIfpack(  )
{
    tableOfIfpacks().purge();
}



#endif /* HAVE_CTRILINOS_IFPACK */


