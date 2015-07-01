
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

#ifdef HAVE_CTRILINOS_PLIRIS

#ifdef HAVE_MPI

#include "CTrilinos_enums.h"
#include "CPliris.h"
#include "CPliris_Cpp.hpp"
#include "Pliris.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_SerialDenseVector_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Pliris */
Table<Pliris>& tableOfPliriss()
{
    static Table<Pliris> loc_tableOfPliriss(CT_Pliris_ID);
    return loc_tableOfPliriss;
}


} // namespace


//
// Definitions from CPliris.h
//


extern "C" {


CT_Pliris_ID_t Pliris_Create ( 
  CT_Epetra_Vector_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_Vector> A = CEpetra::getVector(AID);
    const Teuchos::RCP<Epetra_MultiVector> X = CEpetra::getMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> B = CEpetra::getMultiVector(BID);
    return CPliris::storeNewPliris(new Pliris(A.getRawPtr(), X.getRawPtr(), 
        B.getRawPtr()));
}

CT_Pliris_ID_t Pliris_Create_Default (  )
{
    return CPliris::storeNewPliris(new Pliris());
}

void Pliris_Destroy ( CT_Pliris_ID_t * selfID )
{
    CPliris::removePliris(selfID);
}

int Pliris_SetLHS ( 
  CT_Pliris_ID_t selfID, CT_Epetra_MultiVector_ID_t XID )
{
    const Teuchos::RCP<Epetra_MultiVector> X = CEpetra::getMultiVector(XID);
    return CPliris::getPliris(selfID)->SetLHS(X.getRawPtr());
}

int Pliris_SetRHS ( 
  CT_Pliris_ID_t selfID, CT_Epetra_MultiVector_ID_t BID )
{
    const Teuchos::RCP<Epetra_MultiVector> B = CEpetra::getMultiVector(BID);
    return CPliris::getPliris(selfID)->SetRHS(B.getRawPtr());
}

int Pliris_SetMatrix ( 
  CT_Pliris_ID_t selfID, CT_Epetra_Vector_ID_t AID )
{
    const Teuchos::RCP<Epetra_Vector> A = CEpetra::getVector(AID);
    return CPliris::getPliris(selfID)->SetMatrix(A.getRawPtr());
}

int Pliris_SetMatrix_Serial ( 
  CT_Pliris_ID_t selfID, CT_Epetra_SerialDenseVector_ID_t AID )
{
    const Teuchos::RCP<Epetra_SerialDenseVector> A = 
        CEpetra::getSerialDenseVector(AID);
    return CPliris::getPliris(selfID)->SetMatrix(A.getRawPtr());
}

int Pliris_GetDistribution ( 
  CT_Pliris_ID_t selfID, int * nprocs_row, int * number_of_unknowns, 
  int * nrhs, int * my_rows, int * my_cols, int * my_first_row, 
  int * my_first_col, int * my_rhs, int * my_row, int * my_col )
{
    return CPliris::getPliris(selfID)->GetDistribution(nprocs_row, 
        number_of_unknowns, nrhs, my_rows, my_cols, my_first_row, 
        my_first_col, my_rhs, my_row, my_col);
}

int Pliris_FactorSolve ( 
  CT_Pliris_ID_t selfID, CT_Epetra_Vector_ID_t AID, int my_rows, 
  int my_cols, int * matrix_size, int * num_procsr, int * num_rhs, 
  double * secs )
{
    const Teuchos::RCP<Epetra_Vector> A = CEpetra::getVector(AID);
    return CPliris::getPliris(selfID)->FactorSolve(A.getRawPtr(), my_rows, 
        my_cols, matrix_size, num_procsr, num_rhs, secs);
}

int Pliris_FactorSolve_Serial ( 
  CT_Pliris_ID_t selfID, CT_Epetra_SerialDenseVector_ID_t AAID, 
  int my_rows, int my_cols, int * matrix_size, int * num_procsr, 
  int * num_rhs, double * secs )
{
    const Teuchos::RCP<Epetra_SerialDenseVector> AA = 
        CEpetra::getSerialDenseVector(AAID);
    return CPliris::getPliris(selfID)->FactorSolve(AA.getRawPtr(), my_rows, 
        my_cols, matrix_size, num_procsr, num_rhs, secs);
}

int Pliris_Factor ( 
  CT_Pliris_ID_t selfID, CT_Epetra_Vector_ID_t AID, 
  int * matrix_size, int * num_procsr, int * permute, 
  double * secs )
{
    const Teuchos::RCP<Epetra_Vector> A = CEpetra::getVector(AID);
    return CPliris::getPliris(selfID)->Factor(A.getRawPtr(), matrix_size, 
        num_procsr, permute, secs);
}

int Pliris_Solve ( 
  CT_Pliris_ID_t selfID, int * permute, int * num_rhs )
{
    return CPliris::getPliris(selfID)->Solve(permute, num_rhs);
}


} // extern "C"


//
// Definitions from CPliris_Cpp.hpp
//


/* get Pliris from non-const table using CT_Pliris_ID */
const Teuchos::RCP<Pliris>
CPliris::getPliris( CT_Pliris_ID_t id )
{
        return tableOfPliriss().get<Pliris>(
        CTrilinos::abstractType<CT_Pliris_ID_t>(id));
}

/* get Pliris from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Pliris>
CPliris::getPliris( CTrilinos_Universal_ID_t id )
{
        return tableOfPliriss().get<Pliris>(id);
}

/* get const Pliris from either the const or non-const table
 * using CT_Pliris_ID */
const Teuchos::RCP<const Pliris>
CPliris::getConstPliris( CT_Pliris_ID_t id )
{
        return tableOfPliriss().getConst<Pliris>(
        CTrilinos::abstractType<CT_Pliris_ID_t>(id));
}

/* get const Pliris from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Pliris>
CPliris::getConstPliris( CTrilinos_Universal_ID_t id )
{
        return tableOfPliriss().getConst<Pliris>(id);
}

/* store Pliris (owned) in non-const table */
CT_Pliris_ID_t
CPliris::storeNewPliris( Pliris *pobj )
{
    return CTrilinos::concreteType<CT_Pliris_ID_t>(
        tableOfPliriss().store<Pliris>(pobj, true));
}

/* store Pliris in non-const table */
CT_Pliris_ID_t
CPliris::storePliris( Pliris *pobj )
{
    return CTrilinos::concreteType<CT_Pliris_ID_t>(
        tableOfPliriss().store<Pliris>(pobj, false));
}

/* store const Pliris in const table */
CT_Pliris_ID_t
CPliris::storeConstPliris( const Pliris *pobj )
{
    return CTrilinos::concreteType<CT_Pliris_ID_t>(
        tableOfPliriss().store<Pliris>(pobj, false));
}

/* remove Pliris from table using CT_Pliris_ID */
void
CPliris::removePliris( CT_Pliris_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Pliris_ID_t>(*id);
        tableOfPliriss().remove(&aid);
    *id = CTrilinos::concreteType<CT_Pliris_ID_t>(aid);
}

/* remove Pliris from table using CTrilinos_Universal_ID_t */
void
CPliris::removePliris( CTrilinos_Universal_ID_t *aid )
{
        tableOfPliriss().remove(aid);
}

/* purge Pliris table */
void
CPliris::purgePliris(  )
{
    tableOfPliriss().purge();
}


#endif /* HAVE_MPI */

#endif /* HAVE_CTRILINOS_PLIRIS */


