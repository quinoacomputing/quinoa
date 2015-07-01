
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
#include "CEpetra_JadMatrix.h"
#include "CEpetra_JadMatrix_Cpp.hpp"
#include "Epetra_JadMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_JadMatrix */
Table<Epetra_JadMatrix>& tableOfJadMatrixs()
{
    static Table<Epetra_JadMatrix> loc_tableOfJadMatrixs(CT_Epetra_JadMatrix_ID);
    return loc_tableOfJadMatrixs;
}


} // namespace


//
// Definitions from CEpetra_JadMatrix.h
//


extern "C" {


CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_JadMatrix_Generalize ( 
  CT_Epetra_JadMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id);
}

CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Create ( 
  CT_Epetra_RowMatrix_ID_t MatrixID )
{
    const Teuchos::RCP<const Epetra_RowMatrix> Matrix = 
        CEpetra::getConstRowMatrix(MatrixID);
    return CEpetra::storeNewJadMatrix(new Epetra_JadMatrix(*Matrix));
}

void Epetra_JadMatrix_Destroy ( CT_Epetra_JadMatrix_ID_t * selfID )
{
    CEpetra::removeJadMatrix(selfID);
}

int Epetra_JadMatrix_UpdateValues ( 
  CT_Epetra_JadMatrix_ID_t selfID, 
  CT_Epetra_RowMatrix_ID_t MatrixID, boolean CheckStructure )
{
    const Teuchos::RCP<const Epetra_RowMatrix> Matrix = 
        CEpetra::getConstRowMatrix(MatrixID);
    return CEpetra::getJadMatrix(selfID)->UpdateValues(*Matrix, ((
        CheckStructure) != FALSE ? true : false));
}

int Epetra_JadMatrix_ExtractMyRowCopy ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CEpetra::getConstJadMatrix(selfID)->ExtractMyRowCopy(MyRow, Length, 
        *NumEntries, Values, Indices);
}

int Epetra_JadMatrix_ExtractMyEntryView ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, double * * Value, 
  int * RowIndex, int * ColIndex )
{
    return CEpetra::getJadMatrix(selfID)->ExtractMyEntryView(CurEntry, *Value, 
        *RowIndex, *ColIndex);
}

int Epetra_JadMatrix_ExtractMyEntryView_Const ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, 
  double const ** Value, int * RowIndex, int * ColIndex )
{
    return CEpetra::getConstJadMatrix(selfID)->ExtractMyEntryView(CurEntry, 
        *Value, *RowIndex, *ColIndex);
}

int Epetra_JadMatrix_NumMyRowEntries ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CEpetra::getConstJadMatrix(selfID)->NumMyRowEntries(MyRow, 
        *NumEntries);
}

int Epetra_JadMatrix_Multiply ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstJadMatrix(selfID)->Multiply(((TransA) != 
        FALSE ? true : false), *X, *Y);
}

int Epetra_JadMatrix_Solve ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstJadMatrix(selfID)->Solve(((Upper) != 
        FALSE ? true : false), ((Trans) != FALSE ? true : false), ((
        UnitDiagonal) != FALSE ? true : false), *X, *Y);
}


} // extern "C"


//
// Definitions from CEpetra_JadMatrix_Cpp.hpp
//


/* get Epetra_JadMatrix from non-const table using CT_Epetra_JadMatrix_ID */
const Teuchos::RCP<Epetra_JadMatrix>
CEpetra::getJadMatrix( CT_Epetra_JadMatrix_ID_t id )
{
    if (tableOfJadMatrixs().isType(id.table))
        return tableOfJadMatrixs().get<Epetra_JadMatrix>(
        CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_JadMatrix>(
        CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id));
}

/* get Epetra_JadMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_JadMatrix>
CEpetra::getJadMatrix( CTrilinos_Universal_ID_t id )
{
    if (tableOfJadMatrixs().isType(id.table))
        return tableOfJadMatrixs().get<Epetra_JadMatrix>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_JadMatrix>(id);
}

/* get const Epetra_JadMatrix from either the const or non-const table
 * using CT_Epetra_JadMatrix_ID */
const Teuchos::RCP<const Epetra_JadMatrix>
CEpetra::getConstJadMatrix( CT_Epetra_JadMatrix_ID_t id )
{
    if (tableOfJadMatrixs().isType(id.table))
        return tableOfJadMatrixs().getConst<Epetra_JadMatrix>(
        CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_JadMatrix>(
        CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id));
}

/* get const Epetra_JadMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_JadMatrix>
CEpetra::getConstJadMatrix( CTrilinos_Universal_ID_t id )
{
    if (tableOfJadMatrixs().isType(id.table))
        return tableOfJadMatrixs().getConst<Epetra_JadMatrix>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_JadMatrix>(id);
}

/* store Epetra_JadMatrix (owned) in non-const table */
CT_Epetra_JadMatrix_ID_t
CEpetra::storeNewJadMatrix( Epetra_JadMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(
        tableOfJadMatrixs().store<Epetra_JadMatrix>(pobj, true));
}

/* store Epetra_JadMatrix in non-const table */
CT_Epetra_JadMatrix_ID_t
CEpetra::storeJadMatrix( Epetra_JadMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(
        tableOfJadMatrixs().store<Epetra_JadMatrix>(pobj, false));
}

/* store const Epetra_JadMatrix in const table */
CT_Epetra_JadMatrix_ID_t
CEpetra::storeConstJadMatrix( const Epetra_JadMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(
        tableOfJadMatrixs().store<Epetra_JadMatrix>(pobj, false));
}

/* remove Epetra_JadMatrix from table using CT_Epetra_JadMatrix_ID */
void
CEpetra::removeJadMatrix( CT_Epetra_JadMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(*id);
    if (tableOfJadMatrixs().isType(aid.table))
        tableOfJadMatrixs().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(aid);
}

/* remove Epetra_JadMatrix from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeJadMatrix( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfJadMatrixs().isType(aid->table))
        tableOfJadMatrixs().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_JadMatrix table */
void
CEpetra::purgeJadMatrix(  )
{
    tableOfJadMatrixs().purge();
}

/* store Epetra_JadMatrix in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasJadMatrix( const Teuchos::RCP< Epetra_JadMatrix > & robj )
{
    return tableOfJadMatrixs().alias(robj);
}

/* store const Epetra_JadMatrix in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstJadMatrix( const Teuchos::RCP< const Epetra_JadMatrix > & robj )
{
    return tableOfJadMatrixs().alias(robj);
}



