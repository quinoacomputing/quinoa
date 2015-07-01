
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
#include "CEpetra_RowMatrix.h"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "Epetra_RowMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_RowMatrix */
Table<Epetra_RowMatrix>& tableOfRowMatrixs()
{
    static Table<Epetra_RowMatrix> loc_tableOfRowMatrixs(CT_Epetra_RowMatrix_ID);
    return loc_tableOfRowMatrixs;
}


} // namespace


//
// Definitions from CEpetra_RowMatrix.h
//


extern "C" {


CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_RowMatrix_Generalize ( 
  CT_Epetra_RowMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id);
}

void Epetra_RowMatrix_Destroy ( CT_Epetra_RowMatrix_ID_t * selfID )
{
    CEpetra::removeRowMatrix(selfID);
}

int Epetra_RowMatrix_NumMyRowEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyRowEntries(MyRow, 
        *NumEntries);
}

int Epetra_RowMatrix_MaxNumEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->MaxNumEntries();
}

int Epetra_RowMatrix_ExtractMyRowCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CEpetra::getConstRowMatrix(selfID)->ExtractMyRowCopy(MyRow, Length, 
        *NumEntries, Values, Indices);
}

int Epetra_RowMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID )
{
    const Teuchos::RCP<Epetra_Vector> Diagonal = CEpetra::getVector(
        DiagonalID);
    return CEpetra::getConstRowMatrix(selfID)->ExtractDiagonalCopy(*Diagonal);
}

int Epetra_RowMatrix_Multiply ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstRowMatrix(selfID)->Multiply(((TransA) != 
        FALSE ? true : false), *X, *Y);
}

int Epetra_RowMatrix_Solve ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstRowMatrix(selfID)->Solve(((Upper) != 
        FALSE ? true : false), ((Trans) != FALSE ? true : false), ((
        UnitDiagonal) != FALSE ? true : false), *X, *Y);
}

int Epetra_RowMatrix_InvRowSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = CEpetra::getVector(xID);
    return CEpetra::getConstRowMatrix(selfID)->InvRowSums(*x);
}

int Epetra_RowMatrix_LeftScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<const Epetra_Vector> x = CEpetra::getConstVector(xID);
    return CEpetra::getRowMatrix(selfID)->LeftScale(*x);
}

int Epetra_RowMatrix_InvColSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = CEpetra::getVector(xID);
    return CEpetra::getConstRowMatrix(selfID)->InvColSums(*x);
}

int Epetra_RowMatrix_RightScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<const Epetra_Vector> x = CEpetra::getConstVector(xID);
    return CEpetra::getRowMatrix(selfID)->RightScale(*x);
}

boolean Epetra_RowMatrix_Filled ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return ((CEpetra::getConstRowMatrix(selfID)->Filled()) ? TRUE : FALSE);
}

double Epetra_RowMatrix_NormInf ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NormInf();
}

double Epetra_RowMatrix_NormOne ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NormOne();
}

int Epetra_RowMatrix_NumGlobalNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumGlobalNonzeros();
}

int Epetra_RowMatrix_NumGlobalRows ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumGlobalRows();
}

int Epetra_RowMatrix_NumGlobalCols ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumGlobalCols();
}

int Epetra_RowMatrix_NumGlobalDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumGlobalDiagonals();
}

int Epetra_RowMatrix_NumMyNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyNonzeros();
}

int Epetra_RowMatrix_NumMyRows ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyRows();
}

int Epetra_RowMatrix_NumMyCols ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyCols();
}

int Epetra_RowMatrix_NumMyDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyDiagonals();
}

boolean Epetra_RowMatrix_LowerTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return ((CEpetra::getConstRowMatrix(
        selfID)->LowerTriangular()) ? TRUE : FALSE);
}

boolean Epetra_RowMatrix_UpperTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return ((CEpetra::getConstRowMatrix(
        selfID)->UpperTriangular()) ? TRUE : FALSE);
}

CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixRowMap ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstRowMatrix(
        selfID)->RowMatrixRowMap() ));
}

CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixColMap ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstRowMatrix(
        selfID)->RowMatrixColMap() ));
}

CT_Epetra_Import_ID_t Epetra_RowMatrix_RowMatrixImporter ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::storeConstImport(CEpetra::getConstRowMatrix(
        selfID)->RowMatrixImporter());
}


} // extern "C"


//
// Definitions from CEpetra_RowMatrix_Cpp.hpp
//


/* get Epetra_RowMatrix from non-const table using CT_Epetra_RowMatrix_ID */
const Teuchos::RCP<Epetra_RowMatrix>
CEpetra::getRowMatrix( CT_Epetra_RowMatrix_ID_t id )
{
    if (tableOfRowMatrixs().isType(id.table))
        return tableOfRowMatrixs().get<Epetra_RowMatrix>(
        CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_RowMatrix>(
        CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id));
}

/* get Epetra_RowMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_RowMatrix>
CEpetra::getRowMatrix( CTrilinos_Universal_ID_t id )
{
    if (tableOfRowMatrixs().isType(id.table))
        return tableOfRowMatrixs().get<Epetra_RowMatrix>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_RowMatrix>(id);
}

/* get const Epetra_RowMatrix from either the const or non-const table
 * using CT_Epetra_RowMatrix_ID */
const Teuchos::RCP<const Epetra_RowMatrix>
CEpetra::getConstRowMatrix( CT_Epetra_RowMatrix_ID_t id )
{
    if (tableOfRowMatrixs().isType(id.table))
        return tableOfRowMatrixs().getConst<Epetra_RowMatrix>(
        CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_RowMatrix>(
        CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id));
}

/* get const Epetra_RowMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_RowMatrix>
CEpetra::getConstRowMatrix( CTrilinos_Universal_ID_t id )
{
    if (tableOfRowMatrixs().isType(id.table))
        return tableOfRowMatrixs().getConst<Epetra_RowMatrix>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_RowMatrix>(id);
}

/* store Epetra_RowMatrix (owned) in non-const table */
CT_Epetra_RowMatrix_ID_t
CEpetra::storeNewRowMatrix( Epetra_RowMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(
        tableOfRowMatrixs().store<Epetra_RowMatrix>(pobj, true));
}

/* store Epetra_RowMatrix in non-const table */
CT_Epetra_RowMatrix_ID_t
CEpetra::storeRowMatrix( Epetra_RowMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(
        tableOfRowMatrixs().store<Epetra_RowMatrix>(pobj, false));
}

/* store const Epetra_RowMatrix in const table */
CT_Epetra_RowMatrix_ID_t
CEpetra::storeConstRowMatrix( const Epetra_RowMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(
        tableOfRowMatrixs().store<Epetra_RowMatrix>(pobj, false));
}

/* remove Epetra_RowMatrix from table using CT_Epetra_RowMatrix_ID */
void
CEpetra::removeRowMatrix( CT_Epetra_RowMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(*id);
    if (tableOfRowMatrixs().isType(aid.table))
        tableOfRowMatrixs().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(aid);
}

/* remove Epetra_RowMatrix from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeRowMatrix( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfRowMatrixs().isType(aid->table))
        tableOfRowMatrixs().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_RowMatrix table */
void
CEpetra::purgeRowMatrix(  )
{
    tableOfRowMatrixs().purge();
}

/* store Epetra_RowMatrix in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasRowMatrix( const Teuchos::RCP< Epetra_RowMatrix > & robj )
{
    return tableOfRowMatrixs().alias(robj);
}

/* store const Epetra_RowMatrix in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstRowMatrix( const Teuchos::RCP< const Epetra_RowMatrix > & robj )
{
    return tableOfRowMatrixs().alias(robj);
}



