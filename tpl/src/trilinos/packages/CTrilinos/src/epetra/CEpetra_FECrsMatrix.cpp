
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
#include "CEpetra_FECrsMatrix.h"
#include "CEpetra_FECrsMatrix_Cpp.hpp"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_IntSerialDenseVector_Cpp.hpp"
#include "CEpetra_SerialDenseMatrix_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_FECrsMatrix */
Table<Epetra_FECrsMatrix>& tableOfFECrsMatrixs()
{
    static Table<Epetra_FECrsMatrix> loc_tableOfFECrsMatrixs(CT_Epetra_FECrsMatrix_ID);
    return loc_tableOfFECrsMatrixs;
}


} // namespace


//
// Definitions from CEpetra_FECrsMatrix.h
//


extern "C" {


CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_FECrsMatrix_Generalize ( 
  CT_Epetra_FECrsMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(id);
}

CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_Var ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  int * NumEntriesPerRow, boolean ignoreNonLocalEntries )
{
    const Teuchos::RCP<const Epetra_Map> RowMap = CEpetra::getConstMap(
        RowMapID);
    return CEpetra::storeNewFECrsMatrix(new Epetra_FECrsMatrix(
        (Epetra_DataAccess) CV, *RowMap, NumEntriesPerRow, ((
        ignoreNonLocalEntries) != FALSE ? true : false)));
}

CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  int NumEntriesPerRow, boolean ignoreNonLocalEntries )
{
    const Teuchos::RCP<const Epetra_Map> RowMap = CEpetra::getConstMap(
        RowMapID);
    return CEpetra::storeNewFECrsMatrix(new Epetra_FECrsMatrix(
        (Epetra_DataAccess) CV, *RowMap, NumEntriesPerRow, ((
        ignoreNonLocalEntries) != FALSE ? true : false)));
}

CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_WithColMap_Var ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int * NumEntriesPerRow, 
  boolean ignoreNonLocalEntries )
{
    const Teuchos::RCP<const Epetra_Map> RowMap = CEpetra::getConstMap(
        RowMapID);
    const Teuchos::RCP<const Epetra_Map> ColMap = CEpetra::getConstMap(
        ColMapID);
    return CEpetra::storeNewFECrsMatrix(new Epetra_FECrsMatrix(
        (Epetra_DataAccess) CV, *RowMap, *ColMap, NumEntriesPerRow, ((
        ignoreNonLocalEntries) != FALSE ? true : false)));
}

CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, 
  boolean ignoreNonLocalEntries )
{
    const Teuchos::RCP<const Epetra_Map> RowMap = CEpetra::getConstMap(
        RowMapID);
    const Teuchos::RCP<const Epetra_Map> ColMap = CEpetra::getConstMap(
        ColMapID);
    return CEpetra::storeNewFECrsMatrix(new Epetra_FECrsMatrix(
        (Epetra_DataAccess) CV, *RowMap, *ColMap, NumEntriesPerRow, ((
        ignoreNonLocalEntries) != FALSE ? true : false)));
}

CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_FromGraph ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_CrsGraph_ID_t GraphID, 
  boolean ignoreNonLocalEntries )
{
    const Teuchos::RCP<const Epetra_CrsGraph> Graph = CEpetra::getConstCrsGraph(
        GraphID);
    return CEpetra::storeNewFECrsMatrix(new Epetra_FECrsMatrix(
        (Epetra_DataAccess) CV, *Graph, ((ignoreNonLocalEntries) != 
        FALSE ? true : false)));
}

CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Duplicate ( 
  CT_Epetra_FECrsMatrix_ID_t srcID )
{
    const Teuchos::RCP<const Epetra_FECrsMatrix> src = 
        CEpetra::getConstFECrsMatrix(srcID);
    return CEpetra::storeNewFECrsMatrix(new Epetra_FECrsMatrix(*src));
}

void Epetra_FECrsMatrix_Destroy ( 
  CT_Epetra_FECrsMatrix_ID_t * selfID )
{
    CEpetra::removeFECrsMatrix(selfID);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(GlobalRow, 
        NumEntries, Values, Indices);
}

int Epetra_FECrsMatrix_InsertGlobalValues ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(GlobalRow, 
        NumEntries, Values, Indices);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(GlobalRow, 
        NumEntries, Values, Indices);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues_Ftable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(numIndices, 
        indices, values, format);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues_Ftable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(numRows, rows, 
        numCols, cols, values, format);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues_Ctable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double* const * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(numIndices, 
        indices, values, format);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues_Ctable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double* const * values, 
  int format )
{
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(numRows, rows, 
        numCols, cols, values, format);
}

int Epetra_FECrsMatrix_InsertGlobalValues_Ftable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(numIndices, 
        indices, values, format);
}

int Epetra_FECrsMatrix_InsertGlobalValues_Ftable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(numRows, rows, 
        numCols, cols, values, format);
}

int Epetra_FECrsMatrix_InsertGlobalValues_Ctable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double* const * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(numIndices, 
        indices, values, format);
}

int Epetra_FECrsMatrix_InsertGlobalValues_Ctable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double* const * values, 
  int format )
{
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(numRows, rows, 
        numCols, cols, values, format);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues_Ftable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(numIndices, 
        indices, values, format);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues_Ftable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(numRows, rows, 
        numCols, cols, values, format);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues_Ctable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double* const * values, int format )
{
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(numIndices, 
        indices, values, format);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues_Ctable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double* const * values, 
  int format )
{
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(numRows, rows, 
        numCols, cols, values, format);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues_SubMatrix_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t indicesID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format )
{
    const Teuchos::RCP<const Epetra_IntSerialDenseVector> indices = 
        CEpetra::getConstIntSerialDenseVector(indicesID);
    const Teuchos::RCP<const Epetra_SerialDenseMatrix> values = 
        CEpetra::getConstSerialDenseMatrix(valuesID);
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(*indices, 
        *values, format);
}

int Epetra_FECrsMatrix_SumIntoGlobalValues_SubMatrix ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t rowsID, 
  CT_Epetra_IntSerialDenseVector_ID_t colsID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format )
{
    const Teuchos::RCP<const Epetra_IntSerialDenseVector> rows = 
        CEpetra::getConstIntSerialDenseVector(rowsID);
    const Teuchos::RCP<const Epetra_IntSerialDenseVector> cols = 
        CEpetra::getConstIntSerialDenseVector(colsID);
    const Teuchos::RCP<const Epetra_SerialDenseMatrix> values = 
        CEpetra::getConstSerialDenseMatrix(valuesID);
    return CEpetra::getFECrsMatrix(selfID)->SumIntoGlobalValues(*rows, *cols, 
        *values, format);
}

int Epetra_FECrsMatrix_InsertGlobalValues_SubMatrix_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t indicesID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format )
{
    const Teuchos::RCP<const Epetra_IntSerialDenseVector> indices = 
        CEpetra::getConstIntSerialDenseVector(indicesID);
    const Teuchos::RCP<const Epetra_SerialDenseMatrix> values = 
        CEpetra::getConstSerialDenseMatrix(valuesID);
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(*indices, 
        *values, format);
}

int Epetra_FECrsMatrix_InsertGlobalValues_SubMatrix ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t rowsID, 
  CT_Epetra_IntSerialDenseVector_ID_t colsID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format )
{
    const Teuchos::RCP<const Epetra_IntSerialDenseVector> rows = 
        CEpetra::getConstIntSerialDenseVector(rowsID);
    const Teuchos::RCP<const Epetra_IntSerialDenseVector> cols = 
        CEpetra::getConstIntSerialDenseVector(colsID);
    const Teuchos::RCP<const Epetra_SerialDenseMatrix> values = 
        CEpetra::getConstSerialDenseMatrix(valuesID);
    return CEpetra::getFECrsMatrix(selfID)->InsertGlobalValues(*rows, *cols, 
        *values, format);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues_SubMatrix_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t indicesID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format )
{
    const Teuchos::RCP<const Epetra_IntSerialDenseVector> indices = 
        CEpetra::getConstIntSerialDenseVector(indicesID);
    const Teuchos::RCP<const Epetra_SerialDenseMatrix> values = 
        CEpetra::getConstSerialDenseMatrix(valuesID);
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(*indices, 
        *values, format);
}

int Epetra_FECrsMatrix_ReplaceGlobalValues_SubMatrix ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t rowsID, 
  CT_Epetra_IntSerialDenseVector_ID_t colsID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format )
{
    const Teuchos::RCP<const Epetra_IntSerialDenseVector> rows = 
        CEpetra::getConstIntSerialDenseVector(rowsID);
    const Teuchos::RCP<const Epetra_IntSerialDenseVector> cols = 
        CEpetra::getConstIntSerialDenseVector(colsID);
    const Teuchos::RCP<const Epetra_SerialDenseMatrix> values = 
        CEpetra::getConstSerialDenseMatrix(valuesID);
    return CEpetra::getFECrsMatrix(selfID)->ReplaceGlobalValues(*rows, *cols, 
        *values, format);
}

int Epetra_FECrsMatrix_GlobalAssemble ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, boolean callFillComplete )
{
    return CEpetra::getFECrsMatrix(selfID)->GlobalAssemble(
        ((callFillComplete) != FALSE ? true : false));
}

int Epetra_FECrsMatrix_GlobalAssemble_WithMaps ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_Map_ID_t domain_mapID, CT_Epetra_Map_ID_t range_mapID, 
  boolean callFillComplete )
{
    const Teuchos::RCP<const Epetra_Map> domain_map = CEpetra::getConstMap(
        domain_mapID);
    const Teuchos::RCP<const Epetra_Map> range_map = CEpetra::getConstMap(
        range_mapID);
    return CEpetra::getFECrsMatrix(selfID)->GlobalAssemble(*domain_map, 
        *range_map, ((callFillComplete) != FALSE ? true : false));
}

void Epetra_FECrsMatrix_setIgnoreNonLocalEntries ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, boolean flag )
{
    CEpetra::getFECrsMatrix(selfID)->setIgnoreNonLocalEntries(((flag) != 
        FALSE ? true : false));
}

void Epetra_FECrsMatrix_Assign ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_FECrsMatrix_ID_t srcID )
{
    Epetra_FECrsMatrix& self = *( CEpetra::getFECrsMatrix(selfID) );

    const Teuchos::RCP<const Epetra_FECrsMatrix> src = 
        CEpetra::getConstFECrsMatrix(srcID);
    self = *src;
}


} // extern "C"


//
// Definitions from CEpetra_FECrsMatrix_Cpp.hpp
//


/* get Epetra_FECrsMatrix from non-const table using CT_Epetra_FECrsMatrix_ID */
const Teuchos::RCP<Epetra_FECrsMatrix>
CEpetra::getFECrsMatrix( CT_Epetra_FECrsMatrix_ID_t id )
{
    if (tableOfFECrsMatrixs().isType(id.table))
        return tableOfFECrsMatrixs().get<Epetra_FECrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_FECrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(id));
}

/* get Epetra_FECrsMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_FECrsMatrix>
CEpetra::getFECrsMatrix( CTrilinos_Universal_ID_t id )
{
    if (tableOfFECrsMatrixs().isType(id.table))
        return tableOfFECrsMatrixs().get<Epetra_FECrsMatrix>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_FECrsMatrix>(id);
}

/* get const Epetra_FECrsMatrix from either the const or non-const table
 * using CT_Epetra_FECrsMatrix_ID */
const Teuchos::RCP<const Epetra_FECrsMatrix>
CEpetra::getConstFECrsMatrix( CT_Epetra_FECrsMatrix_ID_t id )
{
    if (tableOfFECrsMatrixs().isType(id.table))
        return tableOfFECrsMatrixs().getConst<Epetra_FECrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_FECrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(id));
}

/* get const Epetra_FECrsMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_FECrsMatrix>
CEpetra::getConstFECrsMatrix( CTrilinos_Universal_ID_t id )
{
    if (tableOfFECrsMatrixs().isType(id.table))
        return tableOfFECrsMatrixs().getConst<Epetra_FECrsMatrix>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_FECrsMatrix>(id);
}

/* store Epetra_FECrsMatrix (owned) in non-const table */
CT_Epetra_FECrsMatrix_ID_t
CEpetra::storeNewFECrsMatrix( Epetra_FECrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        tableOfFECrsMatrixs().store<Epetra_FECrsMatrix>(pobj, true));
}

/* store Epetra_FECrsMatrix in non-const table */
CT_Epetra_FECrsMatrix_ID_t
CEpetra::storeFECrsMatrix( Epetra_FECrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        tableOfFECrsMatrixs().store<Epetra_FECrsMatrix>(pobj, false));
}

/* store const Epetra_FECrsMatrix in const table */
CT_Epetra_FECrsMatrix_ID_t
CEpetra::storeConstFECrsMatrix( const Epetra_FECrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        tableOfFECrsMatrixs().store<Epetra_FECrsMatrix>(pobj, false));
}

/* remove Epetra_FECrsMatrix from table using CT_Epetra_FECrsMatrix_ID */
void
CEpetra::removeFECrsMatrix( CT_Epetra_FECrsMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(*id);
    if (tableOfFECrsMatrixs().isType(aid.table))
        tableOfFECrsMatrixs().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(aid);
}

/* remove Epetra_FECrsMatrix from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeFECrsMatrix( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfFECrsMatrixs().isType(aid->table))
        tableOfFECrsMatrixs().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_FECrsMatrix table */
void
CEpetra::purgeFECrsMatrix(  )
{
    tableOfFECrsMatrixs().purge();
}

/* store Epetra_FECrsMatrix in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasFECrsMatrix( const Teuchos::RCP< Epetra_FECrsMatrix > & robj )
{
    return tableOfFECrsMatrixs().alias(robj);
}

/* store const Epetra_FECrsMatrix in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstFECrsMatrix( const Teuchos::RCP< const Epetra_FECrsMatrix > & robj )
{
    return tableOfFECrsMatrixs().alias(robj);
}



