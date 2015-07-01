
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
#include "CEpetra_CrsMatrix.h"
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_CrsMatrix */
Table<Epetra_CrsMatrix>& tableOfCrsMatrixs()
{
    static Table<Epetra_CrsMatrix> loc_tableOfCrsMatrixs(CT_Epetra_CrsMatrix_ID);
    return loc_tableOfCrsMatrixs;
}


} // namespace


//
// Definitions from CEpetra_CrsMatrix.h
//


extern "C" {


CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_CrsMatrix_Generalize ( 
  CT_Epetra_CrsMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id);
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  const int * NumEntriesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_Map> RowMap = CEpetra::getConstMap(
        RowMapID);
    return CEpetra::storeNewCrsMatrix(new Epetra_CrsMatrix(
        (Epetra_DataAccess) CV, *RowMap, NumEntriesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  int NumEntriesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_Map> RowMap = CEpetra::getConstMap(
        RowMapID);
    return CEpetra::storeNewCrsMatrix(new Epetra_CrsMatrix(
        (Epetra_DataAccess) CV, *RowMap, NumEntriesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, const int * NumEntriesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_Map> RowMap = CEpetra::getConstMap(
        RowMapID);
    const Teuchos::RCP<const Epetra_Map> ColMap = CEpetra::getConstMap(
        ColMapID);
    return CEpetra::storeNewCrsMatrix(new Epetra_CrsMatrix(
        (Epetra_DataAccess) CV, *RowMap, *ColMap, NumEntriesPerRow, ((
        StaticProfile) != FALSE ? true : false)));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_Map> RowMap = CEpetra::getConstMap(
        RowMapID);
    const Teuchos::RCP<const Epetra_Map> ColMap = CEpetra::getConstMap(
        ColMapID);
    return CEpetra::storeNewCrsMatrix(new Epetra_CrsMatrix(
        (Epetra_DataAccess) CV, *RowMap, *ColMap, NumEntriesPerRow, ((
        StaticProfile) != FALSE ? true : false)));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_FromGraph ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_CrsGraph_ID_t GraphID )
{
    const Teuchos::RCP<const Epetra_CrsGraph> Graph = CEpetra::getConstCrsGraph(
        GraphID);
    return CEpetra::storeNewCrsMatrix(new Epetra_CrsMatrix(
        (Epetra_DataAccess) CV, *Graph));
}

CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Duplicate ( 
  CT_Epetra_CrsMatrix_ID_t MatrixID )
{
    const Teuchos::RCP<const Epetra_CrsMatrix> Matrix = 
        CEpetra::getConstCrsMatrix(MatrixID);
    return CEpetra::storeNewCrsMatrix(new Epetra_CrsMatrix(*Matrix));
}

void Epetra_CrsMatrix_Destroy ( CT_Epetra_CrsMatrix_ID_t * selfID )
{
    CEpetra::removeCrsMatrix(selfID);
}

int Epetra_CrsMatrix_PutScalar ( 
  CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant )
{
    return CEpetra::getCrsMatrix(selfID)->PutScalar(ScalarConstant);
}

int Epetra_CrsMatrix_Scale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant )
{
    return CEpetra::getCrsMatrix(selfID)->Scale(ScalarConstant);
}

int Epetra_CrsMatrix_InsertGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->InsertGlobalValues(GlobalRow, 
        NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ReplaceGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->ReplaceGlobalValues(GlobalRow, 
        NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_SumIntoGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->SumIntoGlobalValues(GlobalRow, 
        NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_InsertMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->InsertMyValues(MyRow, NumEntries, 
        Values, Indices);
}

int Epetra_CrsMatrix_ReplaceMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->ReplaceMyValues(MyRow, NumEntries, 
        Values, Indices);
}

int Epetra_CrsMatrix_SumIntoMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices )
{
    return CEpetra::getCrsMatrix(selfID)->SumIntoMyValues(MyRow, NumEntries, 
        Values, Indices);
}

int Epetra_CrsMatrix_ReplaceDiagonalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID )
{
    const Teuchos::RCP<const Epetra_Vector> Diagonal = CEpetra::getConstVector(
        DiagonalID);
    return CEpetra::getCrsMatrix(selfID)->ReplaceDiagonalValues(*Diagonal);
}

int Epetra_CrsMatrix_FillComplete ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean OptimizeDataStorage )
{
    return CEpetra::getCrsMatrix(selfID)->FillComplete(
        ((OptimizeDataStorage) != FALSE ? true : false));
}

int Epetra_CrsMatrix_FillComplete_UsingMaps ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Map_ID_t DomainMapID, 
  CT_Epetra_Map_ID_t RangeMapID, boolean OptimizeDataStorage )
{
    const Teuchos::RCP<const Epetra_Map> DomainMap = CEpetra::getConstMap(
        DomainMapID);
    const Teuchos::RCP<const Epetra_Map> RangeMap = CEpetra::getConstMap(
        RangeMapID);
    return CEpetra::getCrsMatrix(selfID)->FillComplete(*DomainMap, *RangeMap, ((
        OptimizeDataStorage) != FALSE ? true : false));
}

int Epetra_CrsMatrix_OptimizeStorage ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->OptimizeStorage();
}

int Epetra_CrsMatrix_MakeDataContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->MakeDataContiguous();
}

int Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractGlobalRowCopy(GlobalRow, 
        Length, *NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ExtractMyRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractMyRowCopy(MyRow, Length, 
        *NumEntries, Values, Indices);
}

int Epetra_CrsMatrix_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractGlobalRowCopy(GlobalRow, 
        Length, *NumEntries, Values);
}

int Epetra_CrsMatrix_ExtractMyRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractMyRowCopy(MyRow, Length, 
        *NumEntries, Values);
}

int Epetra_CrsMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID )
{
    const Teuchos::RCP<Epetra_Vector> Diagonal = CEpetra::getVector(
        DiagonalID);
    return CEpetra::getConstCrsMatrix(selfID)->ExtractDiagonalCopy(*Diagonal);
}

int Epetra_CrsMatrix_ExtractGlobalRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values, int ** Indices )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractGlobalRowView(GlobalRow, 
        *NumEntries, *Values, *Indices);
}

int Epetra_CrsMatrix_ExtractMyRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values, int ** Indices )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractMyRowView(MyRow, 
        *NumEntries, *Values, *Indices);
}

int Epetra_CrsMatrix_ExtractGlobalRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractGlobalRowView(GlobalRow, 
        *NumEntries, *Values);
}

int Epetra_CrsMatrix_ExtractMyRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values )
{
    return CEpetra::getConstCrsMatrix(selfID)->ExtractMyRowView(MyRow, 
        *NumEntries, *Values);
}

int Epetra_CrsMatrix_Multiply_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID )
{
    const Teuchos::RCP<const Epetra_Vector> x = CEpetra::getConstVector(xID);
    const Teuchos::RCP<Epetra_Vector> y = CEpetra::getVector(yID);
    return CEpetra::getConstCrsMatrix(selfID)->Multiply(((TransA) != 
        FALSE ? true : false), *x, *y);
}

int Epetra_CrsMatrix_Multiply1_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID )
{
    const Teuchos::RCP<const Epetra_Vector> x = CEpetra::getConstVector(xID);
    const Teuchos::RCP<Epetra_Vector> y = CEpetra::getVector(yID);
    return CEpetra::getConstCrsMatrix(selfID)->Multiply1(((TransA) != 
        FALSE ? true : false), *x, *y);
}

int Epetra_CrsMatrix_Multiply_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstCrsMatrix(selfID)->Multiply(((TransA) != 
        FALSE ? true : false), *X, *Y);
}

int Epetra_CrsMatrix_Multiply1_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstCrsMatrix(selfID)->Multiply1(((TransA) != 
        FALSE ? true : false), *X, *Y);
}

int Epetra_CrsMatrix_Solve_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_Vector_ID_t xID, 
  CT_Epetra_Vector_ID_t yID )
{
    const Teuchos::RCP<const Epetra_Vector> x = CEpetra::getConstVector(xID);
    const Teuchos::RCP<Epetra_Vector> y = CEpetra::getVector(yID);
    return CEpetra::getConstCrsMatrix(selfID)->Solve(((Upper) != 
        FALSE ? true : false), ((Trans) != FALSE ? true : false), ((
        UnitDiagonal) != FALSE ? true : false), *x, *y);
}

int Epetra_CrsMatrix_Solve_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstCrsMatrix(selfID)->Solve(((Upper) != 
        FALSE ? true : false), ((Trans) != FALSE ? true : false), ((
        UnitDiagonal) != FALSE ? true : false), *X, *Y);
}

int Epetra_CrsMatrix_InvRowSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = CEpetra::getVector(xID);
    return CEpetra::getConstCrsMatrix(selfID)->InvRowSums(*x);
}

int Epetra_CrsMatrix_InvRowMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = CEpetra::getVector(xID);
    return CEpetra::getConstCrsMatrix(selfID)->InvRowMaxs(*x);
}

int Epetra_CrsMatrix_LeftScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<const Epetra_Vector> x = CEpetra::getConstVector(xID);
    return CEpetra::getCrsMatrix(selfID)->LeftScale(*x);
}

int Epetra_CrsMatrix_InvColSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = CEpetra::getVector(xID);
    return CEpetra::getConstCrsMatrix(selfID)->InvColSums(*x);
}

int Epetra_CrsMatrix_InvColMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = CEpetra::getVector(xID);
    return CEpetra::getConstCrsMatrix(selfID)->InvColMaxs(*x);
}

int Epetra_CrsMatrix_RightScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<const Epetra_Vector> x = CEpetra::getConstVector(xID);
    return CEpetra::getCrsMatrix(selfID)->RightScale(*x);
}

boolean Epetra_CrsMatrix_Filled ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CEpetra::getConstCrsMatrix(selfID)->Filled()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_StorageOptimized ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CEpetra::getConstCrsMatrix(
        selfID)->StorageOptimized()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_IndicesAreGlobal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CEpetra::getConstCrsMatrix(
        selfID)->IndicesAreGlobal()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_IndicesAreLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CEpetra::getConstCrsMatrix(
        selfID)->IndicesAreLocal()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_IndicesAreContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CEpetra::getConstCrsMatrix(
        selfID)->IndicesAreContiguous()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_LowerTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CEpetra::getConstCrsMatrix(
        selfID)->LowerTriangular()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_UpperTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CEpetra::getConstCrsMatrix(
        selfID)->UpperTriangular()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_NoDiagonal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CEpetra::getConstCrsMatrix(selfID)->NoDiagonal()) ? TRUE : FALSE);
}

double Epetra_CrsMatrix_NormInf ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NormInf();
}

double Epetra_CrsMatrix_NormOne ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NormOne();
}

double Epetra_CrsMatrix_NormFrobenius ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NormFrobenius();
}

int Epetra_CrsMatrix_NumGlobalNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumGlobalNonzeros();
}

int Epetra_CrsMatrix_NumGlobalRows ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumGlobalRows();
}

int Epetra_CrsMatrix_NumGlobalCols ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumGlobalCols();
}

int Epetra_CrsMatrix_NumGlobalDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumGlobalDiagonals();
}

int Epetra_CrsMatrix_NumMyNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumMyNonzeros();
}

int Epetra_CrsMatrix_NumMyRows ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumMyRows();
}

int Epetra_CrsMatrix_NumMyCols ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumMyCols();
}

int Epetra_CrsMatrix_NumMyDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumMyDiagonals();
}

int Epetra_CrsMatrix_NumGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumGlobalEntries(Row);
}

int Epetra_CrsMatrix_NumAllocatedGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumAllocatedGlobalEntries(Row);
}

int Epetra_CrsMatrix_MaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->MaxNumEntries();
}

int Epetra_CrsMatrix_GlobalMaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->GlobalMaxNumEntries();
}

int Epetra_CrsMatrix_NumMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumMyEntries(Row);
}

int Epetra_CrsMatrix_NumAllocatedMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumAllocatedMyEntries(Row);
}

int Epetra_CrsMatrix_IndexBase ( CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->IndexBase();
}

boolean Epetra_CrsMatrix_StaticGraph ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CEpetra::getCrsMatrix(selfID)->StaticGraph()) ? TRUE : FALSE);
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsMatrix_Graph ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstCrsGraph(&( CEpetra::getConstCrsMatrix(
        selfID)->Graph() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstCrsMatrix(
        selfID)->RowMap() ));
}

int Epetra_CrsMatrix_ReplaceRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> newmap = 
        CEpetra::getConstBlockMap(newmapID);
    return CEpetra::getCrsMatrix(selfID)->ReplaceRowMap(*newmap);
}

boolean Epetra_CrsMatrix_HaveColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CEpetra::getConstCrsMatrix(selfID)->HaveColMap()) ? TRUE : FALSE);
}

int Epetra_CrsMatrix_ReplaceColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> newmap = 
        CEpetra::getConstBlockMap(newmapID);
    return CEpetra::getCrsMatrix(selfID)->ReplaceColMap(*newmap);
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_ColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstCrsMatrix(
        selfID)->ColMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_DomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstCrsMatrix(
        selfID)->DomainMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstCrsMatrix(
        selfID)->RangeMap() ));
}

CT_Epetra_Import_ID_t Epetra_CrsMatrix_Importer ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstImport(CEpetra::getConstCrsMatrix(
        selfID)->Importer());
}

CT_Epetra_Export_ID_t Epetra_CrsMatrix_Exporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstExport(CEpetra::getConstCrsMatrix(
        selfID)->Exporter());
}

CT_Epetra_Comm_ID_t Epetra_CrsMatrix_Comm ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstComm(&( CEpetra::getConstCrsMatrix(
        selfID)->Comm() ));
}

int Epetra_CrsMatrix_LRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in )
{
    return CEpetra::getConstCrsMatrix(selfID)->LRID(GRID_in);
}

int Epetra_CrsMatrix_GRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in )
{
    return CEpetra::getConstCrsMatrix(selfID)->GRID(LRID_in);
}

int Epetra_CrsMatrix_LCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in )
{
    return CEpetra::getConstCrsMatrix(selfID)->LCID(GCID_in);
}

int Epetra_CrsMatrix_GCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in )
{
    return CEpetra::getConstCrsMatrix(selfID)->GCID(LCID_in);
}

boolean Epetra_CrsMatrix_MyGRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in )
{
    return ((CEpetra::getConstCrsMatrix(selfID)->MyGRID(
        GRID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_MyLRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in )
{
    return ((CEpetra::getConstCrsMatrix(selfID)->MyLRID(
        LRID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_MyGCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in )
{
    return ((CEpetra::getConstCrsMatrix(selfID)->MyGCID(
        GCID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_MyLCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in )
{
    return ((CEpetra::getConstCrsMatrix(selfID)->MyLCID(
        LCID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_MyGlobalRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GID )
{
    return ((CEpetra::getConstCrsMatrix(selfID)->MyGlobalRow(
        GID)) ? TRUE : FALSE);
}

const char * Epetra_CrsMatrix_Label ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getConstCrsMatrix(selfID)->Label();
}

int Epetra_CrsMatrix_SetUseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean UseTranspose_in )
{
    return CEpetra::getCrsMatrix(selfID)->SetUseTranspose(((UseTranspose_in) != 
        FALSE ? true : false));
}

int Epetra_CrsMatrix_Apply ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstCrsMatrix(selfID)->Apply(*X, *Y);
}

int Epetra_CrsMatrix_ApplyInverse ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstCrsMatrix(selfID)->ApplyInverse(*X, *Y);
}

boolean Epetra_CrsMatrix_HasNormInf ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CEpetra::getConstCrsMatrix(selfID)->HasNormInf()) ? TRUE : FALSE);
}

boolean Epetra_CrsMatrix_UseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return ((CEpetra::getConstCrsMatrix(
        selfID)->UseTranspose()) ? TRUE : FALSE);
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorDomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstCrsMatrix(
        selfID)->OperatorDomainMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorRangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstCrsMatrix(
        selfID)->OperatorRangeMap() ));
}

int Epetra_CrsMatrix_NumMyRowEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CEpetra::getConstCrsMatrix(selfID)->NumMyRowEntries(MyRow, 
        *NumEntries);
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstCrsMatrix(
        selfID)->RowMatrixRowMap() ));
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstCrsMatrix(
        selfID)->RowMatrixColMap() ));
}

CT_Epetra_Import_ID_t Epetra_CrsMatrix_RowMatrixImporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstImport(CEpetra::getConstCrsMatrix(
        selfID)->RowMatrixImporter());
}

CT_Epetra_Map_ID_t Epetra_CrsMatrix_ImportMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstCrsMatrix(
        selfID)->ImportMap() ));
}

int Epetra_CrsMatrix_TransformToLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID )
{
    return CEpetra::getCrsMatrix(selfID)->TransformToLocal();
}

int Epetra_CrsMatrix_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Map_ID_t DomainMapID, 
  CT_Epetra_Map_ID_t RangeMapID )
{
    const Teuchos::RCP<const Epetra_Map> DomainMap = CEpetra::getConstMap(
        DomainMapID);
    const Teuchos::RCP<const Epetra_Map> RangeMap = CEpetra::getConstMap(
        RangeMapID);
    return CEpetra::getCrsMatrix(selfID)->TransformToLocal(
        DomainMap.getRawPtr(), RangeMap.getRawPtr());
}

void Epetra_CrsMatrix_Assign ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_CrsMatrix_ID_t srcID )
{
    Epetra_CrsMatrix& self = *( CEpetra::getCrsMatrix(selfID) );

    const Teuchos::RCP<const Epetra_CrsMatrix> src = CEpetra::getConstCrsMatrix(
        srcID);
    self = *src;
}

double * Epetra_CrsMatrix_getRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Loc )
{
    const Epetra_CrsMatrix& self = *( CEpetra::getConstCrsMatrix(selfID) );

    return self[Loc];
}


} // extern "C"


//
// Definitions from CEpetra_CrsMatrix_Cpp.hpp
//


/* get Epetra_CrsMatrix from non-const table using CT_Epetra_CrsMatrix_ID */
const Teuchos::RCP<Epetra_CrsMatrix>
CEpetra::getCrsMatrix( CT_Epetra_CrsMatrix_ID_t id )
{
    if (tableOfCrsMatrixs().isType(id.table))
        return tableOfCrsMatrixs().get<Epetra_CrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_CrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id));
}

/* get Epetra_CrsMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_CrsMatrix>
CEpetra::getCrsMatrix( CTrilinos_Universal_ID_t id )
{
    if (tableOfCrsMatrixs().isType(id.table))
        return tableOfCrsMatrixs().get<Epetra_CrsMatrix>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_CrsMatrix>(id);
}

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CT_Epetra_CrsMatrix_ID */
const Teuchos::RCP<const Epetra_CrsMatrix>
CEpetra::getConstCrsMatrix( CT_Epetra_CrsMatrix_ID_t id )
{
    if (tableOfCrsMatrixs().isType(id.table))
        return tableOfCrsMatrixs().getConst<Epetra_CrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_CrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id));
}

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_CrsMatrix>
CEpetra::getConstCrsMatrix( CTrilinos_Universal_ID_t id )
{
    if (tableOfCrsMatrixs().isType(id.table))
        return tableOfCrsMatrixs().getConst<Epetra_CrsMatrix>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_CrsMatrix>(id);
}

/* store Epetra_CrsMatrix (owned) in non-const table */
CT_Epetra_CrsMatrix_ID_t
CEpetra::storeNewCrsMatrix( Epetra_CrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store<Epetra_CrsMatrix>(pobj, true));
}

/* store Epetra_CrsMatrix in non-const table */
CT_Epetra_CrsMatrix_ID_t
CEpetra::storeCrsMatrix( Epetra_CrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store<Epetra_CrsMatrix>(pobj, false));
}

/* store const Epetra_CrsMatrix in const table */
CT_Epetra_CrsMatrix_ID_t
CEpetra::storeConstCrsMatrix( const Epetra_CrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tableOfCrsMatrixs().store<Epetra_CrsMatrix>(pobj, false));
}

/* remove Epetra_CrsMatrix from table using CT_Epetra_CrsMatrix_ID */
void
CEpetra::removeCrsMatrix( CT_Epetra_CrsMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(*id);
    if (tableOfCrsMatrixs().isType(aid.table))
        tableOfCrsMatrixs().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(aid);
}

/* remove Epetra_CrsMatrix from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeCrsMatrix( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfCrsMatrixs().isType(aid->table))
        tableOfCrsMatrixs().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_CrsMatrix table */
void
CEpetra::purgeCrsMatrix(  )
{
    tableOfCrsMatrixs().purge();
}

/* store Epetra_CrsMatrix in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasCrsMatrix( const Teuchos::RCP< Epetra_CrsMatrix > & robj )
{
    return tableOfCrsMatrixs().alias(robj);
}

/* store const Epetra_CrsMatrix in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstCrsMatrix( const Teuchos::RCP< const Epetra_CrsMatrix > & robj )
{
    return tableOfCrsMatrixs().alias(robj);
}



