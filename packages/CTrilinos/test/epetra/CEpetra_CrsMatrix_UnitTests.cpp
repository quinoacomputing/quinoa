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
#include "CEpetra_BlockMap.h"
#include "CEpetra_Comm.h"
#include "CEpetra_CrsGraph.h"
#include "CEpetra_Export.h"
#include "CEpetra_Import.h"
#include "CEpetra_Map.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "CEpetra_CrsMatrix.h"
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_test_utils.hpp"

#include "CTrilinos_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  const int * NumEntriesPerRow, boolean StaticProfile );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_CrsMatrix , Create_VarPerRow )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(const int NumGlobalElements = 4);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t MapID = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  int NumIndicesPerRow[NumGlobalElements] = {3, 2, 6, 4};
  ECHO(CT_Epetra_DataAccess_E_t CV = CT_Epetra_DataAccess_E_Copy);
  ECHO(CT_Epetra_CrsMatrix_ID_t selfID = Epetra_CrsMatrix_Create_VarPerRow(
       CV, MapID, NumIndicesPerRow, FALSE));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_CrsMatrix_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  int NumEntriesPerRow, boolean StaticProfile );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_CrsMatrix , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 5);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t MapID = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  ECHO(int NumIndicesPerRow = 7);
  ECHO(CT_Epetra_DataAccess_E_t CV = CT_Epetra_DataAccess_E_Copy);
  ECHO(CT_Epetra_CrsMatrix_ID_t selfID = Epetra_CrsMatrix_Create(
       CV, MapID, NumIndicesPerRow, FALSE));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_CrsMatrix_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, const int * NumEntriesPerRow, 
  boolean StaticProfile );
 **********************************************************************/

/**********************************************************************
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, 
  boolean StaticProfile );
 **********************************************************************/

/**********************************************************************
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_FromGraph ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_CrsGraph_ID_t GraphID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Duplicate ( 
  CT_Epetra_CrsMatrix_ID_t MatrixID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_CrsMatrix , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 4);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t MapID = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  /* Create the source matrix */
  ECHO(int NumIndicesPerRow = 4);
  ECHO(CT_Epetra_DataAccess_E_t CV = CT_Epetra_DataAccess_E_Copy);
  ECHO(CT_Epetra_CrsMatrix_ID_t selfID = Epetra_CrsMatrix_Create(
       CV, MapID, NumIndicesPerRow, FALSE));

  /* Initialize the source matrix */
  ECHO(double val = 1.0);
  ECHO(int ret = Epetra_CrsMatrix_PutScalar(selfID, val));
  TEST_EQUALITY(ret, 0);
  ECHO(ret = Epetra_CrsMatrix_FillComplete(selfID, TRUE));
  TEST_EQUALITY(ret, 0);

  /* Call wrapper to duplicate the matrix */
  ECHO(CT_Epetra_CrsMatrix_ID_t dupID = Epetra_CrsMatrix_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.table, CT_Epetra_CrsMatrix_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
void Epetra_CrsMatrix_Destroy ( CT_Epetra_CrsMatrix_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_CrsMatrix , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 6);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t MapID = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  ECHO(int NumIndicesPerRow = 3);
  ECHO(CT_Epetra_DataAccess_E_t CV = CT_Epetra_DataAccess_E_Copy);
  ECHO(CT_Epetra_CrsMatrix_ID_t selfID = Epetra_CrsMatrix_Create(
       CV, MapID, NumIndicesPerRow, FALSE));

  ECHO(Epetra_CrsMatrix_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
void Epetra_CrsMatrix_Assign ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_CrsMatrix_ID_t srcID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_PutScalar ( 
  CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_Scale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_InsertGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ReplaceGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_SumIntoGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_InsertMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ReplaceMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_SumIntoMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ReplaceDiagonalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_FillComplete ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean OptimizeDataStorage );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_FillComplete_UsingMaps ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Map_ID_t DomainMapID, 
  CT_Epetra_Map_ID_t RangeMapID, boolean OptimizeDataStorage );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_OptimizeStorage ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_MakeDataContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ExtractMyRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ExtractMyRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ExtractGlobalRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values, int ** Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ExtractMyRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values, int ** Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ExtractGlobalRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ExtractMyRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_Multiply_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_Multiply1_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_Multiply_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_Multiply1_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_Solve_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_Vector_ID_t xID, 
  CT_Epetra_Vector_ID_t yID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_Solve_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_InvRowSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_InvRowMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_LeftScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_InvColSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_InvColMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_RightScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_Filled ( CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_StorageOptimized ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_IndicesAreGlobal ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_IndicesAreLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_IndicesAreContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_LowerTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_UpperTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_NoDiagonal ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
double Epetra_CrsMatrix_NormInf ( CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
double Epetra_CrsMatrix_NormOne ( CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
double Epetra_CrsMatrix_NormFrobenius ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumGlobalNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumGlobalRows ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumGlobalCols ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumGlobalDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumMyNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumMyRows ( CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumMyCols ( CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumMyDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumAllocatedGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_MaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_GlobalMaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumAllocatedMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_IndexBase ( CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_StaticGraph ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_CrsGraph_ID_t Epetra_CrsMatrix_Graph ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ReplaceRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_HaveColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ReplaceColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_CrsMatrix_ColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_CrsMatrix_DomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_CrsMatrix_RangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Import_ID_t Epetra_CrsMatrix_Importer ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Export_ID_t Epetra_CrsMatrix_Exporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Comm_ID_t Epetra_CrsMatrix_Comm ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_LRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_GRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_LCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_GCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_MyGRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_MyLRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_MyGCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_MyLCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_MyGlobalRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GID );
 **********************************************************************/

/**********************************************************************
const char * Epetra_CrsMatrix_Label ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_SetUseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean UseTranspose_in );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_Apply ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_ApplyInverse ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_HasNormInf ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_CrsMatrix_UseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorDomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorRangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_NumMyRowEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Import_ID_t Epetra_CrsMatrix_RowMatrixImporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
double * Epetra_CrsMatrix_getRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Loc );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_CrsMatrix_ImportMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_TransformToLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_CrsMatrix_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Map_ID_t DomainMapID, 
  CT_Epetra_Map_ID_t RangeMapID );
 **********************************************************************/

/**********************************************************************/

//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( T ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  DEBUG_UNIT_TEST_GROUP( T )


} // namespace

