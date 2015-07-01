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
#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_flex_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_test_utils.hpp"

#include "CTrilinos_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, int NumVectors, boolean zeroOut );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MultiVector , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 7);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(int NumVectors = 4);
  ECHO(boolean zeroOut = FALSE);
  ECHO(CT_Epetra_MultiVector_ID_t selfID = Epetra_MultiVector_Create(MapID.Epetra_BlockMap, NumVectors, zeroOut));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_MultiVector_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( 
  CT_Epetra_MultiVector_ID_t SourceID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MultiVector , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 5);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(int NumVectors = 2);
  ECHO(boolean zeroOut = TRUE);
  ECHO(CT_Epetra_MultiVector_ID_t SourceID = Epetra_MultiVector_Create(MapID.Epetra_BlockMap, NumVectors, zeroOut));

  /* Call the wrapper function */
  ECHO(CT_Epetra_MultiVector_ID_t selfID = Epetra_MultiVector_Duplicate(SourceID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_MultiVector_ID);
  TEST_EQUALITY_CONST(selfID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, SourceID), false);
}

/**********************************************************************
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_From2DA ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, double * A, 
  int MyLDA, int NumVectors );
 **********************************************************************/

/**********************************************************************
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromAOP ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, 
  double ** ArrayOfPointers, int NumVectors );
 **********************************************************************/

/**********************************************************************
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromList ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int * Indices, int NumVectors );
 **********************************************************************/

/**********************************************************************
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromRange ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int StartIndex, int NumVectors );
 **********************************************************************/

/**********************************************************************
void Epetra_MultiVector_Destroy ( 
  CT_Epetra_MultiVector_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MultiVector , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 6);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(int NumVectors = 4);
  ECHO(boolean zeroOut = TRUE);
  ECHO(CT_Epetra_MultiVector_ID_t selfID = Epetra_MultiVector_Create(MapID.Epetra_BlockMap, NumVectors, zeroOut));

  /* Call the wrapper function */
  ECHO(Epetra_MultiVector_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
int Epetra_MultiVector_ReplaceGlobalValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, int VectorIndex, 
  double ScalarValue );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_ReplaceGlobalValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_SumIntoGlobalValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, int VectorIndex, 
  double ScalarValue );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MultiVector , SumIntoGlobalValue )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 7);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(int NumVectors = 4);
  ECHO(boolean zeroOut = TRUE);
  ECHO(CT_Epetra_MultiVector_ID_t selfID = Epetra_MultiVector_Create(MapID.Epetra_BlockMap, NumVectors, zeroOut));

  /* Initialize the multivector state and duplicate it */
  ECHO(Epetra_MultiVector_Random(selfID));
  ECHO(CT_Epetra_MultiVector_ID_t dupID = Epetra_MultiVector_Duplicate(selfID));

  /* Test the wrapper */
  ECHO(int GlobalRow = 2);
  ECHO(int VectorIndex = 1);
  ECHO(double ScalarValue = 3.8);
  ECHO(int ret = Epetra_MultiVector_SumIntoGlobalValue(selfID, GlobalRow, VectorIndex, ScalarValue));
  TEST_EQUALITY(ret, 0);
  for (int c=0; c<NumVectors; c++) {
    double *a1 = Epetra_MultiVector_getArray(selfID, c);
    double *a2 = Epetra_MultiVector_getArray(dupID, c);
    for (int r=0; r<NumGlobalElements; r++) {
      double corr = a2[r]+(r == GlobalRow ? (c == VectorIndex ? ScalarValue : 0.0) : 0.0);
      TEST_EQUALITY(corr, a1[r]);
    }
  }
}

/**********************************************************************
int Epetra_MultiVector_SumIntoGlobalValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_ReplaceMyValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyRow, int VectorIndex, 
  double ScalarValue );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_ReplaceMyValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_SumIntoMyValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyRow, int VectorIndex, 
  double ScalarValue );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_SumIntoMyValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_PutScalar ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarConstant );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_Random ( CT_Epetra_MultiVector_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_ExtractCopy_Fill2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double * A, int MyLDA );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_ExtractCopy_FillAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** ArrayOfPointers );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_ExtractView_Set2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** A, int * MyLDA );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_ExtractView_SetAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double *** ArrayOfPointers );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_Dot ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID, 
  double * Result );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_Abs ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_Reciprocal ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_Scale_Self ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarValue );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_Scale ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_Update_WithA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarThis );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_Update_WithAB ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarB, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_Norm1 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_Norm2 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_NormInf ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_NormWeighted ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t WeightsID, double * Result );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_MinValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_MaxValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_MeanValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_Multiply_Matrix ( 
  CT_Epetra_MultiVector_ID_t selfID, char TransA, char TransB, 
  double ScalarAB, CT_Epetra_MultiVector_ID_t AID, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_Multiply_ByEl ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_ReciprocalMultiply ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_SetSeed ( 
  CT_Epetra_MultiVector_ID_t selfID, unsigned int Seed_in );
 **********************************************************************/

/**********************************************************************
unsigned int Epetra_MultiVector_Seed ( 
  CT_Epetra_MultiVector_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_MultiVector_Assign ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t SourceID );
 **********************************************************************/

/**********************************************************************
double * Epetra_MultiVector_getArray ( 
  CT_Epetra_MultiVector_ID_t selfID, int i );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Vector_ID_t Epetra_MultiVector_getVector ( 
  CT_Epetra_MultiVector_ID_t selfID, int i );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_NumVectors ( 
  CT_Epetra_MultiVector_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MultiVector , NumVectors )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 7);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(int NumVectors = 4);
  ECHO(boolean zeroOut = FALSE);
  ECHO(CT_Epetra_MultiVector_ID_t selfID = Epetra_MultiVector_Create(MapID.Epetra_BlockMap, NumVectors, zeroOut));

  ECHO(int num = Epetra_MultiVector_NumVectors(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(num, NumVectors);
}

/**********************************************************************
int Epetra_MultiVector_MyLength ( 
  CT_Epetra_MultiVector_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MultiVector , MyLength )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 7);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(int NumVectors = 4);
  ECHO(boolean zeroOut = FALSE);
  ECHO(CT_Epetra_MultiVector_ID_t selfID = Epetra_MultiVector_Create(MapID.Epetra_BlockMap, NumVectors, zeroOut));

  ECHO(int len = Epetra_MultiVector_MyLength(selfID));

  /* Now check the result of the call to the wrapper function */
  ECHO(int mine = Epetra_BlockMap_NumMyElements(MapID.Epetra_BlockMap));
  TEST_EQUALITY(len, mine);
}

/**********************************************************************
int Epetra_MultiVector_GlobalLength ( 
  CT_Epetra_MultiVector_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MultiVector , GlobalLength )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 7);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(int NumVectors = 4);
  ECHO(boolean zeroOut = FALSE);
  ECHO(CT_Epetra_MultiVector_ID_t selfID = Epetra_MultiVector_Create(MapID.Epetra_BlockMap, NumVectors, zeroOut));

  ECHO(int len = Epetra_MultiVector_GlobalLength(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(len, NumGlobalElements);
}

/**********************************************************************
int Epetra_MultiVector_Stride ( CT_Epetra_MultiVector_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_MultiVector_ConstantStride ( 
  CT_Epetra_MultiVector_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_MultiVector_ReplaceMap ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID );
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

