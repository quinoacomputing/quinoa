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
#include "CEpetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "CEpetra_Vector.h"
#include "CEpetra_Vector_Cpp.hpp"
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
CT_Epetra_Vector_ID_t Epetra_Vector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 7);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(boolean zeroOut = FALSE);
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Create(MapID.Epetra_BlockMap, zeroOut));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_Vector_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_Vector_ID_t Epetra_Vector_Duplicate ( 
  CT_Epetra_Vector_ID_t SourceID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 8);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(boolean zeroOut = TRUE);
  ECHO(CT_Epetra_Vector_ID_t SourceID = Epetra_Vector_Create(MapID.Epetra_BlockMap, zeroOut));

  /* Call the wrapper function */
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Duplicate(SourceID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_Vector_ID);
  TEST_EQUALITY_CONST(selfID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, SourceID), false);
}

/**********************************************************************
CT_Epetra_Vector_ID_t Epetra_Vector_Create_FromArray ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, double * V );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , Create_FromArray )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(const int NumGlobalElements = 6);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  double V[NumGlobalElements] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Create_FromArray(
       CT_Epetra_DataAccess_E_Copy, MapID.Epetra_BlockMap, V));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_Vector_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_Vector_ID_t Epetra_Vector_FromSource ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int Index );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , FromSource )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 9);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(int NumVectors = 3);
  ECHO(boolean zeroOut = TRUE);
  ECHO(CT_Epetra_MultiVector_ID_t SourceID = Epetra_MultiVector_Create(MapID.Epetra_BlockMap, NumVectors, zeroOut));
  ECHO(int Index = 1);
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_FromSource(CT_Epetra_DataAccess_E_View, SourceID, Index));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_Vector_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
void Epetra_Vector_Destroy ( CT_Epetra_Vector_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 13);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(boolean zeroOut = TRUE);
  ECHO(CT_Epetra_Vector_ID_t selfID = Epetra_Vector_Create(MapID.Epetra_BlockMap, zeroOut));

  /* Call the wrapper function */
  ECHO(Epetra_Vector_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
int Epetra_Vector_ReplaceGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , ReplaceGlobalValues )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 11);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  /* Create a vector using CTrilinos and duplicate it outside */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID.Epetra_BlockMap, TRUE));
  ECHO(Teuchos::RCP<Epetra_Vector> rcpV = CEpetra::getVector(vID));
  ECHO(Epetra_Vector *v2 = new Epetra_Vector(*rcpV));

  /* Set up the problem */
  ECHO(const int NumEntries = 5);
  double vals[NumEntries] = {2.0, 7.4, 21.0, 8.5, 6.7};
  int inds[NumEntries] = {0, 1, 4, 8, 10};

  /* Try out the CTrilinos interface */
  ECHO(Epetra_Vector_ReplaceGlobalValues(vID, NumEntries, vals, inds));

  /* Do the same thing to the control vector */
  ECHO(v2->ReplaceGlobalValues(NumEntries, vals, inds));

  /* Figure out how many elements on this processor */
  ECHO(int NumMyElements = v2->MyLength());

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumMyElements; i++) {
    if ((*rcpV)[i] != (*v2)[i]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_Vector_ReplaceMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , ReplaceMyValues )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumMyElements = 5);
  ECHO(int NumGlobalElements = -1);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create_Linear(
       NumGlobalElements, NumMyElements, IndexBase, CommID));

  /* Create a vector using CTrilinos and duplicate it outside */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID.Epetra_BlockMap, TRUE));
  ECHO(Teuchos::RCP<Epetra_Vector> rcpV = CEpetra::getVector(vID));
  ECHO(Epetra_Vector *v2 = new Epetra_Vector(*rcpV));

  /* Set up the problem */
  ECHO(const int NumEntries = 3);
  double vals[NumEntries] = {1.0, 5.2, 2.1};
  int inds[NumEntries] = {0, 1, 4};

  /* Try out the CTrilinos interface */
  ECHO(Epetra_Vector_ReplaceMyValues(vID, NumEntries, vals, inds));

  /* Do the same thing to the control vector */
  ECHO(v2->ReplaceMyValues(NumEntries, vals, inds));

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumMyElements; i++) {
    if ((*rcpV)[i] != (*v2)[i]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_Vector_SumIntoGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , SumIntoGlobalValues )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 11);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  /* Create a vector using CTrilinos and initialize it as desired */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID.Epetra_BlockMap, TRUE));
  ECHO(Teuchos::RCP<Epetra_Vector> rcpV = CEpetra::getVector(vID));
  ECHO(const int NumEntries0 = 5);
  double vals0[NumEntries0] = {2.0, 7.4, 21.0, 8.5, 6.7};
  int inds0[NumEntries0] = {0, 1, 4, 8, 9};
  ECHO(rcpV->ReplaceGlobalValues(NumEntries0, vals0, inds0));

  /* Duplicate CTrilinos vector outside */
  ECHO(Epetra_Vector *v2 = new Epetra_Vector(*rcpV));

  /* Set up the problem */
  ECHO(const int NumEntries = 4);
  double vals[NumEntries] = {5.1, -2.0, 3.6, 1.1};
  int inds[NumEntries] = {0, 2, 4, 10};

  /* Try out the CTrilinos interface */
  ECHO(Epetra_Vector_SumIntoGlobalValues(vID, NumEntries, vals, inds));

  /* Do the same thing to the control vector */
  ECHO(v2->SumIntoGlobalValues(NumEntries, vals, inds));

  /* Figure out how many elements on this processor */
  ECHO(int NumMyElements = v2->MyLength());

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumMyElements; i++) {
    if ((*rcpV)[i] != (*v2)[i]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_Vector_SumIntoMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , SumIntoMyValues )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumMyElements = 5);
  ECHO(int NumGlobalElements = -1);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create_Linear(
       NumGlobalElements, NumMyElements, IndexBase, CommID));

  /* Create a vector using CTrilinos and initialize it as desired */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID.Epetra_BlockMap, TRUE));
  ECHO(Teuchos::RCP<Epetra_Vector> rcpV = CEpetra::getVector(vID));
  ECHO(const int NumEntries0 = 2);
  double vals0[NumEntries0] = {1.5, 6.9};
  int inds0[NumEntries0] = {1, 4};
  ECHO(rcpV->ReplaceGlobalValues(NumEntries0, vals0, inds0));

  /* Duplicate CTrilinos vector outside */
  ECHO(Epetra_Vector *v2 = new Epetra_Vector(*rcpV));

  /* Set up the problem */
  ECHO(const int NumEntries = 2);
  double vals[NumEntries] = {-2.0, 5.1};
  int inds[NumEntries] = {0, 1};

  /* Try out the CTrilinos interface */
  ECHO(Epetra_Vector_SumIntoMyValues(vID, NumEntries, vals, inds));

  /* Do the same thing to the control vector */
  ECHO(v2->SumIntoMyValues(NumEntries, vals, inds));

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumMyElements; i++) {
    if ((*rcpV)[i] != (*v2)[i]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_Vector_ReplaceGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , ReplaceGlobalValues_BlockPos )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 4);
  ECHO(int ElementSize = 3);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));

  /* Create a vector using CTrilinos and duplicate it outside */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID, TRUE));
  ECHO(Teuchos::RCP<Epetra_Vector> rcpV = CEpetra::getVector(vID));
  ECHO(Epetra_Vector *v2 = new Epetra_Vector(*rcpV));

  /* Set up the problem */
  ECHO(const int NumEntries = 3);
  ECHO(int BlockOffset = 1);
  double vals[NumEntries] = {7.4, 21.0, 8.5};
  int inds[NumEntries] = {0, 1, 3};

  /* Try out the CTrilinos interface */
  ECHO(int ret = Epetra_Vector_ReplaceGlobalValues_BlockPos(
       vID, NumEntries, BlockOffset, vals, inds));
  TEST_EQUALITY_CONST(ret, 0);

  /* Do the same thing to the control vector */
  ECHO(ret = v2->ReplaceGlobalValues(NumEntries, BlockOffset, vals, inds));
  TEST_EQUALITY_CONST(ret, 0);

  /* Figure out how many elements on this processor */
  ECHO(int NumMyElements = v2->MyLength());

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumMyElements; i++) {
    if ((*rcpV)[i] != (*v2)[i]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_Vector_ReplaceMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , ReplaceMyValues_BlockPos )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumMyElements = 2);
  ECHO(int ElementSize = 3);
  ECHO(int NumGlobalElements = -1);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Create_Linear(
       NumGlobalElements, NumMyElements, ElementSize, IndexBase, CommID));

  /* Create a vector using CTrilinos and duplicate it outside */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID, TRUE));
  ECHO(Teuchos::RCP<Epetra_Vector> rcpV = CEpetra::getVector(vID));
  ECHO(Epetra_Vector *v2 = new Epetra_Vector(*rcpV));

  /* Set up the problem */
  ECHO(const int NumEntries = 1);
  ECHO(int BlockOffset = 2);
  double vals[NumEntries] = {5.2};
  int inds[NumEntries] = {0};

  /* Try out the CTrilinos interface */
  ECHO(int ret = Epetra_Vector_ReplaceMyValues_BlockPos(
       vID, NumEntries, BlockOffset, vals, inds));
  TEST_EQUALITY_CONST(ret, 0);

  /* Do the same thing to the control vector */
  ECHO(ret = v2->ReplaceMyValues(NumEntries, BlockOffset, vals, inds));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumMyElements; i++) {
    if ((*rcpV)[i] != (*v2)[i]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_Vector_SumIntoGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , SumIntoGlobalValues_BlockPos )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 5);
  ECHO(int ElementSize = 2);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));

  /* Create a vector using CTrilinos and initialize it as desired */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID, TRUE));
  ECHO(Teuchos::RCP<Epetra_Vector> rcpV = CEpetra::getVector(vID));

  ECHO(const int NumEntries0 = 3);
  ECHO(int BlockOffset0 = 0);
  double vals0[NumEntries0] = {21.0, 8.5, 6.7};
  int inds0[NumEntries0] = {0, 3, 4};
  ECHO(int ret = rcpV->ReplaceGlobalValues(NumEntries0, BlockOffset0, vals0, inds0));
  TEST_EQUALITY_CONST(ret, 0);

  ECHO(const int NumEntries1 = 2);
  ECHO(int BlockOffset1 = 1);
  double vals1[NumEntries1] = {7.2, 0.1};
  int inds1[NumEntries1] = {2, 4};
  ECHO(ret = rcpV->ReplaceGlobalValues(NumEntries1, BlockOffset1, vals1, inds1));
  TEST_EQUALITY_CONST(ret, 0);

  /* Duplicate CTrilinos vector outside */
  ECHO(Epetra_Vector *v2 = new Epetra_Vector(*rcpV));

  /* Set up the problem */
  ECHO(const int NumEntries = 4);
  ECHO(int BlockOffset = 0);
  double vals[NumEntries] = {5.1, -2.0, 3.6, 1.1};
  int inds[NumEntries] = {0, 2, 3, 4};

  /* Try out the CTrilinos interface */
  ECHO(ret = Epetra_Vector_SumIntoGlobalValues_BlockPos(
       vID, NumEntries, BlockOffset, vals, inds));
  TEST_EQUALITY_CONST(ret, 0);

  /* Do the same thing to the control vector */
  ECHO(ret = v2->SumIntoGlobalValues(NumEntries, BlockOffset, vals, inds));
  TEST_EQUALITY_CONST(ret, 0);

  /* Figure out how many elements on this processor */
  ECHO(int NumMyElements = v2->MyLength());

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumMyElements; i++) {
    if ((*rcpV)[i] != (*v2)[i]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_Vector_SumIntoMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , SumIntoMyValues_BlockPos )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumMyElements = 5);
  ECHO(int NumGlobalElements = -1);
  ECHO(int ElementSize = 2);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_BlockMap_ID_t MapID = Epetra_BlockMap_Create_Linear(
       NumGlobalElements, NumMyElements, ElementSize, IndexBase, CommID));

  /* Create a vector using CTrilinos and initialize it as desired */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID, TRUE));
  ECHO(Teuchos::RCP<Epetra_Vector> rcpV = CEpetra::getVector(vID));

  ECHO(const int NumEntries0 = 2);
  ECHO(int BlockOffset0 = 0);
  double vals0[NumEntries0] = {1.5, 6.9};
  int inds0[NumEntries0] = {1, 4};
  ECHO(int ret = rcpV->ReplaceGlobalValues(NumEntries0, BlockOffset0, vals0, inds0));
  TEST_EQUALITY_CONST(ret, 0);

  ECHO(const int NumEntries1 = 3);
  ECHO(int BlockOffset1 = 1);
  double vals1[NumEntries1] = {3.9, 23.2, 7.3};
  int inds1[NumEntries1] = {0, 2, 4};
  ECHO(ret = rcpV->ReplaceGlobalValues(NumEntries1, BlockOffset1, vals1, inds1));
  TEST_EQUALITY_CONST(ret, 0);

  /* Duplicate CTrilinos vector outside */
  ECHO(Epetra_Vector *v2 = new Epetra_Vector(*rcpV));

  /* Set up the problem */
  ECHO(const int NumEntries = 2);
  ECHO(int BlockOffset = 1);
  double vals[NumEntries] = {-2.0, 5.1};
  int inds[NumEntries] = {0, 1};

  /* Try out the CTrilinos interface */
  ECHO(ret = Epetra_Vector_SumIntoMyValues_BlockPos(vID, NumEntries, BlockOffset, vals, inds));
  TEST_EQUALITY_CONST(ret, 0);

  /* Do the same thing to the control vector */
  ECHO(ret = v2->SumIntoMyValues(NumEntries, BlockOffset, vals, inds));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumMyElements; i++) {
    if ((*rcpV)[i] != (*v2)[i]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_Vector_ExtractCopy ( 
  CT_Epetra_Vector_ID_t selfID, double * V );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , ExtractCopy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 7);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  /* Create a vector using CTrilinos */
  ECHO(CT_Epetra_Vector_ID_Flex_t vID);
  ECHO(vID.Epetra_Vector = Epetra_Vector_Create(MapID.Epetra_BlockMap, TRUE));
  ECHO(const int NumEntries = 4);
  double vals[NumEntries] = {7.4, 12.0, 6.7, 21.3};
  int inds[NumEntries] = {0, 2, 3, 6};
  ECHO(Epetra_Vector_ReplaceGlobalValues(vID.Epetra_Vector, NumEntries, vals, inds));

  /* Figure out how many elements on this processor */
  ECHO(int NumMyElements = Epetra_MultiVector_MyLength(vID.Epetra_MultiVector));

  /* Try out the CTrilinos interface */
  ECHO(double *V = (double *)malloc(NumMyElements*sizeof(double)));
  ECHO(int ret = Epetra_Vector_ExtractCopy(vID.Epetra_Vector, V));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumEntries; i++) {
    if (vals[i] != V[inds[i]]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_Vector_ExtractView ( 
  CT_Epetra_Vector_ID_t selfID, double ** V );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , ExtractView )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 7);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  /* Create a vector using CTrilinos */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID.Epetra_BlockMap, TRUE));
  ECHO(const int NumEntries = 4);
  double vals[NumEntries] = {7.4, 12.0, 6.7, 21.3};
  int inds[NumEntries] = {0, 2, 3, 6};
  ECHO(Epetra_Vector_ReplaceGlobalValues(vID, NumEntries, vals, inds));

  /* Try out the CTrilinos interface */
  ECHO(double *V = NULL);
  ECHO(int ret = Epetra_Vector_ExtractView(vID, &V));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumEntries; i++) {
    if (vals[i] != V[inds[i]]) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
double Epetra_Vector_getElement ( 
  CT_Epetra_Vector_ID_t selfID, int index );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Vector , getElement )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 7);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  /* Create a vector using CTrilinos */
  ECHO(CT_Epetra_Vector_ID_t vID = Epetra_Vector_Create(MapID.Epetra_BlockMap, TRUE));
  ECHO(const int NumEntries = 4);
  double vals[NumEntries] = {7.4, 12.0, 6.7, 21.3};
  int inds[NumEntries] = {0, 2, 3, 6};
  ECHO(Epetra_Vector_ReplaceGlobalValues(vID, NumEntries, vals, inds));

  /* Compare the two vectors */
  bool match = true;
  for (int i=0; i<NumEntries; i++) {
    if (vals[i] != Epetra_Vector_getElement(vID, inds[i])) match = false;
  }

  TEST_EQUALITY_CONST(match, true);
}


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

