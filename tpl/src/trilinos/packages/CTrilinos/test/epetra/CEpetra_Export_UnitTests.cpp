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
#include "CEpetra_Comm.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Map.h"
#include "CEpetra_Object.h"
#include "CEpetra_Distributor.h"
#include "Epetra_Export.h"
#include "CEpetra_Export.h"
#include "CEpetra_Export_Cpp.hpp"
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
CT_Epetra_Export_ID_t Epetra_Export_Create ( 
  CT_Epetra_BlockMap_ID_t SourceMapID, 
  CT_Epetra_BlockMap_ID_t TargetMapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Export , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_Flex_t srcID);
  ECHO(srcID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_Flex_t tarID);
  ECHO(tarID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));

  /* Create an exporter */
  ECHO(CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(srcID.Epetra_BlockMap, tarID.Epetra_BlockMap));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_Export_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( 
  CT_Epetra_Export_ID_t ExporterID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Export , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_Flex_t srcID);
  ECHO(srcID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_Flex_t tarID);
  ECHO(tarID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));

  /* Create an exporter */
  ECHO(CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(srcID.Epetra_BlockMap, tarID.Epetra_BlockMap));

  /* Duplicate it */
  ECHO(CT_Epetra_Export_ID_t dupID = Epetra_Export_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.table, CT_Epetra_Export_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
void Epetra_Export_Destroy ( CT_Epetra_Export_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Export , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_Flex_t srcID);
  ECHO(srcID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_Flex_t tarID);
  ECHO(tarID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));

  /* Create an exporter */
  ECHO(CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(srcID.Epetra_BlockMap, tarID.Epetra_BlockMap));

  ECHO(Epetra_Export_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_Export_SourceMap ( 
  CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Export , SourceMap )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_Flex_t srcID);
  ECHO(srcID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_Flex_t tarID);
  ECHO(tarID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));

  /* Create an exporter */
  ECHO(CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(srcID.Epetra_BlockMap, tarID.Epetra_BlockMap));

  /* Get the source map using wrapper */
  ECHO(CT_Epetra_BlockMap_ID_t bsrcID2 = Epetra_Export_SourceMap(selfID));
  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(bsrcID2), NumGlobalElements);
}

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_Export_TargetMap ( 
  CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Export , TargetMap )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_Flex_t srcID);
  ECHO(srcID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_Flex_t tarID);
  ECHO(tarID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));

  /* Create an exporter */
  ECHO(CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(srcID.Epetra_BlockMap, tarID.Epetra_BlockMap));

  /* Get the source map using wrapper */
  ECHO(CT_Epetra_BlockMap_ID_t btarID2 = Epetra_Export_TargetMap(selfID));
  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(btarID2), NumGlobalElements);
}

/**********************************************************************
CT_Epetra_Distributor_ID_t Epetra_Export_Distributor ( 
  CT_Epetra_Export_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Export , Distributor )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_Flex_t srcID);
  ECHO(srcID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_Flex_t tarID);
  ECHO(tarID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));

  /* Create an exporter */
  ECHO(CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(srcID.Epetra_BlockMap, tarID.Epetra_BlockMap));

  /* No distributor will exist if not distributedglobal */
  if (Epetra_BlockMap_DistributedGlobal(srcID.Epetra_BlockMap) == TRUE) {
    ECHO(CT_Epetra_Distributor_ID_t dID = Epetra_Export_Distributor(selfID));

    /* Now check the result of the call to the wrapper function */
    TEST_EQUALITY(dID.table, CT_Epetra_Distributor_ID);
    TEST_EQUALITY_CONST(dID.index, 0);
    TEST_EQUALITY_CONST(CEpetra::getDistributor(dID).is_null(), false);
  }
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

