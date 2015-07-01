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
#include "CEpetra_Directory.h"
#include "CEpetra_Distributor.h"
#include "CEpetra_Map.h"
#include "Epetra_SerialComm.h"
#include "CEpetra_SerialComm.h"
#include "CEpetra_SerialComm_Cpp.hpp"
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
CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Create (  );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_SerialComm_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Duplicate ( 
  CT_Epetra_SerialComm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(CT_Epetra_SerialComm_ID_t dupID = Epetra_SerialComm_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.table, CT_Epetra_SerialComm_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
CT_Epetra_Comm_ID_t Epetra_SerialComm_Clone ( 
  CT_Epetra_SerialComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Clone )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(CT_Epetra_Comm_ID_t dupID = Epetra_SerialComm_Clone(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.table, CT_Epetra_Comm_ID);
  TEST_EQUALITY_CONST(dupID.index, 0);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
void Epetra_SerialComm_Destroy ( 
  CT_Epetra_SerialComm_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(Epetra_SerialComm_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
void Epetra_SerialComm_Barrier ( CT_Epetra_SerialComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Barrier )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(Epetra_SerialComm_Barrier(selfID));
}

/**********************************************************************
int Epetra_SerialComm_Broadcast_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, int Count, 
  int Root );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Broadcast_Double )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());
  ECHO(int Root = Epetra_SerialComm_MyPID(selfID));

  ECHO(const int Count = 6);
  double MyVals[Count] = {4.6, 2.6, 3.1, 7.7, -0.5, 1.0};

  ECHO(int ret = Epetra_SerialComm_Broadcast_Double(selfID, MyVals, Count, Root));
  TEST_EQUALITY_CONST(ret, 0);
}

/**********************************************************************
int Epetra_SerialComm_Broadcast_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int Count, 
  int Root );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Broadcast_Int )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());
  ECHO(int Root = Epetra_SerialComm_MyPID(selfID));

  ECHO(const int Count = 5);
  int MyVals[Count] = {7, 2, 5, 8, 4};

  ECHO(int ret = Epetra_SerialComm_Broadcast_Int(selfID, MyVals, Count, Root));
  TEST_EQUALITY_CONST(ret, 0);
}

/**********************************************************************
int Epetra_SerialComm_Broadcast_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, int Count, 
  int Root );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Broadcast_Long )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());
  ECHO(int Root = Epetra_SerialComm_MyPID(selfID));

  ECHO(const int Count = 4);
  long MyVals[Count] = {27, 22, 25, 24};

  ECHO(int ret = Epetra_SerialComm_Broadcast_Long(selfID, MyVals, Count, Root));
  TEST_EQUALITY_CONST(ret, 0);
}

/**********************************************************************
int Epetra_SerialComm_Broadcast_Char ( 
  CT_Epetra_SerialComm_ID_t selfID, char * MyVals, int Count, 
  int Root );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Broadcast_Char )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());
  ECHO(int Root = Epetra_SerialComm_MyPID(selfID));

  ECHO(const int Count = 4);
  char MyVals[Count] = {'d', 'o', 'n', 'e'};

  ECHO(int ret = Epetra_SerialComm_Broadcast_Char(selfID, MyVals, Count, Root));
  TEST_EQUALITY_CONST(ret, 0);
}

/**********************************************************************
int Epetra_SerialComm_GatherAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , GatherAll_Double )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 6);
  double MyVals[Count] = {4.6, 2.6, 3.1, 7.7, -0.5, 1.0};
  double AllVals[Count];

  ECHO(int ret = Epetra_SerialComm_GatherAll_Double(selfID, MyVals, AllVals, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (MyVals[i] != AllVals[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_GatherAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , GatherAll_Int )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 5);
  int MyVals[Count] = {7, 2, 5, 8, 4};
  int AllVals[Count];

  ECHO(int ret = Epetra_SerialComm_GatherAll_Int(selfID, MyVals, AllVals, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (MyVals[i] != AllVals[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_GatherAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , GatherAll_Long )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 4);
  long MyVals[Count] = {27, 22, 25, 24};
  long AllVals[Count];

  ECHO(int ret = Epetra_SerialComm_GatherAll_Long(selfID, MyVals, AllVals, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (MyVals[i] != AllVals[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_SumAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , SumAll_Double )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 6);
  double PartialSums[Count] = {4.6, 2.6, 3.1, 7.7, -0.5, 1.0};
  double GlobalSums[Count];

  ECHO(int ret = Epetra_SerialComm_SumAll_Double(selfID, PartialSums, GlobalSums, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (PartialSums[i] != GlobalSums[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_SumAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialSums, 
  int * GlobalSums, int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , SumAll_Int )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 5);
  int PartialSums[Count] = {7, 2, 5, 8, 4};
  int GlobalSums[Count];

  ECHO(int ret = Epetra_SerialComm_SumAll_Int(selfID, PartialSums, GlobalSums, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (PartialSums[i] != GlobalSums[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_SumAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , SumAll_Long )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 4);
  long PartialSums[Count] = {27, 22, 25, 24};
  long GlobalSums[Count];

  ECHO(int ret = Epetra_SerialComm_SumAll_Long(selfID, PartialSums, GlobalSums, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (PartialSums[i] != GlobalSums[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_MaxAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , MaxAll_Double )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 6);
  double PartialMaxs[Count] = {4.6, 2.6, 3.1, 7.7, -0.5, 1.0};
  double GlobalMaxs[Count];

  ECHO(int ret = Epetra_SerialComm_MaxAll_Double(selfID, PartialMaxs, GlobalMaxs, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (PartialMaxs[i] != GlobalMaxs[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_MaxAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMaxs, 
  int * GlobalMaxs, int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , MaxAll_Int )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 5);
  int PartialMaxs[Count] = {7, 2, 5, 8, 4};
  int GlobalMaxs[Count];

  ECHO(int ret = Epetra_SerialComm_MaxAll_Int(selfID, PartialMaxs, GlobalMaxs, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (PartialMaxs[i] != GlobalMaxs[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_MaxAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , MaxAll_Long )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 4);
  long PartialMaxs[Count] = {27, 22, 25, 24};
  long GlobalMaxs[Count];

  ECHO(int ret = Epetra_SerialComm_MaxAll_Long(selfID, PartialMaxs, GlobalMaxs, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (PartialMaxs[i] != GlobalMaxs[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_MinAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , MinAll_Double )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 6);
  double PartialMins[Count] = {4.6, 2.6, 3.1, 7.7, -0.5, 1.0};
  double GlobalMins[Count];

  ECHO(int ret = Epetra_SerialComm_MinAll_Double(selfID, PartialMins, GlobalMins, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (PartialMins[i] != GlobalMins[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_MinAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMins, 
  int * GlobalMins, int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , MinAll_Int )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 5);
  int PartialMins[Count] = {7, 2, 5, 8, 4};
  int GlobalMins[Count];

  ECHO(int ret = Epetra_SerialComm_MinAll_Int(selfID, PartialMins, GlobalMins, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (PartialMins[i] != GlobalMins[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_MinAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , MinAll_Long )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 4);
  long PartialMins[Count] = {27, 22, 25, 24};
  long GlobalMins[Count];

  ECHO(int ret = Epetra_SerialComm_MinAll_Long(selfID, PartialMins, GlobalMins, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (PartialMins[i] != GlobalMins[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_ScanSum_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  double * ScanSums, int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , ScanSum_Double )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 6);
  double MyVals[Count] = {4.6, 2.6, 3.1, 7.7, -0.5, 1.0};
  double ScanSums[Count];

  ECHO(int ret = Epetra_SerialComm_ScanSum_Double(selfID, MyVals, ScanSums, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (MyVals[i] != ScanSums[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_ScanSum_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , ScanSum_Int )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 5);
  int MyVals[Count] = {7, 2, 5, 8, 4};
  int ScanSums[Count];

  ECHO(int ret = Epetra_SerialComm_ScanSum_Int(selfID, MyVals, ScanSums, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (MyVals[i] != ScanSums[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_ScanSum_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , ScanSum_Long )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(const int Count = 4);
  long MyVals[Count] = {27, 22, 25, 24};
  long ScanSums[Count];

  ECHO(int ret = Epetra_SerialComm_ScanSum_Long(selfID, MyVals, ScanSums, Count));
  TEST_EQUALITY_CONST(ret, 0);

  /* Compare the two vectors (no change when serial) */
  bool match = true;
  for (int i=0; i<Count; i++) {
    if (MyVals[i] != ScanSums[i]) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
int Epetra_SerialComm_MyPID ( CT_Epetra_SerialComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , MyPID )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(int ret = Epetra_SerialComm_MyPID(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY_CONST(ret, 0);
}

/**********************************************************************
int Epetra_SerialComm_NumProc ( CT_Epetra_SerialComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , NumProc )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(int ret = Epetra_SerialComm_NumProc(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY_CONST(ret, 1);
}

/**********************************************************************
CT_Epetra_Distributor_ID_t Epetra_SerialComm_CreateDistributor ( 
  CT_Epetra_SerialComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , CreateDistributor )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(CT_Epetra_Distributor_ID_t disID = Epetra_SerialComm_CreateDistributor(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(disID.table, CT_Epetra_Distributor_ID);
  TEST_EQUALITY_CONST(disID.index, 0);
}

/**********************************************************************
CT_Epetra_Directory_ID_t Epetra_SerialComm_CreateDirectory ( 
  CT_Epetra_SerialComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , CreateDirectory )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_Flex_t CommID);
  ECHO(CommID.Epetra_SerialComm = Epetra_SerialComm_Create());

  ECHO(int NumGlobalElements = 9);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID.Epetra_Comm));
  ECHO(CT_Epetra_Directory_ID_t dirID = Epetra_SerialComm_CreateDirectory(
      CommID.Epetra_SerialComm, MapID.Epetra_BlockMap));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dirID.table, CT_Epetra_Directory_ID);
  TEST_EQUALITY_CONST(dirID.index, 0);
}

/**********************************************************************
void Epetra_SerialComm_Assign ( 
  CT_Epetra_SerialComm_ID_t selfID, CT_Epetra_SerialComm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_SerialComm , Assign )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_SerialComm_ID_t commID = Epetra_SerialComm_Create());
  ECHO(CT_Epetra_SerialComm_ID_t selfID = Epetra_SerialComm_Create());

  ECHO(Epetra_SerialComm_Assign(selfID, commID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, commID), false);
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

