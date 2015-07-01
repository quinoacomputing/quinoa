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


#ifdef HAVE_MPI

#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Directory.h"
#include "CEpetra_Distributor.h"
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "CEpetra_MpiComm.h"
#include "CEpetra_MpiComm_Cpp.hpp"
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
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create ( MPI_Comm comm );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_MpiComm_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( 
  CT_Epetra_MpiComm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(CT_Epetra_MpiComm_ID_t dupID = Epetra_MpiComm_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.table, CT_Epetra_MpiComm_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
CT_Epetra_Comm_ID_t Epetra_MpiComm_Clone ( 
  CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Clone )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(CT_Epetra_Comm_ID_t dupID = Epetra_MpiComm_Clone(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.table, CT_Epetra_Comm_ID);
  TEST_EQUALITY_CONST(dupID.index, 0);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
void Epetra_MpiComm_Destroy ( CT_Epetra_MpiComm_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(Epetra_MpiComm_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
void Epetra_MpiComm_Barrier ( CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Barrier )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(Epetra_MpiComm_Barrier(selfID));

  /* ??? */
}

/**********************************************************************
int Epetra_MpiComm_Broadcast_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, int Count, 
  int Root );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_Broadcast_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int Count, int Root );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_Broadcast_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, int Count, int Root );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_Broadcast_Char ( 
  CT_Epetra_MpiComm_ID_t selfID, char * MyVals, int Count, int Root );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_GatherAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_GatherAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_GatherAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_SumAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_SumAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialSums, int * GlobalSums, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_SumAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_MaxAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_MaxAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMaxs, int * GlobalMaxs, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_MaxAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_MinAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_MinAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMins, int * GlobalMins, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_MinAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_ScanSum_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * ScanSums, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_ScanSum_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count );
 **********************************************************************/

/**********************************************************************
int Epetra_MpiComm_ScanSum_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count );
 **********************************************************************/

/**********************************************************************
MPI_Comm Epetra_MpiComm_Comm ( CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Comm )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(MPI_Comm c = Epetra_MpiComm_Comm(selfID));

  /* Now check the result of the call to the wrapper function */
  /* ??? */

  (void) c;
}

/**********************************************************************
int Epetra_MpiComm_MyPID ( CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , MyPID )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(int ret = Epetra_MpiComm_MyPID(selfID));

  /* Now check the result of the call to the wrapper function */
  /* ??? */

  (void) ret;
}

/**********************************************************************
int Epetra_MpiComm_NumProc ( CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , NumProc )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(int ret = Epetra_MpiComm_NumProc(selfID));

  /* Now check the result of the call to the wrapper function */
  /* ??? */

  (void) ret;
}

/**********************************************************************
CT_Epetra_Distributor_ID_t Epetra_MpiComm_CreateDistributor ( 
  CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , CreateDistributor )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(CT_Epetra_Distributor_ID_t disID = Epetra_MpiComm_CreateDistributor(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(disID.table, CT_Epetra_Distributor_ID);
  TEST_EQUALITY_CONST(disID.index, 0);
}

/**********************************************************************
CT_Epetra_Directory_ID_t Epetra_MpiComm_CreateDirectory ( 
  CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , CreateDirectory )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_Flex_t CommID);
  ECHO(CommID.Epetra_MpiComm = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(int NumGlobalElements = 9);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID.Epetra_Comm));
  ECHO(CT_Epetra_Directory_ID_t dirID = Epetra_MpiComm_CreateDirectory(CommID.Epetra_MpiComm, MapID.Epetra_BlockMap));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dirID.table, CT_Epetra_Directory_ID);
  TEST_EQUALITY_CONST(dirID.index, 0);
}

/**********************************************************************
int Epetra_MpiComm_GetMpiTag ( CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

/**********************************************************************
MPI_Comm Epetra_MpiComm_GetMpiComm ( CT_Epetra_MpiComm_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_MpiComm_Assign ( 
  CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_MpiComm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_MpiComm , Assign )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_MpiComm_ID_t commID = Epetra_MpiComm_Create(MPI_COMM_WORLD));
  ECHO(CT_Epetra_MpiComm_ID_t selfID = Epetra_MpiComm_Create(MPI_COMM_WORLD));

  ECHO(Epetra_MpiComm_Assign(selfID, commID));

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

#endif /* HAVE_MPI */

