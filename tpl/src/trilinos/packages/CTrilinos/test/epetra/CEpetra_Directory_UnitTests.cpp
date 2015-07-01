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
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Map.h"
#include "Epetra_BasicDirectory.h"
#include "Epetra_Directory.h"
#include "CEpetra_Directory.h"
#include "CEpetra_Directory_Cpp.hpp"
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
void Epetra_Directory_Destroy ( CT_Epetra_Directory_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Directory , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 9);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_Directory_ID_t selfID = Epetra_Comm_CreateDirectory(CommID, MapID.Epetra_BlockMap));

  ECHO(Epetra_Directory_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
int Epetra_Directory_GetDirectoryEntries ( 
  CT_Epetra_Directory_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID, 
  const int NumEntries, const int * GlobalEntries, int * Procs, 
  int * LocalEntries, int * EntrySizes, 
  boolean high_rank_sharing_procs );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Directory , GetDirectoryEntries )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 11);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_Directory_ID_t selfID = Epetra_Comm_CreateDirectory(CommID, MapID.Epetra_BlockMap));

  /* Now get a reference to the real object */
  ECHO(Epetra_Directory &self = *(CEpetra::getDirectory(selfID)));

  /* Now test the wrapper function */
  ECHO(const int NumEntries = 3);
  const int GlobalEntries[NumEntries] = {0, 7, 10};
  int Procs[NumEntries];
  int LocalEntries[NumEntries];
  int EntrySizes[NumEntries];
  ECHO(boolean high_rank = FALSE);
  ECHO(int ret = Epetra_Directory_GetDirectoryEntries(selfID, MapID.Epetra_BlockMap, NumEntries,
       GlobalEntries, Procs, LocalEntries, EntrySizes, high_rank));
  TEST_EQUALITY_CONST(ret, 0);

  /* Now test the original function */
  ECHO(Epetra_BlockMap &bmap = *(CEpetra::getBlockMap(MapID.Epetra_BlockMap)));
  int Procs1[NumEntries];
  int LocalEntries1[NumEntries];
  int EntrySizes1[NumEntries];
  ECHO(bool high_rank1 = false);
  ECHO(ret = self.GetDirectoryEntries(bmap, NumEntries, GlobalEntries,
       Procs1, LocalEntries1, EntrySizes1, high_rank1));
  TEST_EQUALITY_CONST(ret, 0);

  /* And compare the results */
  bool match = true;
  for (int i=0; i<NumEntries; i++) {
    if ((Procs[i] != Procs1[i]) || (LocalEntries[i] != LocalEntries1[i]) ||
        (EntrySizes[i] != EntrySizes1[i])) match = false;
  }
  TEST_EQUALITY_CONST(match, true);
}

/**********************************************************************
boolean Epetra_Directory_GIDsAllUniquelyOwned ( 
  CT_Epetra_Directory_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Directory , GIDsAllUniquelyOwned )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 5);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_Flex_t MapID);
  ECHO(MapID.Epetra_Map = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_Directory_ID_t selfID = Epetra_Comm_CreateDirectory(CommID, MapID.Epetra_BlockMap));

  /* Now get a reference to the real object */
  ECHO(Epetra_Directory &self = *(CEpetra::getDirectory(selfID)));

  /* Now try out the wrapper function */
  ECHO(boolean ret = Epetra_Directory_GIDsAllUniquelyOwned(selfID));
  ECHO(bool ret0 = (ret != FALSE ? true : false));

  /* Now try out the original function */
  ECHO(bool ret1 = self.GIDsAllUniquelyOwned());

  /* And compare the results */
  TEST_EQUALITY(ret0, ret1);
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

