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
#include "CEpetra_Map.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_CrsMatrix.h"
#include "Epetra_Operator.h"
#include "CEpetra_Operator.h"
#include "CEpetra_Operator_Cpp.hpp"
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
void Epetra_Operator_Destroy ( CT_Epetra_Operator_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Operator , Destroy )
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
  ECHO(CT_Epetra_CrsMatrix_ID_Flex_t crsID);
  ECHO(crsID.Epetra_CrsMatrix = Epetra_CrsMatrix_Create(
       CV, MapID, NumIndicesPerRow, FALSE));

  /* Initialize the source matrix */
  ECHO(double val = 1.0);
  ECHO(int ret = Epetra_CrsMatrix_PutScalar(crsID.Epetra_CrsMatrix, val));
  TEST_EQUALITY(ret, 0);
  ECHO(ret = Epetra_CrsMatrix_FillComplete(crsID.Epetra_CrsMatrix, TRUE));
  TEST_EQUALITY(ret, 0);

  ECHO(Epetra_Operator_Destroy(&crsID.Epetra_Operator));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(crsID.Epetra_Operator.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(crsID.Epetra_Operator.index, -1);
}

/**********************************************************************
int Epetra_Operator_SetUseTranspose ( 
  CT_Epetra_Operator_ID_t selfID, boolean UseTranspose );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Operator , SetUseTranspose )
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
  ECHO(CT_Epetra_CrsMatrix_ID_Flex_t crsID);
  ECHO(crsID.Epetra_CrsMatrix = Epetra_CrsMatrix_Create(
       CV, MapID, NumIndicesPerRow, FALSE));

  /* Initialize the source matrix */
  ECHO(double val = 1.0);
  ECHO(int ret = Epetra_CrsMatrix_PutScalar(crsID.Epetra_CrsMatrix, val));
  TEST_EQUALITY(ret, 0);
  ECHO(ret = Epetra_CrsMatrix_FillComplete(crsID.Epetra_CrsMatrix, TRUE));
  TEST_EQUALITY(ret, 0);

  /* Test true */
  ECHO(boolean tr = TRUE);
  ECHO(ret = Epetra_Operator_SetUseTranspose(crsID.Epetra_Operator, tr));
  TEST_EQUALITY_CONST(ret, 0);
  ECHO(boolean tr2 = Epetra_Operator_UseTranspose(crsID.Epetra_Operator));
  TEST_EQUALITY(tr, tr2);

  /* Test false */
  ECHO(tr = FALSE);
  ECHO(ret = Epetra_Operator_SetUseTranspose(crsID.Epetra_Operator, tr));
  TEST_EQUALITY_CONST(ret, 0);
  ECHO(tr2 = Epetra_Operator_UseTranspose(crsID.Epetra_Operator));
  TEST_EQUALITY(tr, tr2);
}

/**********************************************************************
int Epetra_Operator_Apply ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
int Epetra_Operator_ApplyInverse ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
double Epetra_Operator_NormInf ( CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
const char * Epetra_Operator_Label ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_Operator_UseTranspose ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Operator , UseTranspose )
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
  ECHO(CT_Epetra_CrsMatrix_ID_Flex_t crsID);
  ECHO(crsID.Epetra_CrsMatrix = Epetra_CrsMatrix_Create(
       CV, MapID, NumIndicesPerRow, FALSE));

  /* Initialize the source matrix */
  ECHO(double val = 1.0);
  ECHO(int ret = Epetra_CrsMatrix_PutScalar(crsID.Epetra_CrsMatrix, val));
  TEST_EQUALITY(ret, 0);
  ECHO(ret = Epetra_CrsMatrix_FillComplete(crsID.Epetra_CrsMatrix, TRUE));
  TEST_EQUALITY(ret, 0);

  /* Test true */
  ECHO(boolean tr = TRUE);
  ECHO(ret = Epetra_Operator_SetUseTranspose(crsID.Epetra_Operator, tr));
  TEST_EQUALITY_CONST(ret, 0);
  ECHO(boolean tr2 = Epetra_Operator_UseTranspose(crsID.Epetra_Operator));
  TEST_EQUALITY(tr, tr2);

  /* Test false */
  ECHO(tr = FALSE);
  ECHO(ret = Epetra_Operator_SetUseTranspose(crsID.Epetra_Operator, tr));
  TEST_EQUALITY_CONST(ret, 0);
  ECHO(tr2 = Epetra_Operator_UseTranspose(crsID.Epetra_Operator));
  TEST_EQUALITY(tr, tr2);
}

/**********************************************************************
boolean Epetra_Operator_HasNormInf ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Comm_ID_t Epetra_Operator_Comm ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Operator_OperatorDomainMap ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Operator_OperatorRangeMap ( 
  CT_Epetra_Operator_ID_t selfID );
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

