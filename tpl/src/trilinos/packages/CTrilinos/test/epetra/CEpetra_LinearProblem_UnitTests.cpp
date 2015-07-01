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
#include "CTrilinos_config.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_Operator.h"
#include "CEpetra_RowMatrix.h"
#include "CEpetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "CEpetra_LinearProblem.h"
#include "CEpetra_LinearProblem_Cpp.hpp"
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
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create (  );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_LinearProblem , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_LinearProblem_ID_t selfID = Epetra_LinearProblem_Create());

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_LinearProblem_ID);
}

/**********************************************************************
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromMatrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromOperator ( 
  CT_Epetra_Operator_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Duplicate ( 
  CT_Epetra_LinearProblem_ID_t ProblemID );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_Destroy ( 
  CT_Epetra_LinearProblem_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_LinearProblem_CheckInput ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_AssertSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_SetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_ProblemDifficultyLevel_E_t PDL );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_SetOperator_Matrix ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_RowMatrix_ID_t AID );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_SetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Operator_ID_t AID );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_SetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t XID );
 **********************************************************************/

/**********************************************************************
void Epetra_LinearProblem_SetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t BID );
 **********************************************************************/

/**********************************************************************
int Epetra_LinearProblem_LeftScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID );
 **********************************************************************/

/**********************************************************************
int Epetra_LinearProblem_RightScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Operator_ID_t Epetra_LinearProblem_GetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_RowMatrix_ID_t Epetra_LinearProblem_GetMatrix ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_ProblemDifficultyLevel_E_t Epetra_LinearProblem_GetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_LinearProblem_IsOperatorSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID );
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

