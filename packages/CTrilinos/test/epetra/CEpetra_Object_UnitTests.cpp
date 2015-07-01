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

#include <cstring>

#include "Epetra_Object.h"
#include "CEpetra_Object.h"
#include "CEpetra_Object_Cpp.hpp"
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
CT_Epetra_Object_ID_t Epetra_Object_Create ( 
  int TracebackModeIn, boolean set_label );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(int TracebackModeIn = -1);
  ECHO(boolean set_label = TRUE);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create(TracebackModeIn, set_label));

  TEST_EQUALITY(selfID.table, CT_Epetra_Object_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( 
  const char * const Label, int TracebackModeIn );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , Create_WithLabel )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(const char * Label = "blah");
  ECHO(int TracebackModeIn = -1);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create_WithLabel(Label, TracebackModeIn));

  TEST_EQUALITY(selfID.table, CT_Epetra_Object_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( 
  CT_Epetra_Object_ID_t ObjectID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(int TracebackModeIn = -1);
  ECHO(boolean set_label = TRUE);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create(TracebackModeIn, set_label));

  ECHO(CT_Epetra_Object_ID_t dupID = Epetra_Object_Duplicate(selfID));

  TEST_EQUALITY(dupID.table, CT_Epetra_Object_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(int TracebackModeIn = -1);
  ECHO(boolean set_label = TRUE);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create(TracebackModeIn, set_label));

  ECHO(Epetra_Object_Destroy(&selfID));

  TEST_EQUALITY(selfID.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
void Epetra_Object_SetLabel ( 
  CT_Epetra_Object_ID_t selfID, const char * const Label );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , SetLabel )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(int TracebackModeIn = -1);
  ECHO(boolean set_label = TRUE);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create(TracebackModeIn, set_label));

  ECHO(const char * Label = "blah");
  ECHO(Epetra_Object_SetLabel(selfID, Label));

  ECHO(const char * Label2 = Epetra_Object_Label(selfID));
  TEST_EQUALITY_CONST(strcmp(Label, Label2), 0);
}

/**********************************************************************
const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , Label )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(const char * Label = "blah");
  ECHO(int TracebackModeIn = -1);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create_WithLabel(Label, TracebackModeIn));

  ECHO(const char * Label2 = Epetra_Object_Label(selfID));
  TEST_EQUALITY_CONST(strcmp(Label, Label2), 0);
}

/**********************************************************************
int Epetra_Object_ReportError ( 
  CT_Epetra_Object_ID_t selfID, const char * Message, int ErrorCode );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , ReportError )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(const char * Label = "This is ONLY a test object; please ignore errors");
  ECHO(int TracebackModeIn = 1);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create_WithLabel(Label, TracebackModeIn));

  ECHO(const char * msg = "This is ONLY a test error message; please ignore");
  ECHO(int code = -5);
  ECHO(int ret = Epetra_Object_ReportError(selfID, msg, code));

  TEST_EQUALITY(ret, code);
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

