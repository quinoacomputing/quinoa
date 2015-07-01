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
#include "CTeuchos_ParameterList.h"
#include "CTeuchos_any.h"
#include "Teuchos_ParameterEntry.hpp"
#include "CTeuchos_ParameterEntry.h"
#include "CTeuchos_ParameterEntry_Cpp.hpp"
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
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Create (  );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Teuchos_ParameterEntry , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Teuchos_ParameterEntry_ID_t selfID = Teuchos_ParameterEntry_Create());

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Teuchos_ParameterEntry_ID);
}

/**********************************************************************
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Duplicate ( 
  CT_Teuchos_ParameterEntry_ID_t sourceID );
 **********************************************************************/

/**********************************************************************
void Teuchos_ParameterEntry_Destroy ( 
  CT_Teuchos_ParameterEntry_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
void Teuchos_ParameterEntry_Assign ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, 
  CT_Teuchos_ParameterEntry_ID_t sourceID );
 **********************************************************************/

/**********************************************************************
void Teuchos_ParameterEntry_setAnyValue ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, 
  CT_Teuchos_any_ID_t valueID, boolean isDefault );
 **********************************************************************/

/**********************************************************************
void Teuchos_ParameterEntry_setDocString ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, const char docString[] );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterEntry_setList ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean isDefault, 
  const char docString[] );
 **********************************************************************/

/**********************************************************************
double Teuchos_ParameterEntry_getValue_double ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, double * ptr );
 **********************************************************************/

/**********************************************************************
int Teuchos_ParameterEntry_getValue_int ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, int * ptr );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean activeQry );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny_const ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean activeQry );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_ParameterEntry_isUsed ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_ParameterEntry_isList ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_ParameterEntry_isType_double ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_ParameterEntry_isType_int ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_ParameterEntry_isDefault ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );
 **********************************************************************/

/**********************************************************************
const char * Teuchos_ParameterEntry_docString ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );
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

