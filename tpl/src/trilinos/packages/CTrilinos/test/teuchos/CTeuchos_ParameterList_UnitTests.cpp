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
#include "CTeuchos_ParameterEntry.h"
#include "Teuchos_ParameterList.hpp"
#include "CTeuchos_ParameterList.h"
#include "CTeuchos_ParameterList_Cpp.hpp"
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
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create (  );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Teuchos_ParameterList , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Teuchos_ParameterList_ID_t selfID = Teuchos_ParameterList_Create());

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Teuchos_ParameterList_ID);
}

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_WithName ( 
  const char name[] );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_FromSource ( 
  CT_Teuchos_ParameterList_ID_t sourceID );
 **********************************************************************/

/**********************************************************************
void Teuchos_ParameterList_Destroy ( 
  CT_Teuchos_ParameterList_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setName ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
void Teuchos_ParameterList_Assign ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParametersNotAlreadySet ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_disableRecursiveValidation ( 
  CT_Teuchos_ParameterList_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  CT_Teuchos_ParameterList_ID_t valueID, char const docString[] );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  CT_Teuchos_ParameterEntry_ID_t entryID );
 **********************************************************************/

/**********************************************************************
double Teuchos_ParameterList_get_double_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double def_value );
 **********************************************************************/

/**********************************************************************
int Teuchos_ParameterList_get_int_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int def_value );
 **********************************************************************/

/**********************************************************************
const char * Teuchos_ParameterList_get_char_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  char def_value[] );
 **********************************************************************/

/**********************************************************************
const char * Teuchos_ParameterList_get_const_char_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  const char def_value[] );
 **********************************************************************/

/**********************************************************************
double Teuchos_ParameterList_get_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
int Teuchos_ParameterList_get_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
double Teuchos_ParameterList_get_double_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
int Teuchos_ParameterList_get_int_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
double * Teuchos_ParameterList_getPtr_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
int * Teuchos_ParameterList_getPtr_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
const double * Teuchos_ParameterList_getPtr_double_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
const int * Teuchos_ParameterList_getPtr_int_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_ParameterList_remove ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  boolean throwIfNotExists );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  boolean mustAlreadyExist, const char docString[] );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist_existing ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
const char * Teuchos_ParameterList_name_it ( 
  CT_Teuchos_ParameterList_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_ParameterList_isParameter ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_ParameterList_isSublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_ParameterList_isType_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_ParameterList_isType_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_ParameterList_isType_double_type ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double * ptr );
 **********************************************************************/

/**********************************************************************
boolean Teuchos_ParameterList_isType_int_type ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int * ptr );
 **********************************************************************/

/**********************************************************************
const char * Teuchos_ParameterList_currentParametersString ( 
  CT_Teuchos_ParameterList_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Teuchos_ParameterList_validateParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth, 
  const CT_EValidateUsed_E_t validateUsed, 
  const CT_EValidateDefaults_E_t validateDefaults );
 **********************************************************************/

/**********************************************************************
void Teuchos_ParameterList_validateParametersAndSetDefaults ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth );
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

