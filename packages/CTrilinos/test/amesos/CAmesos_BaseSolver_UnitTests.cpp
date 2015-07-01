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


#ifdef HAVE_CTRILINOS_AMESOS

#include "CEpetra_Comm.h"
#include "CEpetra_LinearProblem.h"
#include "CTeuchos_ParameterList.h"
#include "Amesos_BaseSolver.h"
#include "CAmesos_BaseSolver.h"
#include "CAmesos_BaseSolver_Cpp.hpp"
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
void Amesos_BaseSolver_Destroy ( CT_Amesos_BaseSolver_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
int Amesos_BaseSolver_SymbolicFactorization ( 
  CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Amesos_BaseSolver_NumericFactorization ( 
  CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Amesos_BaseSolver_Solve ( CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Amesos_BaseSolver_SetUseTranspose ( 
  CT_Amesos_BaseSolver_ID_t selfID, boolean UseTranspose );
 **********************************************************************/

/**********************************************************************
boolean Amesos_BaseSolver_UseTranspose ( 
  CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Amesos_BaseSolver_SetParameters ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t ParameterListID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_LinearProblem_ID_t Amesos_BaseSolver_GetProblem ( 
  CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Amesos_BaseSolver_MatrixShapeOK ( 
  CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Comm_ID_t Amesos_BaseSolver_Comm ( 
  CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Amesos_BaseSolver_NumSymbolicFact ( 
  CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Amesos_BaseSolver_NumNumericFact ( 
  CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Amesos_BaseSolver_NumSolve ( CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Amesos_BaseSolver_PrintStatus ( 
  CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Amesos_BaseSolver_PrintTiming ( 
  CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Amesos_BaseSolver_setParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t paramListID );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_getNonconstParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_unsetParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID );
 **********************************************************************/

/**********************************************************************
void Amesos_BaseSolver_GetTiming ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t TimingParameterListID );
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

#endif /* HAVE_CTRILINOS_AMESOS */

