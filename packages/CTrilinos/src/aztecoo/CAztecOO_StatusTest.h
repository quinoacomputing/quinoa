#ifndef CAZTECOO_STATUSTEST_H
#define CAZTECOO_STATUSTEST_H

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


#ifdef HAVE_CTRILINOS_AZTECOO



/*! @file CAztecOO_StatusTest.h
 * @brief Wrappers for AztecOO_StatusTest */

/* True C header file! */

#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_AztecOO_StatusTest_ID_t AztecOO_StatusTest_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t AztecOO_StatusTest_Generalize ( 
  CT_AztecOO_StatusTest_ID_t id );

/*@}*/

/*! @name AztecOO_StatusTest destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual AztecOO_StatusTest::~AztecOO_StatusTest()
*/
void AztecOO_StatusTest_Destroy ( 
  CT_AztecOO_StatusTest_ID_t * selfID );

/*@}*/

/*! @name AztecOO_StatusTest member wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual bool AztecOO_StatusTest::ResidualVectorRequired() const = 0
*/
boolean AztecOO_StatusTest_ResidualVectorRequired ( 
  CT_AztecOO_StatusTest_ID_t selfID );

/*! @brief Wrapper for 
   virtual AztecOO_StatusType AztecOO_StatusTest::CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst, bool SolutionUpdated) = 0
*/
CT_AztecOO_StatusType_E_t AztecOO_StatusTest_CheckStatus ( 
  CT_AztecOO_StatusTest_ID_t selfID, int CurrentIter, 
  CT_Epetra_MultiVector_ID_t CurrentResVectorID, 
  double CurrentResNormEst, boolean SolutionUpdated );

/*! @brief Wrapper for 
   virtual AztecOO_StatusType AztecOO_StatusTest::GetStatus() const = 0
*/
CT_AztecOO_StatusType_E_t AztecOO_StatusTest_GetStatus ( 
  CT_AztecOO_StatusTest_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* HAVE_CTRILINOS_AZTECOO */
#endif /* CAZTECOO_STATUSTEST_H */

