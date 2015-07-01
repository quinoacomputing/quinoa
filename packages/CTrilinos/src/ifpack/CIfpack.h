#ifndef CIFPACK_H
#define CIFPACK_H

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


#ifdef HAVE_CTRILINOS_IFPACK



/*! @file CIfpack.h
 * @brief Wrappers for Ifpack */

/* True C header file! */

#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name Ifpack constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Ifpack::Ifpack()
*/
CT_Ifpack_ID_t Ifpack_Create (  );

/*@}*/

/*! @name Ifpack destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Ifpack::~Ifpack()
*/
void Ifpack_Destroy ( CT_Ifpack_ID_t * selfID );

/*@}*/

/*! @name Ifpack member wrappers */
/*@{*/

/*! @brief Wrapper for 
   Ifpack_Preconditioner* Ifpack::Create(const string PrecType, Epetra_RowMatrix* Matrix, const int overlap = 0)
*/
CT_Ifpack_Preconditioner_ID_t Ifpack_CreatePreconditioner_UsingName ( 
  CT_Ifpack_ID_t selfID, const char PrecType[], 
  CT_Epetra_RowMatrix_ID_t MatrixID, const int overlap );

/*! @brief Wrapper for 
   int Ifpack::SetParameters(int argc, char* argv[], Teuchos::ParameterList& List, string& PrecType, int& Overlap)
*/
int Ifpack_SetParameters ( 
  CT_Ifpack_ID_t selfID, int argc, char * argv[], 
  CT_Teuchos_ParameterList_ID_t ListID, char * PrecType[], 
  int * Overlap );

/*@}*/

/*! @name Ifpack static function wrappers */
/*@{*/

/*! @brief Wrapper for 
   static const char* Ifpack::toString(const EPrecType precType)
*/
const char * Ifpack_toString ( const CT_EPrecType_E_t precType );

/*! @brief Wrapper for 
   static Ifpack_Preconditioner* Ifpack::Create( EPrecType PrecType, Epetra_RowMatrix* Matrix, const int overlap = 0 )
*/
CT_Ifpack_Preconditioner_ID_t Ifpack_CreatePreconditioner_UsingType ( 
  CT_EPrecType_E_t PrecType, CT_Epetra_RowMatrix_ID_t MatrixID, 
  const int overlap );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* HAVE_CTRILINOS_IFPACK */
#endif /* CIFPACK_H */

