#ifndef CEPETRA_LINEARPROBLEM_H
#define CEPETRA_LINEARPROBLEM_H

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


/*! @file CEpetra_LinearProblem.h
 * @brief Wrappers for Epetra_LinearProblem */

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
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_LinearProblem_Generalize ( 
  CT_Epetra_LinearProblem_ID_t id );

/*@}*/

/*! @name Epetra_LinearProblem constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_LinearProblem::Epetra_LinearProblem(void)
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create (  );

/*! @brief Wrapper for 
   Epetra_LinearProblem::Epetra_LinearProblem(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B)
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromMatrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   Epetra_LinearProblem::Epetra_LinearProblem(Epetra_Operator * A, Epetra_MultiVector * X, Epetra_MultiVector * B)
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromOperator ( 
  CT_Epetra_Operator_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   Epetra_LinearProblem::Epetra_LinearProblem(const Epetra_LinearProblem& Problem)
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Duplicate ( 
  CT_Epetra_LinearProblem_ID_t ProblemID );

/*@}*/

/*! @name Epetra_LinearProblem destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_LinearProblem::~Epetra_LinearProblem(void)
*/
void Epetra_LinearProblem_Destroy ( 
  CT_Epetra_LinearProblem_ID_t * selfID );

/*@}*/

/*! @name Epetra_LinearProblem member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_LinearProblem::CheckInput() const
*/
int Epetra_LinearProblem_CheckInput ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   void Epetra_LinearProblem::AssertSymmetric()
*/
void Epetra_LinearProblem_AssertSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   void Epetra_LinearProblem::SetPDL(ProblemDifficultyLevel PDL)
*/
void Epetra_LinearProblem_SetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_ProblemDifficultyLevel_E_t PDL );

/*! @brief Wrapper for 
   void Epetra_LinearProblem::SetOperator(Epetra_RowMatrix * A)
*/
void Epetra_LinearProblem_SetOperator_Matrix ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_RowMatrix_ID_t AID );

/*! @brief Wrapper for 
   void Epetra_LinearProblem::SetOperator(Epetra_Operator * A)
*/
void Epetra_LinearProblem_SetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Operator_ID_t AID );

/*! @brief Wrapper for 
   void Epetra_LinearProblem::SetLHS(Epetra_MultiVector * X)
*/
void Epetra_LinearProblem_SetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t XID );

/*! @brief Wrapper for 
   void Epetra_LinearProblem::SetRHS(Epetra_MultiVector * B)
*/
void Epetra_LinearProblem_SetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   int Epetra_LinearProblem::LeftScale(const Epetra_Vector & D)
*/
int Epetra_LinearProblem_LeftScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID );

/*! @brief Wrapper for 
   int Epetra_LinearProblem::RightScale(const Epetra_Vector & D)
*/
int Epetra_LinearProblem_RightScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID );

/*! @brief Wrapper for 
   Epetra_Operator * Epetra_LinearProblem::GetOperator() const
*/
CT_Epetra_Operator_ID_t Epetra_LinearProblem_GetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_RowMatrix * Epetra_LinearProblem::GetMatrix() const
*/
CT_Epetra_RowMatrix_ID_t Epetra_LinearProblem_GetMatrix ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_MultiVector * Epetra_LinearProblem::GetLHS() const
*/
CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_MultiVector * Epetra_LinearProblem::GetRHS() const
*/
CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   ProblemDifficultyLevel Epetra_LinearProblem::GetPDL() const
*/
CT_ProblemDifficultyLevel_E_t Epetra_LinearProblem_GetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_LinearProblem::IsOperatorSymmetric() const
*/
boolean Epetra_LinearProblem_IsOperatorSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_LINEARPROBLEM_H */

