#ifndef CAMESOS_BASESOLVER_H
#define CAMESOS_BASESOLVER_H

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



/*! @file CAmesos_BaseSolver.h
 * @brief Wrappers for Amesos_BaseSolver */

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
CT_Amesos_BaseSolver_ID_t Amesos_BaseSolver_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Amesos_BaseSolver_Generalize ( 
  CT_Amesos_BaseSolver_ID_t id );

/*@}*/

/*! @name Amesos_BaseSolver destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Amesos_BaseSolver::~Amesos_BaseSolver()
*/
void Amesos_BaseSolver_Destroy ( CT_Amesos_BaseSolver_ID_t * selfID );

/*@}*/

/*! @name Amesos_BaseSolver member wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual int Amesos_BaseSolver::SymbolicFactorization() = 0
*/
int Amesos_BaseSolver_SymbolicFactorization ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Amesos_BaseSolver::NumericFactorization() = 0
*/
int Amesos_BaseSolver_NumericFactorization ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Amesos_BaseSolver::Solve() = 0
*/
int Amesos_BaseSolver_Solve ( CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Amesos_BaseSolver::SetUseTranspose(bool UseTranspose) = 0
*/
int Amesos_BaseSolver_SetUseTranspose ( 
  CT_Amesos_BaseSolver_ID_t selfID, boolean UseTranspose );

/*! @brief Wrapper for 
   virtual bool Amesos_BaseSolver::UseTranspose() const = 0
*/
boolean Amesos_BaseSolver_UseTranspose ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Amesos_BaseSolver::SetParameters( Teuchos::ParameterList &ParameterList ) = 0
*/
int Amesos_BaseSolver_SetParameters ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t ParameterListID );

/*! @brief Wrapper for 
   virtual const Epetra_LinearProblem* Amesos_BaseSolver::GetProblem() const = 0
*/
CT_Epetra_LinearProblem_ID_t Amesos_BaseSolver_GetProblem ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual bool Amesos_BaseSolver::MatrixShapeOK() const = 0
*/
boolean Amesos_BaseSolver_MatrixShapeOK ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual const Epetra_Comm & Amesos_BaseSolver::Comm() const = 0
*/
CT_Epetra_Comm_ID_t Amesos_BaseSolver_Comm ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Amesos_BaseSolver::NumSymbolicFact() const = 0
*/
int Amesos_BaseSolver_NumSymbolicFact ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Amesos_BaseSolver::NumNumericFact() const = 0
*/
int Amesos_BaseSolver_NumNumericFact ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Amesos_BaseSolver::NumSolve() const = 0
*/
int Amesos_BaseSolver_NumSolve ( CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual void Amesos_BaseSolver::PrintStatus() const = 0
*/
void Amesos_BaseSolver_PrintStatus ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual void Amesos_BaseSolver::PrintTiming() const = 0
*/
void Amesos_BaseSolver_PrintTiming ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual void Amesos_BaseSolver::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
*/
void Amesos_BaseSolver_setParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t paramListID );

/*! @brief Wrapper for 
   virtual Teuchos::RCP<Teuchos::ParameterList> Amesos_BaseSolver::getNonconstParameterList()
*/
CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_getNonconstParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual Teuchos::RCP<Teuchos::ParameterList> Amesos_BaseSolver::unsetParameterList()
*/
CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_unsetParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID );

/*! @brief Wrapper for 
   virtual void Amesos_BaseSolver::GetTiming( Teuchos::ParameterList &TimingParameterList ) const
*/
void Amesos_BaseSolver_GetTiming ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t TimingParameterListID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* HAVE_CTRILINOS_AMESOS */
#endif /* CAMESOS_BASESOLVER_H */

