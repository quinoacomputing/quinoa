#ifndef CEPETRA_OPERATOR_H
#define CEPETRA_OPERATOR_H

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


/*! @file CEpetra_Operator.h
 * @brief Wrappers for Epetra_Operator */

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
CT_Epetra_Operator_ID_t Epetra_Operator_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_Operator_Generalize ( 
  CT_Epetra_Operator_ID_t id );

/*@}*/

/*! @name Epetra_Operator destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Operator::~Epetra_Operator()
*/
void Epetra_Operator_Destroy ( CT_Epetra_Operator_ID_t * selfID );

/*@}*/

/*! @name Epetra_Operator member wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual int Epetra_Operator::SetUseTranspose(bool UseTranspose) = 0
*/
int Epetra_Operator_SetUseTranspose ( 
  CT_Epetra_Operator_ID_t selfID, boolean UseTranspose );

/*! @brief Wrapper for 
   virtual int Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0
*/
int Epetra_Operator_Apply ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/*! @brief Wrapper for 
   virtual int Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0
*/
int Epetra_Operator_ApplyInverse ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/*! @brief Wrapper for 
   virtual double Epetra_Operator::NormInf() const = 0
*/
double Epetra_Operator_NormInf ( CT_Epetra_Operator_ID_t selfID );

/*! @brief Wrapper for 
   virtual const char * Epetra_Operator::Label() const = 0
*/
const char * Epetra_Operator_Label ( CT_Epetra_Operator_ID_t selfID );

/*! @brief Wrapper for 
   virtual bool Epetra_Operator::UseTranspose() const = 0
*/
boolean Epetra_Operator_UseTranspose ( 
  CT_Epetra_Operator_ID_t selfID );

/*! @brief Wrapper for 
   virtual bool Epetra_Operator::HasNormInf() const = 0
*/
boolean Epetra_Operator_HasNormInf ( CT_Epetra_Operator_ID_t selfID );

/*! @brief Wrapper for 
   virtual const Epetra_Comm & Epetra_Operator::Comm() const = 0
*/
CT_Epetra_Comm_ID_t Epetra_Operator_Comm ( 
  CT_Epetra_Operator_ID_t selfID );

/*! @brief Wrapper for 
   virtual const Epetra_Map & Epetra_Operator::OperatorDomainMap() const = 0
*/
CT_Epetra_Map_ID_t Epetra_Operator_OperatorDomainMap ( 
  CT_Epetra_Operator_ID_t selfID );

/*! @brief Wrapper for 
   virtual const Epetra_Map & Epetra_Operator::OperatorRangeMap() const = 0
*/
CT_Epetra_Map_ID_t Epetra_Operator_OperatorRangeMap ( 
  CT_Epetra_Operator_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_OPERATOR_H */

