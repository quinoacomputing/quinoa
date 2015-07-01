#ifndef CEPETRA_COMPOBJECT_H
#define CEPETRA_COMPOBJECT_H

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


/*! @file CEpetra_CompObject.h
 * @brief Wrappers for Epetra_CompObject */

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
CT_Epetra_CompObject_ID_t Epetra_CompObject_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_CompObject_Generalize ( 
  CT_Epetra_CompObject_ID_t id );

/*@}*/

/*! @name Epetra_CompObject constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_CompObject::Epetra_CompObject()
*/
CT_Epetra_CompObject_ID_t Epetra_CompObject_Create (  );

/*! @brief Wrapper for 
   Epetra_CompObject::Epetra_CompObject(const Epetra_CompObject& Source)
*/
CT_Epetra_CompObject_ID_t Epetra_CompObject_Duplicate ( 
  CT_Epetra_CompObject_ID_t SourceID );

/*@}*/

/*! @name Epetra_CompObject destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_CompObject::~Epetra_CompObject()
*/
void Epetra_CompObject_Destroy ( CT_Epetra_CompObject_ID_t * selfID );

/*@}*/

/*! @name Epetra_CompObject member wrappers */
/*@{*/

/*! @brief Wrapper for 
   void Epetra_CompObject::SetFlopCounter(const Epetra_Flops & FlopCounter_in)
*/
void Epetra_CompObject_SetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_Flops_ID_t FlopCounter_inID );

/*! @brief Wrapper for 
   void Epetra_CompObject::SetFlopCounter(const Epetra_CompObject & CompObject)
*/
void Epetra_CompObject_SetFlopCounter_Matching ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_CompObject_ID_t CompObjectID );

/*! @brief Wrapper for 
   void Epetra_CompObject::UnsetFlopCounter()
*/
void Epetra_CompObject_UnsetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_Flops * Epetra_CompObject::GetFlopCounter() const
*/
CT_Epetra_Flops_ID_t Epetra_CompObject_GetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID );

/*! @brief Wrapper for 
   void Epetra_CompObject::ResetFlops() const
*/
void Epetra_CompObject_ResetFlops ( 
  CT_Epetra_CompObject_ID_t selfID );

/*! @brief Wrapper for 
   double Epetra_CompObject::Flops() const
*/
double Epetra_CompObject_Flops ( CT_Epetra_CompObject_ID_t selfID );

/*! @brief Wrapper for 
   void Epetra_CompObject::UpdateFlops(int Flops_in) const
*/
void Epetra_CompObject_UpdateFlops_Int ( 
  CT_Epetra_CompObject_ID_t selfID, int Flops_in );

/*! @brief Wrapper for 
   void Epetra_CompObject::UpdateFlops(long int Flops_in) const
*/
void Epetra_CompObject_UpdateFlops_Long ( 
  CT_Epetra_CompObject_ID_t selfID, long int Flops_in );

/*! @brief Wrapper for 
   void Epetra_CompObject::UpdateFlops(double Flops_in) const
*/
void Epetra_CompObject_UpdateFlops_Double ( 
  CT_Epetra_CompObject_ID_t selfID, double Flops_in );

/*! @brief Wrapper for 
   void Epetra_CompObject::UpdateFlops(float Flops_in) const
*/
void Epetra_CompObject_UpdateFlops_Float ( 
  CT_Epetra_CompObject_ID_t selfID, float Flops_in );

/*@}*/

/*! @name Epetra_CompObject operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_CompObject& Epetra_CompObject::operator=(const Epetra_CompObject& src)
*/
void Epetra_CompObject_Assign ( 
  CT_Epetra_CompObject_ID_t selfID, CT_Epetra_CompObject_ID_t srcID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_COMPOBJECT_H */

