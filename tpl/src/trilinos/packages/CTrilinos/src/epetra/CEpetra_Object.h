#ifndef CEPETRA_OBJECT_H
#define CEPETRA_OBJECT_H

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


/*! @file CEpetra_Object.h
 * @brief Wrappers for Epetra_Object */

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
CT_Epetra_Object_ID_t Epetra_Object_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_Object_Generalize ( 
  CT_Epetra_Object_ID_t id );

/*@}*/

/*! @name Epetra_Object constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Object::Epetra_Object(int TracebackModeIn = -1, bool set_label = true)
*/
CT_Epetra_Object_ID_t Epetra_Object_Create ( 
  int TracebackModeIn, boolean set_label );

/*! @brief Wrapper for 
   Epetra_Object::Epetra_Object(const char * const Label, int TracebackModeIn = -1)
*/
CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( 
  const char * const Label, int TracebackModeIn );

/*! @brief Wrapper for 
   Epetra_Object::Epetra_Object(const Epetra_Object& Object)
*/
CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( 
  CT_Epetra_Object_ID_t ObjectID );

/*@}*/

/*! @name Epetra_Object destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Object::~Epetra_Object()
*/
void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID );

/*@}*/

/*! @name Epetra_Object member wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual void Epetra_Object::SetLabel(const char * const Label)
*/
void Epetra_Object_SetLabel ( 
  CT_Epetra_Object_ID_t selfID, const char * const Label );

/*! @brief Wrapper for 
   virtual const char * Epetra_Object::Label() const
*/
const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_Object::ReportError(const string Message, int ErrorCode) const
*/
int Epetra_Object_ReportError ( 
  CT_Epetra_Object_ID_t selfID, const char Message[], int ErrorCode );

/*@}*/

/*! @name Epetra_Object static function wrappers */
/*@{*/

/*! @brief Wrapper for 
   static void Epetra_Object::SetTracebackMode(int TracebackModeValue)
*/
void Epetra_Object_SetTracebackMode ( int TracebackModeValue );

/*! @brief Wrapper for 
   static int Epetra_Object::GetTracebackMode()
*/
int Epetra_Object_GetTracebackMode (  );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_OBJECT_H */

