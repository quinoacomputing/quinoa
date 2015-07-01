#ifndef CEPETRA_INTSERIALDENSEVECTOR_H
#define CEPETRA_INTSERIALDENSEVECTOR_H

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


/*! @file CEpetra_IntSerialDenseVector.h
 * @brief Wrappers for Epetra_IntSerialDenseVector */

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
CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_IntSerialDenseVector_Generalize ( 
  CT_Epetra_IntSerialDenseVector_ID_t id );

/*@}*/

/*! @name Epetra_IntSerialDenseVector constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector()
*/
CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Create_Empty ( 
   );

/*! @brief Wrapper for 
   Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(int Length_in)
*/
CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Create ( 
  int Length_in );

/*! @brief Wrapper for 
   Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(Epetra_DataAccess CV_in, int* Values_in, int Length_in)
*/
CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Create_FromArray ( 
  CT_Epetra_DataAccess_E_t CV_in, int * Values_in, int Length_in );

/*! @brief Wrapper for 
   Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(const Epetra_IntSerialDenseVector& Source)
*/
CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Duplicate ( 
  CT_Epetra_IntSerialDenseVector_ID_t SourceID );

/*@}*/

/*! @name Epetra_IntSerialDenseVector destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_IntSerialDenseVector::~Epetra_IntSerialDenseVector()
*/
void Epetra_IntSerialDenseVector_Destroy ( 
  CT_Epetra_IntSerialDenseVector_ID_t * selfID );

/*@}*/

/*! @name Epetra_IntSerialDenseVector member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_IntSerialDenseVector::Size(int Length_in)
*/
int Epetra_IntSerialDenseVector_Size ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Length_in );

/*! @brief Wrapper for 
   int Epetra_IntSerialDenseVector::Resize(int Length_in)
*/
int Epetra_IntSerialDenseVector_Resize ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Length_in );

/*! @brief Wrapper for 
   int Epetra_IntSerialDenseVector::Random()
*/
int Epetra_IntSerialDenseVector_Random ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_IntSerialDenseVector::Length() const
*/
int Epetra_IntSerialDenseVector_Length ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   int* Epetra_IntSerialDenseVector::Values()
*/
int * Epetra_IntSerialDenseVector_Values ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   const int* Epetra_IntSerialDenseVector::Values() const
*/
const int * Epetra_IntSerialDenseVector_Values_Const ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_DataAccess Epetra_IntSerialDenseVector::CV() const
*/
CT_Epetra_DataAccess_E_t Epetra_IntSerialDenseVector_CV ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_IntSerialDenseVector::MakeViewOf(const Epetra_IntSerialDenseVector& Source)
*/
int Epetra_IntSerialDenseVector_MakeViewOf ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t SourceID );

/*@}*/

/*! @name Epetra_IntSerialDenseVector operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   int& Epetra_IntSerialDenseVector::operator() (int Index)
*/
void Epetra_IntSerialDenseVector_setElement ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Index, 
  int * value );

/*! @brief Wrapper for 
   const int& Epetra_IntSerialDenseVector::operator() (int Index) const
*/
int Epetra_IntSerialDenseVector_getElement ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Index );

/*! @brief Wrapper for 
   int& Epetra_IntSerialDenseVector::operator[] (int Index)
*/
void Epetra_IntSerialDenseVector_setElement_Bracket ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Index, 
  int * value );

/*! @brief Wrapper for 
   const int& Epetra_IntSerialDenseVector::operator[] (int Index) const
*/
int Epetra_IntSerialDenseVector_getElement_Bracket ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Index );

/*! @brief Wrapper for 
   Epetra_IntSerialDenseVector& Epetra_IntSerialDenseVector::operator= (const Epetra_IntSerialDenseVector& Source)
*/
void Epetra_IntSerialDenseVector_Assign ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t SourceID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_INTSERIALDENSEVECTOR_H */

