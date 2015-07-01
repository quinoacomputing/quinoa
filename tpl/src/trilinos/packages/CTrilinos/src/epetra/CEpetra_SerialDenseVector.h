#ifndef CEPETRA_SERIALDENSEVECTOR_H
#define CEPETRA_SERIALDENSEVECTOR_H

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


/*! @file CEpetra_SerialDenseVector.h
 * @brief Wrappers for Epetra_SerialDenseVector */

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
CT_Epetra_SerialDenseVector_ID_t Epetra_SerialDenseVector_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_SerialDenseVector_Generalize ( 
  CT_Epetra_SerialDenseVector_ID_t id );

/*@}*/

/*! @name Epetra_SerialDenseVector constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_SerialDenseVector::Epetra_SerialDenseVector()
*/
CT_Epetra_SerialDenseVector_ID_t Epetra_SerialDenseVector_Create_Empty ( 
   );

/*! @brief Wrapper for 
   Epetra_SerialDenseVector::Epetra_SerialDenseVector(int Length)
*/
CT_Epetra_SerialDenseVector_ID_t Epetra_SerialDenseVector_Create ( 
  int Length );

/*! @brief Wrapper for 
   Epetra_SerialDenseVector::Epetra_SerialDenseVector(Epetra_DataAccess CV, double* Values, int Length)
*/
CT_Epetra_SerialDenseVector_ID_t Epetra_SerialDenseVector_Create_FromArray ( 
  CT_Epetra_DataAccess_E_t CV, double * Values, int Length );

/*! @brief Wrapper for 
   Epetra_SerialDenseVector::Epetra_SerialDenseVector(const Epetra_SerialDenseVector& Source)
*/
CT_Epetra_SerialDenseVector_ID_t Epetra_SerialDenseVector_Duplicate ( 
  CT_Epetra_SerialDenseVector_ID_t SourceID );

/*@}*/

/*! @name Epetra_SerialDenseVector destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_SerialDenseVector::~Epetra_SerialDenseVector()
*/
void Epetra_SerialDenseVector_Destroy ( 
  CT_Epetra_SerialDenseVector_ID_t * selfID );

/*@}*/

/*! @name Epetra_SerialDenseVector member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_SerialDenseVector::Size(int Length_in)
*/
int Epetra_SerialDenseVector_Size ( 
  CT_Epetra_SerialDenseVector_ID_t selfID, int Length_in );

/*! @brief Wrapper for 
   int Epetra_SerialDenseVector::Resize(int Length_in)
*/
int Epetra_SerialDenseVector_Resize ( 
  CT_Epetra_SerialDenseVector_ID_t selfID, int Length_in );

/*! @brief Wrapper for 
   int Epetra_SerialDenseVector::Random()
*/
int Epetra_SerialDenseVector_Random ( 
  CT_Epetra_SerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   double Epetra_SerialDenseVector::Dot(const Epetra_SerialDenseVector & x) const
*/
double Epetra_SerialDenseVector_Dot ( 
  CT_Epetra_SerialDenseVector_ID_t selfID, 
  CT_Epetra_SerialDenseVector_ID_t xID );

/*! @brief Wrapper for 
   double Epetra_SerialDenseVector::Norm1() const
*/
double Epetra_SerialDenseVector_Norm1 ( 
  CT_Epetra_SerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   double Epetra_SerialDenseVector::Norm2() const
*/
double Epetra_SerialDenseVector_Norm2 ( 
  CT_Epetra_SerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   double Epetra_SerialDenseVector::NormInf() const
*/
double Epetra_SerialDenseVector_NormInf ( 
  CT_Epetra_SerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_SerialDenseVector::Length() const
*/
int Epetra_SerialDenseVector_Length ( 
  CT_Epetra_SerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   double* Epetra_SerialDenseVector::Values() const
*/
double * Epetra_SerialDenseVector_Values ( 
  CT_Epetra_SerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_DataAccess Epetra_SerialDenseVector::CV() const
*/
CT_Epetra_DataAccess_E_t Epetra_SerialDenseVector_CV ( 
  CT_Epetra_SerialDenseVector_ID_t selfID );

/*@}*/

/*! @name Epetra_SerialDenseVector operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_SerialDenseVector& Epetra_SerialDenseVector::operator= (const Epetra_SerialDenseVector& Source)
*/
void Epetra_SerialDenseVector_Assign ( 
  CT_Epetra_SerialDenseVector_ID_t selfID, 
  CT_Epetra_SerialDenseVector_ID_t SourceID );

/*! @brief Wrapper for 
   double& Epetra_SerialDenseVector::operator() (int Index)
*/
void Epetra_SerialDenseVector_setElement ( 
  CT_Epetra_SerialDenseVector_ID_t selfID, int Index, 
  double * value );

/*! @brief Wrapper for 
   const double& Epetra_SerialDenseVector::operator() (int Index) const
*/
double Epetra_SerialDenseVector_getElement ( 
  CT_Epetra_SerialDenseVector_ID_t selfID, int Index );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_SERIALDENSEVECTOR_H */

