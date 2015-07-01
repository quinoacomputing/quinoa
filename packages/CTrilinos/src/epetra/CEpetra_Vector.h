#ifndef CEPETRA_VECTOR_H
#define CEPETRA_VECTOR_H

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


/*! @file CEpetra_Vector.h
 * @brief Wrappers for Epetra_Vector */

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
CT_Epetra_Vector_ID_t Epetra_Vector_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_Vector_Generalize ( 
  CT_Epetra_Vector_ID_t id );

/*@}*/

/*! @name Epetra_Vector constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Vector::Epetra_Vector(const Epetra_BlockMap& Map, bool zeroOut = true)
*/
CT_Epetra_Vector_ID_t Epetra_Vector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut );

/*! @brief Wrapper for 
   Epetra_Vector::Epetra_Vector(const Epetra_Vector& Source)
*/
CT_Epetra_Vector_ID_t Epetra_Vector_Duplicate ( 
  CT_Epetra_Vector_ID_t SourceID );

/*! @brief Wrapper for 
   Epetra_Vector::Epetra_Vector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double *V)
*/
CT_Epetra_Vector_ID_t Epetra_Vector_Create_FromArray ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, 
  double * V );

/*! @brief Wrapper for 
   Epetra_Vector::Epetra_Vector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, int Index)
*/
CT_Epetra_Vector_ID_t Epetra_Vector_FromSource ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int Index );

/*@}*/

/*! @name Epetra_Vector destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Vector::~Epetra_Vector()
*/
void Epetra_Vector_Destroy ( CT_Epetra_Vector_ID_t * selfID );

/*@}*/

/*! @name Epetra_Vector member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_Vector::ReplaceGlobalValues(int NumEntries, double * Values, int * Indices)
*/
int Epetra_Vector_ReplaceGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::ReplaceMyValues(int NumEntries, double * Values, int * Indices)
*/
int Epetra_Vector_ReplaceMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::SumIntoGlobalValues(int NumEntries, double * Values, int * Indices)
*/
int Epetra_Vector_SumIntoGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::SumIntoMyValues(int NumEntries, double * Values, int * Indices)
*/
int Epetra_Vector_SumIntoMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::ReplaceGlobalValues(int NumEntries, int BlockOffset, double * Values, int * Indices)
*/
int Epetra_Vector_ReplaceGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::ReplaceMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices)
*/
int Epetra_Vector_ReplaceMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::SumIntoGlobalValues(int NumEntries, int BlockOffset, double * Values, int * Indices)
*/
int Epetra_Vector_SumIntoGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::SumIntoMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices)
*/
int Epetra_Vector_SumIntoMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::ExtractCopy(double *V) const
*/
int Epetra_Vector_ExtractCopy ( 
  CT_Epetra_Vector_ID_t selfID, double * V );

/*! @brief Wrapper for 
   int Epetra_Vector::ExtractView(double **V) const
*/
int Epetra_Vector_ExtractView ( 
  CT_Epetra_Vector_ID_t selfID, double ** V );

/*@}*/

/*! @name Epetra_Vector operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   const double& Epetra_Vector::operator[] (int index) const
*/
double Epetra_Vector_getElement ( 
  CT_Epetra_Vector_ID_t selfID, int index );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_VECTOR_H */

