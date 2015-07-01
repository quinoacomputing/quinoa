#ifndef CEPETRA_MULTIVECTOR_H
#define CEPETRA_MULTIVECTOR_H

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


/*! @file CEpetra_MultiVector.h
 * @brief Wrappers for Epetra_MultiVector */

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
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_MultiVector_Generalize ( 
  CT_Epetra_MultiVector_ID_t id );

/*@}*/

/*! @name Epetra_MultiVector constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_MultiVector::Epetra_MultiVector(const Epetra_BlockMap& Map, int NumVectors, bool zeroOut = true)
*/
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, int NumVectors, boolean zeroOut );

/*! @brief Wrapper for 
   Epetra_MultiVector::Epetra_MultiVector(const Epetra_MultiVector& Source)
*/
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( 
  CT_Epetra_MultiVector_ID_t SourceID );

/*! @brief Wrapper for 
   Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double *A, int MyLDA, int NumVectors)
*/
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_From2DA ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, 
  double * A, int MyLDA, int NumVectors );

/*! @brief Wrapper for 
   Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double **ArrayOfPointers, int NumVectors)
*/
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromAOP ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, 
  double ** ArrayOfPointers, int NumVectors );

/*! @brief Wrapper for 
   Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, int *Indices, int NumVectors)
*/
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromList ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int * Indices, int NumVectors );

/*! @brief Wrapper for 
   Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, int StartIndex, int NumVectors)
*/
CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromRange ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int StartIndex, int NumVectors );

/*@}*/

/*! @name Epetra_MultiVector destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_MultiVector::~Epetra_MultiVector()
*/
void Epetra_MultiVector_Destroy ( 
  CT_Epetra_MultiVector_ID_t * selfID );

/*@}*/

/*! @name Epetra_MultiVector member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_MultiVector::ReplaceGlobalValue(int GlobalRow, int VectorIndex, double ScalarValue)
*/
int Epetra_MultiVector_ReplaceGlobalValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, int VectorIndex, 
  double ScalarValue );

/*! @brief Wrapper for 
   int Epetra_MultiVector::ReplaceGlobalValue(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue)
*/
int Epetra_MultiVector_ReplaceGlobalValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue );

/*! @brief Wrapper for 
   int Epetra_MultiVector::SumIntoGlobalValue(int GlobalRow, int VectorIndex, double ScalarValue)
*/
int Epetra_MultiVector_SumIntoGlobalValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, int VectorIndex, 
  double ScalarValue );

/*! @brief Wrapper for 
   int Epetra_MultiVector::SumIntoGlobalValue(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue)
*/
int Epetra_MultiVector_SumIntoGlobalValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue );

/*! @brief Wrapper for 
   int Epetra_MultiVector::ReplaceMyValue(int MyRow, int VectorIndex, double ScalarValue)
*/
int Epetra_MultiVector_ReplaceMyValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyRow, int VectorIndex, 
  double ScalarValue );

/*! @brief Wrapper for 
   int Epetra_MultiVector::ReplaceMyValue(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue)
*/
int Epetra_MultiVector_ReplaceMyValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue );

/*! @brief Wrapper for 
   int Epetra_MultiVector::SumIntoMyValue(int MyRow, int VectorIndex, double ScalarValue)
*/
int Epetra_MultiVector_SumIntoMyValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyRow, int VectorIndex, 
  double ScalarValue );

/*! @brief Wrapper for 
   int Epetra_MultiVector::SumIntoMyValue(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue)
*/
int Epetra_MultiVector_SumIntoMyValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue );

/*! @brief Wrapper for 
   int Epetra_MultiVector::PutScalar(double ScalarConstant)
*/
int Epetra_MultiVector_PutScalar ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarConstant );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Random()
*/
int Epetra_MultiVector_Random ( CT_Epetra_MultiVector_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_MultiVector::ExtractCopy(double *A, int MyLDA) const
*/
int Epetra_MultiVector_ExtractCopy_Fill2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double * A, int MyLDA );

/*! @brief Wrapper for 
   int Epetra_MultiVector::ExtractCopy(double **ArrayOfPointers) const
*/
int Epetra_MultiVector_ExtractCopy_FillAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** ArrayOfPointers );

/*! @brief Wrapper for 
   int Epetra_MultiVector::ExtractView(double **A, int *MyLDA) const
*/
int Epetra_MultiVector_ExtractView_Set2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** A, int * MyLDA );

/*! @brief Wrapper for 
   int Epetra_MultiVector::ExtractView(double ***ArrayOfPointers) const
*/
int Epetra_MultiVector_ExtractView_SetAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double *** ArrayOfPointers );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Dot(const Epetra_MultiVector& A, double *Result) const
*/
int Epetra_MultiVector_Dot ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID, 
  double * Result );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Abs(const Epetra_MultiVector& A)
*/
int Epetra_MultiVector_Abs ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Reciprocal(const Epetra_MultiVector& A)
*/
int Epetra_MultiVector_Reciprocal ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Scale(double ScalarValue)
*/
int Epetra_MultiVector_Scale_Self ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarValue );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Scale(double ScalarA, const Epetra_MultiVector& A)
*/
int Epetra_MultiVector_Scale ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Update(double ScalarA, const Epetra_MultiVector& A, double ScalarThis)
*/
int Epetra_MultiVector_Update_WithA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarThis );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Update(double ScalarA, const Epetra_MultiVector& A, double ScalarB, const Epetra_MultiVector& B, double ScalarThis)
*/
int Epetra_MultiVector_Update_WithAB ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarB, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Norm1(double * Result) const
*/
int Epetra_MultiVector_Norm1 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Norm2(double * Result) const
*/
int Epetra_MultiVector_Norm2 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );

/*! @brief Wrapper for 
   int Epetra_MultiVector::NormInf(double * Result) const
*/
int Epetra_MultiVector_NormInf ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );

/*! @brief Wrapper for 
   int Epetra_MultiVector::NormWeighted(const Epetra_MultiVector& Weights, double * Result) const
*/
int Epetra_MultiVector_NormWeighted ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t WeightsID, double * Result );

/*! @brief Wrapper for 
   int Epetra_MultiVector::MinValue(double * Result) const
*/
int Epetra_MultiVector_MinValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );

/*! @brief Wrapper for 
   int Epetra_MultiVector::MaxValue(double * Result) const
*/
int Epetra_MultiVector_MaxValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );

/*! @brief Wrapper for 
   int Epetra_MultiVector::MeanValue(double * Result) const
*/
int Epetra_MultiVector_MeanValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Multiply(char TransA, char TransB, double ScalarAB, const Epetra_MultiVector& A, const Epetra_MultiVector& B, double ScalarThis )
*/
int Epetra_MultiVector_Multiply_Matrix ( 
  CT_Epetra_MultiVector_ID_t selfID, char TransA, char TransB, 
  double ScalarAB, CT_Epetra_MultiVector_ID_t AID, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Multiply(double ScalarAB, const Epetra_MultiVector& A, const Epetra_MultiVector& B, double ScalarThis )
*/
int Epetra_MultiVector_Multiply_ByEl ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis );

/*! @brief Wrapper for 
   int Epetra_MultiVector::ReciprocalMultiply(double ScalarAB, const Epetra_MultiVector& A, const Epetra_MultiVector& B, double ScalarThis )
*/
int Epetra_MultiVector_ReciprocalMultiply ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis );

/*! @brief Wrapper for 
   int Epetra_MultiVector::SetSeed(unsigned int Seed_in)
*/
int Epetra_MultiVector_SetSeed ( 
  CT_Epetra_MultiVector_ID_t selfID, unsigned int Seed_in );

/*! @brief Wrapper for 
   unsigned int Epetra_MultiVector::Seed()
*/
unsigned int Epetra_MultiVector_Seed ( 
  CT_Epetra_MultiVector_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_MultiVector::NumVectors() const
*/
int Epetra_MultiVector_NumVectors ( 
  CT_Epetra_MultiVector_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_MultiVector::MyLength() const
*/
int Epetra_MultiVector_MyLength ( CT_Epetra_MultiVector_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_MultiVector::GlobalLength() const
*/
int Epetra_MultiVector_GlobalLength ( 
  CT_Epetra_MultiVector_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_MultiVector::Stride() const
*/
int Epetra_MultiVector_Stride ( CT_Epetra_MultiVector_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_MultiVector::ConstantStride() const
*/
boolean Epetra_MultiVector_ConstantStride ( 
  CT_Epetra_MultiVector_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_MultiVector::ReplaceMap(const Epetra_BlockMap& map)
*/
int Epetra_MultiVector_ReplaceMap ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID );

/*@}*/

/*! @name Epetra_MultiVector operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_MultiVector& Epetra_MultiVector::operator= (const Epetra_MultiVector& Source)
*/
void Epetra_MultiVector_Assign ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t SourceID );

/*! @brief Wrapper for 
   double * const & Epetra_MultiVector::operator[] (int i) const
*/
double * Epetra_MultiVector_getArray ( 
  CT_Epetra_MultiVector_ID_t selfID, int i );

/*! @brief Wrapper for 
   const Epetra_Vector * & Epetra_MultiVector::operator() (int i) const
*/
CT_Epetra_Vector_ID_t Epetra_MultiVector_getVector ( 
  CT_Epetra_MultiVector_ID_t selfID, int i );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_MULTIVECTOR_H */

