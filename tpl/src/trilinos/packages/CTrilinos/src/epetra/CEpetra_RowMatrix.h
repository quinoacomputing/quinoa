#ifndef CEPETRA_ROWMATRIX_H
#define CEPETRA_ROWMATRIX_H

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


/*! @file CEpetra_RowMatrix.h
 * @brief Wrappers for Epetra_RowMatrix */

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
CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_RowMatrix_Generalize ( 
  CT_Epetra_RowMatrix_ID_t id );

/*@}*/

/*! @name Epetra_RowMatrix destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_RowMatrix::~Epetra_RowMatrix()
*/
void Epetra_RowMatrix_Destroy ( CT_Epetra_RowMatrix_ID_t * selfID );

/*@}*/

/*! @name Epetra_RowMatrix member wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::NumMyRowEntries(int MyRow, int & NumEntries) const = 0
*/
int Epetra_RowMatrix_NumMyRowEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int * NumEntries );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::MaxNumEntries() const = 0
*/
int Epetra_RowMatrix_MaxNumEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const = 0
*/
int Epetra_RowMatrix_ExtractMyRowCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const = 0
*/
int Epetra_RowMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0
*/
int Epetra_RowMatrix_Multiply ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0
*/
int Epetra_RowMatrix_Solve ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::InvRowSums(Epetra_Vector& x) const = 0
*/
int Epetra_RowMatrix_InvRowSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::LeftScale(const Epetra_Vector& x) = 0
*/
int Epetra_RowMatrix_LeftScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::InvColSums(Epetra_Vector& x) const = 0
*/
int Epetra_RowMatrix_InvColSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::RightScale(const Epetra_Vector& x) = 0
*/
int Epetra_RowMatrix_RightScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/*! @brief Wrapper for 
   virtual bool Epetra_RowMatrix::Filled() const = 0
*/
boolean Epetra_RowMatrix_Filled ( CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual double Epetra_RowMatrix::NormInf() const = 0
*/
double Epetra_RowMatrix_NormInf ( CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual double Epetra_RowMatrix::NormOne() const = 0
*/
double Epetra_RowMatrix_NormOne ( CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::NumGlobalNonzeros() const = 0
*/
int Epetra_RowMatrix_NumGlobalNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::NumGlobalRows() const = 0
*/
int Epetra_RowMatrix_NumGlobalRows ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::NumGlobalCols() const= 0
*/
int Epetra_RowMatrix_NumGlobalCols ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::NumGlobalDiagonals() const = 0
*/
int Epetra_RowMatrix_NumGlobalDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::NumMyNonzeros() const = 0
*/
int Epetra_RowMatrix_NumMyNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::NumMyRows() const = 0
*/
int Epetra_RowMatrix_NumMyRows ( CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::NumMyCols() const = 0
*/
int Epetra_RowMatrix_NumMyCols ( CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_RowMatrix::NumMyDiagonals() const = 0
*/
int Epetra_RowMatrix_NumMyDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual bool Epetra_RowMatrix::LowerTriangular() const = 0
*/
boolean Epetra_RowMatrix_LowerTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual bool Epetra_RowMatrix::UpperTriangular() const = 0
*/
boolean Epetra_RowMatrix_UpperTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual const Epetra_Map & Epetra_RowMatrix::RowMatrixRowMap() const = 0
*/
CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixRowMap ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual const Epetra_Map & Epetra_RowMatrix::RowMatrixColMap() const = 0
*/
CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixColMap ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual const Epetra_Import * Epetra_RowMatrix::RowMatrixImporter() const = 0
*/
CT_Epetra_Import_ID_t Epetra_RowMatrix_RowMatrixImporter ( 
  CT_Epetra_RowMatrix_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_ROWMATRIX_H */

