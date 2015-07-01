#ifndef CEPETRA_JADMATRIX_H
#define CEPETRA_JADMATRIX_H

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


/*! @file CEpetra_JadMatrix.h
 * @brief Wrappers for Epetra_JadMatrix */

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
CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_JadMatrix_Generalize ( 
  CT_Epetra_JadMatrix_ID_t id );

/*@}*/

/*! @name Epetra_JadMatrix constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_JadMatrix::Epetra_JadMatrix(const Epetra_RowMatrix & Matrix)
*/
CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Create ( 
  CT_Epetra_RowMatrix_ID_t MatrixID );

/*@}*/

/*! @name Epetra_JadMatrix destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_JadMatrix::~Epetra_JadMatrix()
*/
void Epetra_JadMatrix_Destroy ( CT_Epetra_JadMatrix_ID_t * selfID );

/*@}*/

/*! @name Epetra_JadMatrix member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_JadMatrix::UpdateValues(const Epetra_RowMatrix & Matrix, bool CheckStructure = false)
*/
int Epetra_JadMatrix_UpdateValues ( 
  CT_Epetra_JadMatrix_ID_t selfID, 
  CT_Epetra_RowMatrix_ID_t MatrixID, boolean CheckStructure );

/*! @brief Wrapper for 
   int Epetra_JadMatrix::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const
*/
int Epetra_JadMatrix_ExtractMyRowCopy ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_JadMatrix::ExtractMyEntryView(int CurEntry, double * &Value, int & RowIndex, int & ColIndex)
*/
int Epetra_JadMatrix_ExtractMyEntryView ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, double * * Value, 
  int * RowIndex, int * ColIndex );

/*! @brief Wrapper for 
   int Epetra_JadMatrix::ExtractMyEntryView(int CurEntry, double const * & Value, int & RowIndex, int & ColIndex) const
*/
int Epetra_JadMatrix_ExtractMyEntryView_Const ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, 
  double const ** Value, int * RowIndex, int * ColIndex );

/*! @brief Wrapper for 
   int Epetra_JadMatrix::NumMyRowEntries(int MyRow, int & NumEntries) const
*/
int Epetra_JadMatrix_NumMyRowEntries ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int * NumEntries );

/*! @brief Wrapper for 
   int Epetra_JadMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
*/
int Epetra_JadMatrix_Multiply ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

/*! @brief Wrapper for 
   int Epetra_JadMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
*/
int Epetra_JadMatrix_Solve ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_JADMATRIX_H */

