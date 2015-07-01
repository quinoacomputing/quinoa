#ifndef CEPETRA_FECRSMATRIX_H
#define CEPETRA_FECRSMATRIX_H

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


/*! @file CEpetra_FECrsMatrix.h
 * @brief Wrappers for Epetra_FECrsMatrix */

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
CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_FECrsMatrix_Generalize ( 
  CT_Epetra_FECrsMatrix_ID_t id );

/*@}*/

/*! @name Epetra_FECrsMatrix constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, int* NumEntriesPerRow, bool ignoreNonLocalEntries=false)
*/
CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_Var ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  int * NumEntriesPerRow, boolean ignoreNonLocalEntries );

/*! @brief Wrapper for 
   Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, int NumEntriesPerRow, bool ignoreNonLocalEntries=false)
*/
CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  int NumEntriesPerRow, boolean ignoreNonLocalEntries );

/*! @brief Wrapper for 
   Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, int* NumEntriesPerRow, bool ignoreNonLocalEntries=false)
*/
CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_WithColMap_Var ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int * NumEntriesPerRow, 
  boolean ignoreNonLocalEntries );

/*! @brief Wrapper for 
   Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, int NumEntriesPerRow, bool ignoreNonLocalEntries=false)
*/
CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, 
  boolean ignoreNonLocalEntries );

/*! @brief Wrapper for 
   Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph& Graph, bool ignoreNonLocalEntries=false)
*/
CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_FromGraph ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_CrsGraph_ID_t GraphID, 
  boolean ignoreNonLocalEntries );

/*! @brief Wrapper for 
   Epetra_FECrsMatrix::Epetra_FECrsMatrix(const Epetra_FECrsMatrix& src)
*/
CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Duplicate ( 
  CT_Epetra_FECrsMatrix_ID_t srcID );

/*@}*/

/*! @name Epetra_FECrsMatrix destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_FECrsMatrix::~Epetra_FECrsMatrix()
*/
void Epetra_FECrsMatrix_Destroy ( 
  CT_Epetra_FECrsMatrix_ID_t * selfID );

/*@}*/

/*! @name Epetra_FECrsMatrix member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
*/
int Epetra_FECrsMatrix_SumIntoGlobalValues ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
*/
int Epetra_FECrsMatrix_InsertGlobalValues ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::ReplaceGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
*/
int Epetra_FECrsMatrix_ReplaceGlobalValues ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const int* indices, const double* values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
*/
int Epetra_FECrsMatrix_SumIntoGlobalValues_Ftable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double * values, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const int* rows, int numCols, const int* cols, const double* values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
*/
int Epetra_FECrsMatrix_SumIntoGlobalValues_Ftable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double * values, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const int* indices, const double* const* values, int format=Epetra_FECrsMatrix::ROW_MAJOR)
*/
int Epetra_FECrsMatrix_SumIntoGlobalValues_Ctable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double* const * values, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const int* rows, int numCols, const int* cols, const double* const* values, int format=Epetra_FECrsMatrix::ROW_MAJOR)
*/
int Epetra_FECrsMatrix_SumIntoGlobalValues_Ctable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double* const * values, 
  int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::InsertGlobalValues(int numIndices, const int* indices, const double* values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
*/
int Epetra_FECrsMatrix_InsertGlobalValues_Ftable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double * values, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::InsertGlobalValues(int numRows, const int* rows, int numCols, const int* cols, const double* values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
*/
int Epetra_FECrsMatrix_InsertGlobalValues_Ftable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double * values, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::InsertGlobalValues(int numIndices, const int* indices, const double* const* values, int format=Epetra_FECrsMatrix::ROW_MAJOR)
*/
int Epetra_FECrsMatrix_InsertGlobalValues_Ctable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double* const * values, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::InsertGlobalValues(int numRows, const int* rows, int numCols, const int* cols, const double* const* values, int format=Epetra_FECrsMatrix::ROW_MAJOR)
*/
int Epetra_FECrsMatrix_InsertGlobalValues_Ctable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double* const * values, 
  int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const int* indices, const double* values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
*/
int Epetra_FECrsMatrix_ReplaceGlobalValues_Ftable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double * values, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const int* rows, int numCols, const int* cols, const double* values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
*/
int Epetra_FECrsMatrix_ReplaceGlobalValues_Ftable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double * values, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const int* indices, const double* const* values, int format=Epetra_FECrsMatrix::ROW_MAJOR)
*/
int Epetra_FECrsMatrix_ReplaceGlobalValues_Ctable_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numIndices, 
  const int * indices, const double* const * values, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const int* rows, int numCols, const int* cols, const double* const* values, int format=Epetra_FECrsMatrix::ROW_MAJOR)
*/
int Epetra_FECrsMatrix_ReplaceGlobalValues_Ctable ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, int numRows, const int * rows, 
  int numCols, const int * cols, const double* const * values, 
  int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::SumIntoGlobalValues(const Epetra_IntSerialDenseVector& indices, const Epetra_SerialDenseMatrix& values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
*/
int Epetra_FECrsMatrix_SumIntoGlobalValues_SubMatrix_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t indicesID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::SumIntoGlobalValues(const Epetra_IntSerialDenseVector& rows, const Epetra_IntSerialDenseVector& cols, const Epetra_SerialDenseMatrix& values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
*/
int Epetra_FECrsMatrix_SumIntoGlobalValues_SubMatrix ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t rowsID, 
  CT_Epetra_IntSerialDenseVector_ID_t colsID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::InsertGlobalValues(const Epetra_IntSerialDenseVector& indices, const Epetra_SerialDenseMatrix& values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
*/
int Epetra_FECrsMatrix_InsertGlobalValues_SubMatrix_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t indicesID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::InsertGlobalValues(const Epetra_IntSerialDenseVector& rows, const Epetra_IntSerialDenseVector& cols, const Epetra_SerialDenseMatrix& values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
*/
int Epetra_FECrsMatrix_InsertGlobalValues_SubMatrix ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t rowsID, 
  CT_Epetra_IntSerialDenseVector_ID_t colsID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::ReplaceGlobalValues(const Epetra_IntSerialDenseVector& indices, const Epetra_SerialDenseMatrix& values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
*/
int Epetra_FECrsMatrix_ReplaceGlobalValues_SubMatrix_Square ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t indicesID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::ReplaceGlobalValues(const Epetra_IntSerialDenseVector& rows, const Epetra_IntSerialDenseVector& cols, const Epetra_SerialDenseMatrix& values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
*/
int Epetra_FECrsMatrix_ReplaceGlobalValues_SubMatrix ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t rowsID, 
  CT_Epetra_IntSerialDenseVector_ID_t colsID, 
  CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::GlobalAssemble(bool callFillComplete=true)
*/
int Epetra_FECrsMatrix_GlobalAssemble ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, boolean callFillComplete );

/*! @brief Wrapper for 
   int Epetra_FECrsMatrix::GlobalAssemble(const Epetra_Map& domain_map, const Epetra_Map& range_map, bool callFillComplete=true)
*/
int Epetra_FECrsMatrix_GlobalAssemble_WithMaps ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_Map_ID_t domain_mapID, CT_Epetra_Map_ID_t range_mapID, 
  boolean callFillComplete );

/*! @brief Wrapper for 
   void Epetra_FECrsMatrix::setIgnoreNonLocalEntries(bool flag)
*/
void Epetra_FECrsMatrix_setIgnoreNonLocalEntries ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, boolean flag );

/*@}*/

/*! @name Epetra_FECrsMatrix operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_FECrsMatrix& Epetra_FECrsMatrix::operator=(const Epetra_FECrsMatrix& src)
*/
void Epetra_FECrsMatrix_Assign ( 
  CT_Epetra_FECrsMatrix_ID_t selfID, 
  CT_Epetra_FECrsMatrix_ID_t srcID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_FECRSMATRIX_H */

