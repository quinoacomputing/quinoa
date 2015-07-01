#ifndef CEPETRA_CRSMATRIX_H
#define CEPETRA_CRSMATRIX_H

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


/*! @file CEpetra_CrsMatrix.h
 * @brief Wrappers for Epetra_CrsMatrix */

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
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_CrsMatrix_Generalize ( 
  CT_Epetra_CrsMatrix_ID_t id );

/*@}*/

/*! @name Epetra_CrsMatrix constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const int* NumEntriesPerRow, bool StaticProfile = false)
*/
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  const int * NumEntriesPerRow, boolean StaticProfile );

/*! @brief Wrapper for 
   Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, int NumEntriesPerRow, bool StaticProfile = false)
*/
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  int NumEntriesPerRow, boolean StaticProfile );

/*! @brief Wrapper for 
   Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, const int* NumEntriesPerRow, bool StaticProfile = false)
*/
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, const int * NumEntriesPerRow, 
  boolean StaticProfile );

/*! @brief Wrapper for 
   Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, int NumEntriesPerRow, bool StaticProfile = false)
*/
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_Map_ID_t RowMapID, 
  CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, 
  boolean StaticProfile );

/*! @brief Wrapper for 
   Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph& Graph)
*/
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_FromGraph ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_CrsGraph_ID_t GraphID );

/*! @brief Wrapper for 
   Epetra_CrsMatrix::Epetra_CrsMatrix(const Epetra_CrsMatrix& Matrix)
*/
CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Duplicate ( 
  CT_Epetra_CrsMatrix_ID_t MatrixID );

/*@}*/

/*! @name Epetra_CrsMatrix destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_CrsMatrix::~Epetra_CrsMatrix()
*/
void Epetra_CrsMatrix_Destroy ( CT_Epetra_CrsMatrix_ID_t * selfID );

/*@}*/

/*! @name Epetra_CrsMatrix member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::PutScalar(double ScalarConstant)
*/
int Epetra_CrsMatrix_PutScalar ( 
  CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::Scale(double ScalarConstant)
*/
int Epetra_CrsMatrix_Scale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant );

/*! @brief Wrapper for 
   virtual int Epetra_CrsMatrix::InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
*/
int Epetra_CrsMatrix_InsertGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   virtual int Epetra_CrsMatrix::ReplaceGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
*/
int Epetra_CrsMatrix_ReplaceGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   virtual int Epetra_CrsMatrix::SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
*/
int Epetra_CrsMatrix_SumIntoGlobalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int NumEntries, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::InsertMyValues(int MyRow, int NumEntries, double* Values, int* Indices)
*/
int Epetra_CrsMatrix_InsertMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ReplaceMyValues(int MyRow, int NumEntries, double* Values, int* Indices)
*/
int Epetra_CrsMatrix_ReplaceMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::SumIntoMyValues(int MyRow, int NumEntries, double* Values, int* Indices)
*/
int Epetra_CrsMatrix_SumIntoMyValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int NumEntries, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ReplaceDiagonalValues(const Epetra_Vector& Diagonal)
*/
int Epetra_CrsMatrix_ReplaceDiagonalValues ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::FillComplete(bool OptimizeDataStorage = true)
*/
int Epetra_CrsMatrix_FillComplete ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean OptimizeDataStorage );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::FillComplete(const Epetra_Map& DomainMap, const Epetra_Map& RangeMap, bool OptimizeDataStorage = true)
*/
int Epetra_CrsMatrix_FillComplete_UsingMaps ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Map_ID_t DomainMapID, 
  CT_Epetra_Map_ID_t RangeMapID, boolean OptimizeDataStorage );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::OptimizeStorage()
*/
int Epetra_CrsMatrix_OptimizeStorage ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::MakeDataContiguous()
*/
int Epetra_CrsMatrix_MakeDataContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ExtractGlobalRowCopy(int GlobalRow, int Length, int& NumEntries, double* Values, int* Indices) const
*/
int Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, double* Values, int* Indices) const
*/
int Epetra_CrsMatrix_ExtractMyRowCopy_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ExtractGlobalRowCopy(int GlobalRow, int Length, int& NumEntries, double* Values) const
*/
int Epetra_CrsMatrix_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int Length, 
  int * NumEntries, double * Values );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, double* Values) const
*/
int Epetra_CrsMatrix_ExtractMyRowCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ExtractDiagonalCopy(Epetra_Vector& Diagonal) const
*/
int Epetra_CrsMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ExtractGlobalRowView(int GlobalRow, int& NumEntries, double*& Values, int*& Indices) const
*/
int Epetra_CrsMatrix_ExtractGlobalRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values, int ** Indices );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ExtractMyRowView(int MyRow, int& NumEntries, double*& Values, int*& Indices) const
*/
int Epetra_CrsMatrix_ExtractMyRowView_WithIndices ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values, int ** Indices );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ExtractGlobalRowView(int GlobalRow, int& NumEntries, double*& Values) const
*/
int Epetra_CrsMatrix_ExtractGlobalRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, int * NumEntries, 
  double ** Values );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ExtractMyRowView(int MyRow, int& NumEntries, double*& Values) const
*/
int Epetra_CrsMatrix_ExtractMyRowView ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries, 
  double ** Values );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::Multiply(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const
*/
int Epetra_CrsMatrix_Multiply_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::Multiply1(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const
*/
int Epetra_CrsMatrix_Multiply1_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
*/
int Epetra_CrsMatrix_Multiply_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::Multiply1(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
*/
int Epetra_CrsMatrix_Multiply1_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector& y) const
*/
int Epetra_CrsMatrix_Solve_Vector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_Vector_ID_t xID, 
  CT_Epetra_Vector_ID_t yID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
*/
int Epetra_CrsMatrix_Solve_MultiVector ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::InvRowSums(Epetra_Vector& x) const
*/
int Epetra_CrsMatrix_InvRowSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::InvRowMaxs(Epetra_Vector& x) const
*/
int Epetra_CrsMatrix_InvRowMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::LeftScale(const Epetra_Vector& x)
*/
int Epetra_CrsMatrix_LeftScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::InvColSums(Epetra_Vector& x) const
*/
int Epetra_CrsMatrix_InvColSums ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::InvColMaxs(Epetra_Vector& x) const
*/
int Epetra_CrsMatrix_InvColMaxs ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::RightScale(const Epetra_Vector& x)
*/
int Epetra_CrsMatrix_RightScale ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::Filled() const
*/
boolean Epetra_CrsMatrix_Filled ( CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::StorageOptimized() const
*/
boolean Epetra_CrsMatrix_StorageOptimized ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::IndicesAreGlobal() const
*/
boolean Epetra_CrsMatrix_IndicesAreGlobal ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::IndicesAreLocal() const
*/
boolean Epetra_CrsMatrix_IndicesAreLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::IndicesAreContiguous() const
*/
boolean Epetra_CrsMatrix_IndicesAreContiguous ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::LowerTriangular() const
*/
boolean Epetra_CrsMatrix_LowerTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::UpperTriangular() const
*/
boolean Epetra_CrsMatrix_UpperTriangular ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::NoDiagonal() const
*/
boolean Epetra_CrsMatrix_NoDiagonal ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   double Epetra_CrsMatrix::NormInf() const
*/
double Epetra_CrsMatrix_NormInf ( CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   double Epetra_CrsMatrix::NormOne() const
*/
double Epetra_CrsMatrix_NormOne ( CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   double Epetra_CrsMatrix::NormFrobenius() const
*/
double Epetra_CrsMatrix_NormFrobenius ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumGlobalNonzeros() const
*/
int Epetra_CrsMatrix_NumGlobalNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumGlobalRows() const
*/
int Epetra_CrsMatrix_NumGlobalRows ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumGlobalCols() const
*/
int Epetra_CrsMatrix_NumGlobalCols ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumGlobalDiagonals() const
*/
int Epetra_CrsMatrix_NumGlobalDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumMyNonzeros() const
*/
int Epetra_CrsMatrix_NumMyNonzeros ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumMyRows() const
*/
int Epetra_CrsMatrix_NumMyRows ( CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumMyCols() const
*/
int Epetra_CrsMatrix_NumMyCols ( CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumMyDiagonals() const
*/
int Epetra_CrsMatrix_NumMyDiagonals ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumGlobalEntries(int Row) const
*/
int Epetra_CrsMatrix_NumGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumAllocatedGlobalEntries(int Row) const
*/
int Epetra_CrsMatrix_NumAllocatedGlobalEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::MaxNumEntries() const
*/
int Epetra_CrsMatrix_MaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::GlobalMaxNumEntries() const
*/
int Epetra_CrsMatrix_GlobalMaxNumEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumMyEntries(int Row) const
*/
int Epetra_CrsMatrix_NumMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumAllocatedMyEntries(int Row) const
*/
int Epetra_CrsMatrix_NumAllocatedMyEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Row );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::IndexBase() const
*/
int Epetra_CrsMatrix_IndexBase ( CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::StaticGraph()
*/
boolean Epetra_CrsMatrix_StaticGraph ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_CrsGraph& Epetra_CrsMatrix::Graph() const
*/
CT_Epetra_CrsGraph_ID_t Epetra_CrsMatrix_Graph ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Map& Epetra_CrsMatrix::RowMap() const
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ReplaceRowMap(const Epetra_BlockMap& newmap)
*/
int Epetra_CrsMatrix_ReplaceRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::HaveColMap() const
*/
boolean Epetra_CrsMatrix_HaveColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ReplaceColMap(const Epetra_BlockMap& newmap)
*/
int Epetra_CrsMatrix_ReplaceColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID );

/*! @brief Wrapper for 
   const Epetra_Map& Epetra_CrsMatrix::ColMap() const
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_ColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Map& Epetra_CrsMatrix::DomainMap() const
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_DomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Map& Epetra_CrsMatrix::RangeMap() const
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_RangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Import* Epetra_CrsMatrix::Importer() const
*/
CT_Epetra_Import_ID_t Epetra_CrsMatrix_Importer ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Export* Epetra_CrsMatrix::Exporter() const
*/
CT_Epetra_Export_ID_t Epetra_CrsMatrix_Exporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Comm& Epetra_CrsMatrix::Comm() const
*/
CT_Epetra_Comm_ID_t Epetra_CrsMatrix_Comm ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::LRID( int GRID_in) const
*/
int Epetra_CrsMatrix_LRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::GRID( int LRID_in) const
*/
int Epetra_CrsMatrix_GRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::LCID( int GCID_in) const
*/
int Epetra_CrsMatrix_LCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::GCID( int LCID_in) const
*/
int Epetra_CrsMatrix_GCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::MyGRID(int GRID_in) const
*/
boolean Epetra_CrsMatrix_MyGRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::MyLRID(int LRID_in) const
*/
boolean Epetra_CrsMatrix_MyLRID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::MyGCID(int GCID_in) const
*/
boolean Epetra_CrsMatrix_MyGCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::MyLCID(int LCID_in) const
*/
boolean Epetra_CrsMatrix_MyLCID ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::MyGlobalRow(int GID) const
*/
boolean Epetra_CrsMatrix_MyGlobalRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int GID );

/*! @brief Wrapper for 
   const char* Epetra_CrsMatrix::Label() const
*/
const char * Epetra_CrsMatrix_Label ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::SetUseTranspose(bool UseTranspose_in)
*/
int Epetra_CrsMatrix_SetUseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID, boolean UseTranspose_in );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
*/
int Epetra_CrsMatrix_Apply ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
*/
int Epetra_CrsMatrix_ApplyInverse ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::HasNormInf() const
*/
boolean Epetra_CrsMatrix_HasNormInf ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsMatrix::UseTranspose() const
*/
boolean Epetra_CrsMatrix_UseTranspose ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Map& Epetra_CrsMatrix::OperatorDomainMap() const
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorDomainMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Map& Epetra_CrsMatrix::OperatorRangeMap() const
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorRangeMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::NumMyRowEntries(int MyRow, int& NumEntries) const
*/
int Epetra_CrsMatrix_NumMyRowEntries ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, int * NumEntries );

/*! @brief Wrapper for 
   const Epetra_Map& Epetra_CrsMatrix::RowMatrixRowMap() const
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixRowMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Map& Epetra_CrsMatrix::RowMatrixColMap() const
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixColMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Import* Epetra_CrsMatrix::RowMatrixImporter() const
*/
CT_Epetra_Import_ID_t Epetra_CrsMatrix_RowMatrixImporter ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Map& Epetra_CrsMatrix::ImportMap() const
*/
CT_Epetra_Map_ID_t Epetra_CrsMatrix_ImportMap ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::TransformToLocal()
*/
int Epetra_CrsMatrix_TransformToLocal ( 
  CT_Epetra_CrsMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsMatrix::TransformToLocal(const Epetra_Map* DomainMap, const Epetra_Map* RangeMap)
*/
int Epetra_CrsMatrix_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Map_ID_t DomainMapID, 
  CT_Epetra_Map_ID_t RangeMapID );

/*@}*/

/*! @name Epetra_CrsMatrix operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_CrsMatrix& Epetra_CrsMatrix::operator=(const Epetra_CrsMatrix& src)
*/
void Epetra_CrsMatrix_Assign ( 
  CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_CrsMatrix_ID_t srcID );

/*! @brief Wrapper for 
   inline double* Epetra_CrsMatrix::operator[] (int Loc) const
*/
double * Epetra_CrsMatrix_getRow ( 
  CT_Epetra_CrsMatrix_ID_t selfID, int Loc );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_CRSMATRIX_H */

