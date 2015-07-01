#ifndef CEPETRA_CRSGRAPH_H
#define CEPETRA_CRSGRAPH_H

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


/*! @file CEpetra_CrsGraph.h
 * @brief Wrappers for Epetra_CrsGraph */

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
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_CrsGraph_Generalize ( 
  CT_Epetra_CrsGraph_ID_t id );

/*@}*/

/*! @name Epetra_CrsGraph constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const int* NumIndicesPerRow, bool StaticProfile = false)
*/
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  const int * NumIndicesPerRow, boolean StaticProfile );

/*! @brief Wrapper for 
   Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int NumIndicesPerRow, bool StaticProfile = false)
*/
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  int NumIndicesPerRow, boolean StaticProfile );

/*! @brief Wrapper for 
   Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const Epetra_BlockMap& ColMap, const int* NumIndicesPerRow, bool StaticProfile = false)
*/
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, const int * NumIndicesPerRow, 
  boolean StaticProfile );

/*! @brief Wrapper for 
   Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const Epetra_BlockMap& ColMap, int NumIndicesPerRow, bool StaticProfile = false)
*/
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_With_ColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, int NumIndicesPerRow, 
  boolean StaticProfile );

/*! @brief Wrapper for 
   Epetra_CrsGraph::Epetra_CrsGraph(const Epetra_CrsGraph& Graph)
*/
CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Duplicate ( 
  CT_Epetra_CrsGraph_ID_t GraphID );

/*@}*/

/*! @name Epetra_CrsGraph destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_CrsGraph::~Epetra_CrsGraph()
*/
void Epetra_CrsGraph_Destroy ( CT_Epetra_CrsGraph_ID_t * selfID );

/*@}*/

/*! @name Epetra_CrsGraph member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_CrsGraph::InsertGlobalIndices(int GlobalRow, int NumIndices, int* Indices)
*/
int Epetra_CrsGraph_InsertGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::RemoveGlobalIndices(int GlobalRow, int NumIndices, int* Indices)
*/
int Epetra_CrsGraph_RemoveGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::RemoveGlobalIndices(int Row)
*/
int Epetra_CrsGraph_RemoveGlobalIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::InsertMyIndices(int LocalRow, int NumIndices, int* Indices)
*/
int Epetra_CrsGraph_InsertMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::RemoveMyIndices(int LocalRow, int NumIndices, int* Indices)
*/
int Epetra_CrsGraph_RemoveMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::RemoveMyIndices(int Row)
*/
int Epetra_CrsGraph_RemoveMyIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::FillComplete()
*/
int Epetra_CrsGraph_FillComplete ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::FillComplete(const Epetra_BlockMap& DomainMap, const Epetra_BlockMap& RangeMap)
*/
int Epetra_CrsGraph_FillComplete_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, 
  CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::OptimizeStorage()
*/
int Epetra_CrsGraph_OptimizeStorage ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::ExtractGlobalRowCopy(int GlobalRow, int LenOfIndices, int& NumIndices, int* Indices) const
*/
int Epetra_CrsGraph_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int LenOfIndices, 
  int * NumIndices, int * Indices );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::ExtractMyRowCopy(int LocalRow, int LenOfIndices, int& NumIndices, int* Indices) const
*/
int Epetra_CrsGraph_ExtractMyRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int LenOfIndices, 
  int * NumIndices, int * Indices );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::ExtractGlobalRowView(int GlobalRow, int& NumIndices, int*& Indices) const
*/
int Epetra_CrsGraph_ExtractGlobalRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int * NumIndices, 
  int ** Indices );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::ExtractMyRowView(int LocalRow, int& NumIndices, int*& Indices) const
*/
int Epetra_CrsGraph_ExtractMyRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int * NumIndices, 
  int ** Indices );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::Filled() const
*/
boolean Epetra_CrsGraph_Filled ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::StorageOptimized() const
*/
boolean Epetra_CrsGraph_StorageOptimized ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::IndicesAreGlobal() const
*/
boolean Epetra_CrsGraph_IndicesAreGlobal ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::IndicesAreLocal() const
*/
boolean Epetra_CrsGraph_IndicesAreLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::LowerTriangular() const
*/
boolean Epetra_CrsGraph_LowerTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::UpperTriangular() const
*/
boolean Epetra_CrsGraph_UpperTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::NoDiagonal() const
*/
boolean Epetra_CrsGraph_NoDiagonal ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::MyGlobalRow(int GID) const
*/
boolean Epetra_CrsGraph_MyGlobalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GID );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::HaveColMap() const
*/
boolean Epetra_CrsGraph_HaveColMap ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumMyRows() const
*/
int Epetra_CrsGraph_NumMyRows ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumGlobalRows() const
*/
int Epetra_CrsGraph_NumGlobalRows ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumMyCols() const
*/
int Epetra_CrsGraph_NumMyCols ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumGlobalCols() const
*/
int Epetra_CrsGraph_NumGlobalCols ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumGlobalNonzeros() const
*/
int Epetra_CrsGraph_NumGlobalNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumGlobalDiagonals() const
*/
int Epetra_CrsGraph_NumGlobalDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumMyDiagonals() const
*/
int Epetra_CrsGraph_NumMyDiagonals ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumMyBlockRows() const
*/
int Epetra_CrsGraph_NumMyBlockRows ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumGlobalBlockRows() const
*/
int Epetra_CrsGraph_NumGlobalBlockRows ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumMyBlockCols() const
*/
int Epetra_CrsGraph_NumMyBlockCols ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumGlobalBlockCols() const
*/
int Epetra_CrsGraph_NumGlobalBlockCols ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumMyBlockDiagonals() const
*/
int Epetra_CrsGraph_NumMyBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumGlobalBlockDiagonals() const
*/
int Epetra_CrsGraph_NumGlobalBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumGlobalEntries() const
*/
int Epetra_CrsGraph_NumGlobalEntries ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumMyEntries() const
*/
int Epetra_CrsGraph_NumMyEntries ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::MaxRowDim() const
*/
int Epetra_CrsGraph_MaxRowDim ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::GlobalMaxRowDim() const
*/
int Epetra_CrsGraph_GlobalMaxRowDim ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::MaxColDim() const
*/
int Epetra_CrsGraph_MaxColDim ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::GlobalMaxColDim() const
*/
int Epetra_CrsGraph_GlobalMaxColDim ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumMyNonzeros() const
*/
int Epetra_CrsGraph_NumMyNonzeros ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumGlobalIndices(int Row) const
*/
int Epetra_CrsGraph_NumGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumAllocatedGlobalIndices(int Row) const
*/
int Epetra_CrsGraph_NumAllocatedGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::MaxNumIndices() const
*/
int Epetra_CrsGraph_MaxNumIndices ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::GlobalMaxNumIndices() const
*/
int Epetra_CrsGraph_GlobalMaxNumIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::MaxNumNonzeros() const
*/
int Epetra_CrsGraph_MaxNumNonzeros ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::GlobalMaxNumNonzeros() const
*/
int Epetra_CrsGraph_GlobalMaxNumNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumMyIndices(int Row) const
*/
int Epetra_CrsGraph_NumMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::NumAllocatedMyIndices(int Row) const
*/
int Epetra_CrsGraph_NumAllocatedMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::IndexBase() const
*/
int Epetra_CrsGraph_IndexBase ( CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_BlockMap& Epetra_CrsGraph::RowMap() const
*/
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::ReplaceRowMap(const Epetra_BlockMap& newmap)
*/
int Epetra_CrsGraph_ReplaceRowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::ReplaceColMap(const Epetra_BlockMap& newmap)
*/
int Epetra_CrsGraph_ReplaceColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID );

/*! @brief Wrapper for 
   const Epetra_BlockMap& Epetra_CrsGraph::ColMap() const
*/
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_BlockMap& Epetra_CrsGraph::DomainMap() const
*/
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_DomainMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_BlockMap& Epetra_CrsGraph::RangeMap() const
*/
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RangeMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Import* Epetra_CrsGraph::Importer() const
*/
CT_Epetra_Import_ID_t Epetra_CrsGraph_Importer ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Export* Epetra_CrsGraph::Exporter() const
*/
CT_Epetra_Export_ID_t Epetra_CrsGraph_Exporter ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Comm& Epetra_CrsGraph::Comm() const
*/
CT_Epetra_Comm_ID_t Epetra_CrsGraph_Comm ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::LRID(int GRID_in) const
*/
int Epetra_CrsGraph_LRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::GRID(int LRID_in) const
*/
int Epetra_CrsGraph_GRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::LCID(int GCID_in) const
*/
int Epetra_CrsGraph_LCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::GCID(int LCID_in) const
*/
int Epetra_CrsGraph_GCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::MyGRID(int GRID_in) const
*/
boolean Epetra_CrsGraph_MyGRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::MyLRID(int LRID_in) const
*/
boolean Epetra_CrsGraph_MyLRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::MyGCID(int GCID_in) const
*/
boolean Epetra_CrsGraph_MyGCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in );

/*! @brief Wrapper for 
   bool Epetra_CrsGraph::MyLCID(int LCID_in) const
*/
boolean Epetra_CrsGraph_MyLCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in );

/*! @brief Wrapper for 
   const Epetra_BlockMap& Epetra_CrsGraph::ImportMap() const
*/
CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ImportMap ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::TransformToLocal()
*/
int Epetra_CrsGraph_TransformToLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_CrsGraph::TransformToLocal(const Epetra_BlockMap* DomainMap, const Epetra_BlockMap* RangeMap)
*/
int Epetra_CrsGraph_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, 
  CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID );

/*@}*/

/*! @name Epetra_CrsGraph operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   inline int* Epetra_CrsGraph::operator[]( int Loc ) const
*/
int * Epetra_CrsGraph_getRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Loc );

/*! @brief Wrapper for 
   Epetra_CrsGraph& Epetra_CrsGraph::operator= (const Epetra_CrsGraph& Source)
*/
void Epetra_CrsGraph_Assign ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_CrsGraph_ID_t SourceID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_CRSGRAPH_H */

