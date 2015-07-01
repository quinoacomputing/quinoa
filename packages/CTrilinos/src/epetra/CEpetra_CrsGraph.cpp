
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

#include "CTrilinos_enums.h"
#include "CEpetra_CrsGraph.h"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "Epetra_CrsGraph.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_CrsGraph */
Table<Epetra_CrsGraph>& tableOfCrsGraphs()
{
    static Table<Epetra_CrsGraph> loc_tableOfCrsGraphs(CT_Epetra_CrsGraph_ID);
    return loc_tableOfCrsGraphs;
}


} // namespace


//
// Definitions from CEpetra_CrsGraph.h
//


extern "C" {


CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_CrsGraph_Generalize ( 
  CT_Epetra_CrsGraph_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(id);
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  const int * NumIndicesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_BlockMap> RowMap = 
        CEpetra::getConstBlockMap(RowMapID);
    return CEpetra::storeNewCrsGraph(new Epetra_CrsGraph(
        (Epetra_DataAccess) CV, *RowMap, NumIndicesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  int NumIndicesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_BlockMap> RowMap = 
        CEpetra::getConstBlockMap(RowMapID);
    return CEpetra::storeNewCrsGraph(new Epetra_CrsGraph(
        (Epetra_DataAccess) CV, *RowMap, NumIndicesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, const int * NumIndicesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_BlockMap> RowMap = 
        CEpetra::getConstBlockMap(RowMapID);
    const Teuchos::RCP<const Epetra_BlockMap> ColMap = 
        CEpetra::getConstBlockMap(ColMapID);
    return CEpetra::storeNewCrsGraph(new Epetra_CrsGraph(
        (Epetra_DataAccess) CV, *RowMap, *ColMap, NumIndicesPerRow, ((
        StaticProfile) != FALSE ? true : false)));
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_With_ColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, int NumIndicesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_BlockMap> RowMap = 
        CEpetra::getConstBlockMap(RowMapID);
    const Teuchos::RCP<const Epetra_BlockMap> ColMap = 
        CEpetra::getConstBlockMap(ColMapID);
    return CEpetra::storeNewCrsGraph(new Epetra_CrsGraph(
        (Epetra_DataAccess) CV, *RowMap, *ColMap, NumIndicesPerRow, ((
        StaticProfile) != FALSE ? true : false)));
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Duplicate ( 
  CT_Epetra_CrsGraph_ID_t GraphID )
{
    const Teuchos::RCP<const Epetra_CrsGraph> Graph = CEpetra::getConstCrsGraph(
        GraphID);
    return CEpetra::storeNewCrsGraph(new Epetra_CrsGraph(*Graph));
}

void Epetra_CrsGraph_Destroy ( CT_Epetra_CrsGraph_ID_t * selfID )
{
    CEpetra::removeCrsGraph(selfID);
}

int Epetra_CrsGraph_InsertGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->InsertGlobalIndices(GlobalRow, 
        NumIndices, Indices);
}

int Epetra_CrsGraph_RemoveGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->RemoveGlobalIndices(GlobalRow, 
        NumIndices, Indices);
}

int Epetra_CrsGraph_RemoveGlobalIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getCrsGraph(selfID)->RemoveGlobalIndices(Row);
}

int Epetra_CrsGraph_InsertMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->InsertMyIndices(LocalRow, NumIndices, 
        Indices);
}

int Epetra_CrsGraph_RemoveMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->RemoveMyIndices(LocalRow, NumIndices, 
        Indices);
}

int Epetra_CrsGraph_RemoveMyIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getCrsGraph(selfID)->RemoveMyIndices(Row);
}

int Epetra_CrsGraph_FillComplete ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->FillComplete();
}

int Epetra_CrsGraph_FillComplete_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, 
  CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> DomainMap = 
        CEpetra::getConstBlockMap(DomainMapID);
    const Teuchos::RCP<const Epetra_BlockMap> RangeMap = 
        CEpetra::getConstBlockMap(RangeMapID);
    return CEpetra::getCrsGraph(selfID)->FillComplete(*DomainMap, *RangeMap);
}

int Epetra_CrsGraph_OptimizeStorage ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->OptimizeStorage();
}

int Epetra_CrsGraph_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int LenOfIndices, 
  int * NumIndices, int * Indices )
{
    return CEpetra::getConstCrsGraph(selfID)->ExtractGlobalRowCopy(GlobalRow, 
        LenOfIndices, *NumIndices, Indices);
}

int Epetra_CrsGraph_ExtractMyRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int LenOfIndices, 
  int * NumIndices, int * Indices )
{
    return CEpetra::getConstCrsGraph(selfID)->ExtractMyRowCopy(LocalRow, 
        LenOfIndices, *NumIndices, Indices);
}

int Epetra_CrsGraph_ExtractGlobalRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int * NumIndices, 
  int ** Indices )
{
    return CEpetra::getConstCrsGraph(selfID)->ExtractGlobalRowView(GlobalRow, 
        *NumIndices, *Indices);
}

int Epetra_CrsGraph_ExtractMyRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int * NumIndices, 
  int ** Indices )
{
    return CEpetra::getConstCrsGraph(selfID)->ExtractMyRowView(LocalRow, 
        *NumIndices, *Indices);
}

boolean Epetra_CrsGraph_Filled ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(selfID)->Filled()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_StorageOptimized ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(
        selfID)->StorageOptimized()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_IndicesAreGlobal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(
        selfID)->IndicesAreGlobal()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_IndicesAreLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(
        selfID)->IndicesAreLocal()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_LowerTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(
        selfID)->LowerTriangular()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_UpperTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(
        selfID)->UpperTriangular()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_NoDiagonal ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(selfID)->NoDiagonal()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_MyGlobalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GID )
{
    return ((CEpetra::getConstCrsGraph(selfID)->MyGlobalRow(
        GID)) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_HaveColMap ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(selfID)->HaveColMap()) ? TRUE : FALSE);
}

int Epetra_CrsGraph_NumMyRows ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyRows();
}

int Epetra_CrsGraph_NumGlobalRows ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalRows();
}

int Epetra_CrsGraph_NumMyCols ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyCols();
}

int Epetra_CrsGraph_NumGlobalCols ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalCols();
}

int Epetra_CrsGraph_NumGlobalNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalNonzeros();
}

int Epetra_CrsGraph_NumGlobalDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalDiagonals();
}

int Epetra_CrsGraph_NumMyDiagonals ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyDiagonals();
}

int Epetra_CrsGraph_NumMyBlockRows ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyBlockRows();
}

int Epetra_CrsGraph_NumGlobalBlockRows ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalBlockRows();
}

int Epetra_CrsGraph_NumMyBlockCols ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyBlockCols();
}

int Epetra_CrsGraph_NumGlobalBlockCols ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalBlockCols();
}

int Epetra_CrsGraph_NumMyBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyBlockDiagonals();
}

int Epetra_CrsGraph_NumGlobalBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalBlockDiagonals();
}

int Epetra_CrsGraph_NumGlobalEntries ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalEntries();
}

int Epetra_CrsGraph_NumMyEntries ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyEntries();
}

int Epetra_CrsGraph_MaxRowDim ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->MaxRowDim();
}

int Epetra_CrsGraph_GlobalMaxRowDim ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->GlobalMaxRowDim();
}

int Epetra_CrsGraph_MaxColDim ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->MaxColDim();
}

int Epetra_CrsGraph_GlobalMaxColDim ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->GlobalMaxColDim();
}

int Epetra_CrsGraph_NumMyNonzeros ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyNonzeros();
}

int Epetra_CrsGraph_NumGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalIndices(Row);
}

int Epetra_CrsGraph_NumAllocatedGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsGraph(selfID)->NumAllocatedGlobalIndices(Row);
}

int Epetra_CrsGraph_MaxNumIndices ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->MaxNumIndices();
}

int Epetra_CrsGraph_GlobalMaxNumIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->GlobalMaxNumIndices();
}

int Epetra_CrsGraph_MaxNumNonzeros ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->MaxNumNonzeros();
}

int Epetra_CrsGraph_GlobalMaxNumNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->GlobalMaxNumNonzeros();
}

int Epetra_CrsGraph_NumMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyIndices(Row);
}

int Epetra_CrsGraph_NumAllocatedMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsGraph(selfID)->NumAllocatedMyIndices(Row);
}

int Epetra_CrsGraph_IndexBase ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->IndexBase();
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstCrsGraph(
        selfID)->RowMap() ));
}

int Epetra_CrsGraph_ReplaceRowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> newmap = 
        CEpetra::getConstBlockMap(newmapID);
    return CEpetra::getCrsGraph(selfID)->ReplaceRowMap(*newmap);
}

int Epetra_CrsGraph_ReplaceColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> newmap = 
        CEpetra::getConstBlockMap(newmapID);
    return CEpetra::getCrsGraph(selfID)->ReplaceColMap(*newmap);
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstCrsGraph(
        selfID)->ColMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_DomainMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstCrsGraph(
        selfID)->DomainMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RangeMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstCrsGraph(
        selfID)->RangeMap() ));
}

CT_Epetra_Import_ID_t Epetra_CrsGraph_Importer ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstImport(CEpetra::getConstCrsGraph(
        selfID)->Importer());
}

CT_Epetra_Export_ID_t Epetra_CrsGraph_Exporter ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstExport(CEpetra::getConstCrsGraph(
        selfID)->Exporter());
}

CT_Epetra_Comm_ID_t Epetra_CrsGraph_Comm ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstComm(&( CEpetra::getConstCrsGraph(
        selfID)->Comm() ));
}

int Epetra_CrsGraph_LRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in )
{
    return CEpetra::getConstCrsGraph(selfID)->LRID(GRID_in);
}

int Epetra_CrsGraph_GRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in )
{
    return CEpetra::getConstCrsGraph(selfID)->GRID(LRID_in);
}

int Epetra_CrsGraph_LCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in )
{
    return CEpetra::getConstCrsGraph(selfID)->LCID(GCID_in);
}

int Epetra_CrsGraph_GCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in )
{
    return CEpetra::getConstCrsGraph(selfID)->GCID(LCID_in);
}

boolean Epetra_CrsGraph_MyGRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in )
{
    return ((CEpetra::getConstCrsGraph(selfID)->MyGRID(
        GRID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_MyLRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in )
{
    return ((CEpetra::getConstCrsGraph(selfID)->MyLRID(
        LRID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_MyGCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in )
{
    return ((CEpetra::getConstCrsGraph(selfID)->MyGCID(
        GCID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_MyLCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in )
{
    return ((CEpetra::getConstCrsGraph(selfID)->MyLCID(
        LCID_in)) ? TRUE : FALSE);
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ImportMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstCrsGraph(
        selfID)->ImportMap() ));
}

int Epetra_CrsGraph_TransformToLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->TransformToLocal();
}

int Epetra_CrsGraph_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, 
  CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> DomainMap = 
        CEpetra::getConstBlockMap(DomainMapID);
    const Teuchos::RCP<const Epetra_BlockMap> RangeMap = 
        CEpetra::getConstBlockMap(RangeMapID);
    return CEpetra::getCrsGraph(selfID)->TransformToLocal(
        DomainMap.getRawPtr(), RangeMap.getRawPtr());
}

int * Epetra_CrsGraph_getRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Loc )
{
    const Epetra_CrsGraph& self = *( CEpetra::getConstCrsGraph(selfID) );

    return self[Loc];
}

void Epetra_CrsGraph_Assign ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_CrsGraph_ID_t SourceID )
{
    Epetra_CrsGraph& self = *( CEpetra::getCrsGraph(selfID) );

    const Teuchos::RCP<const Epetra_CrsGraph> Source = 
        CEpetra::getConstCrsGraph(SourceID);
    self = *Source;
}


} // extern "C"


//
// Definitions from CEpetra_CrsGraph_Cpp.hpp
//


/* get Epetra_CrsGraph from non-const table using CT_Epetra_CrsGraph_ID */
const Teuchos::RCP<Epetra_CrsGraph>
CEpetra::getCrsGraph( CT_Epetra_CrsGraph_ID_t id )
{
    if (tableOfCrsGraphs().isType(id.table))
        return tableOfCrsGraphs().get<Epetra_CrsGraph>(
        CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_CrsGraph>(
        CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(id));
}

/* get Epetra_CrsGraph from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_CrsGraph>
CEpetra::getCrsGraph( CTrilinos_Universal_ID_t id )
{
    if (tableOfCrsGraphs().isType(id.table))
        return tableOfCrsGraphs().get<Epetra_CrsGraph>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_CrsGraph>(id);
}

/* get const Epetra_CrsGraph from either the const or non-const table
 * using CT_Epetra_CrsGraph_ID */
const Teuchos::RCP<const Epetra_CrsGraph>
CEpetra::getConstCrsGraph( CT_Epetra_CrsGraph_ID_t id )
{
    if (tableOfCrsGraphs().isType(id.table))
        return tableOfCrsGraphs().getConst<Epetra_CrsGraph>(
        CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_CrsGraph>(
        CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(id));
}

/* get const Epetra_CrsGraph from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_CrsGraph>
CEpetra::getConstCrsGraph( CTrilinos_Universal_ID_t id )
{
    if (tableOfCrsGraphs().isType(id.table))
        return tableOfCrsGraphs().getConst<Epetra_CrsGraph>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_CrsGraph>(id);
}

/* store Epetra_CrsGraph (owned) in non-const table */
CT_Epetra_CrsGraph_ID_t
CEpetra::storeNewCrsGraph( Epetra_CrsGraph *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(
        tableOfCrsGraphs().store<Epetra_CrsGraph>(pobj, true));
}

/* store Epetra_CrsGraph in non-const table */
CT_Epetra_CrsGraph_ID_t
CEpetra::storeCrsGraph( Epetra_CrsGraph *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(
        tableOfCrsGraphs().store<Epetra_CrsGraph>(pobj, false));
}

/* store const Epetra_CrsGraph in const table */
CT_Epetra_CrsGraph_ID_t
CEpetra::storeConstCrsGraph( const Epetra_CrsGraph *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(
        tableOfCrsGraphs().store<Epetra_CrsGraph>(pobj, false));
}

/* remove Epetra_CrsGraph from table using CT_Epetra_CrsGraph_ID */
void
CEpetra::removeCrsGraph( CT_Epetra_CrsGraph_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(*id);
    if (tableOfCrsGraphs().isType(aid.table))
        tableOfCrsGraphs().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(aid);
}

/* remove Epetra_CrsGraph from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeCrsGraph( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfCrsGraphs().isType(aid->table))
        tableOfCrsGraphs().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_CrsGraph table */
void
CEpetra::purgeCrsGraph(  )
{
    tableOfCrsGraphs().purge();
}

/* store Epetra_CrsGraph in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasCrsGraph( const Teuchos::RCP< Epetra_CrsGraph > & robj )
{
    return tableOfCrsGraphs().alias(robj);
}

/* store const Epetra_CrsGraph in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstCrsGraph( const Teuchos::RCP< const Epetra_CrsGraph > & robj )
{
    return tableOfCrsGraphs().alias(robj);
}



