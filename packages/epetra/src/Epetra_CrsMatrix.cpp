/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"
#include "Epetra_OffsetIndex.h"
#include "Epetra_BLAS_wrappers.h"


#include <cstdlib>

#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
# include "Teuchos_VerboseObject.hpp"
bool Epetra_CrsMatrixTraceDumpMultiply = false;
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY

#ifdef HAVE_EPETRA_TEUCHOS
// Define this macro to see some timers for some of these functions
# define EPETRA_CRSMATRIX_TEUCHOS_TIMERS
#endif

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
#  include "Teuchos_TimeMonitor.hpp"
#endif

//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& rowMap, const int* NumEntriesPerRow, bool StaticProfile) 
  : Epetra_DistObject(rowMap, "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(CV, rowMap, NumEntriesPerRow, StaticProfile),
    Allocated_(false),
    StaticGraph_(false),
    UseTranspose_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    StorageOptimized_(false),
    Values_(0),
    Values_alloc_lengths_(0),
    All_Values_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    NumMyRows_(rowMap.NumMyPoints()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  Allocate();
}

//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& rowMap, int NumEntriesPerRow, bool StaticProfile) 
  : Epetra_DistObject(rowMap, "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(CV, rowMap, NumEntriesPerRow, StaticProfile),
    Allocated_(false),
    StaticGraph_(false),
    UseTranspose_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    StorageOptimized_(false),
    Values_(0),
    Values_alloc_lengths_(0),
    All_Values_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    NumMyRows_(rowMap.NumMyPoints()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  Allocate();
}
//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& rowMap, 
				   const Epetra_Map& colMap, const int* NumEntriesPerRow, bool StaticProfile) 
  : Epetra_DistObject(rowMap, "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(CV, rowMap, colMap, NumEntriesPerRow, StaticProfile),
    Allocated_(false),
    StaticGraph_(false),
    UseTranspose_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    StorageOptimized_(false),
    Values_(0),
    Values_alloc_lengths_(0),
    All_Values_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    NumMyRows_(rowMap.NumMyPoints()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  Allocate();
}

//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& rowMap, 
				   const Epetra_Map& colMap, int NumEntriesPerRow, bool StaticProfile) 
  : Epetra_DistObject(rowMap, "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(CV, rowMap, colMap,  NumEntriesPerRow, StaticProfile),
    Allocated_(false),
    StaticGraph_(false),
    UseTranspose_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    StorageOptimized_(false),
    Values_(0),
    Values_alloc_lengths_(0),
    All_Values_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    NumMyRows_(rowMap.NumMyPoints()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  Allocate();
}
//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph& graph) 
  : Epetra_DistObject(graph.Map(), "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(graph),
    Allocated_(false),
    StaticGraph_(true),
    UseTranspose_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    StorageOptimized_(false),
    Values_(0),
    Values_alloc_lengths_(0),
    All_Values_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    NumMyRows_(graph.NumMyRows()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  constructedWithFilledGraph_ = graph.Filled();
  InitializeDefaults();
  Allocate();
}

//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(const Epetra_CrsMatrix& Matrix) 
  : Epetra_DistObject(Matrix),
    Epetra_CompObject(Matrix),
    Epetra_BLAS(),
    Graph_(Matrix.Graph()),
    Allocated_(false),
    StaticGraph_(true),
    UseTranspose_(Matrix.UseTranspose_),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    StorageOptimized_(false),
    Values_(0),
    Values_alloc_lengths_(0),
    All_Values_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    NumMyRows_(Matrix.NumMyRows()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(Copy),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  operator=(Matrix);
}

//==============================================================================
Epetra_CrsMatrix& Epetra_CrsMatrix::operator=(const Epetra_CrsMatrix& src)
{
  if (this == &src) {
    return( *this );
  }

  if (!src.Filled()) throw ReportError("Copying an Epetra_CrsMatrix requires source matrix to have Filled()==true", -1);

  Graph_ = src.Graph_; // Copy graph

  DeleteMemory();

  StaticGraph_ = true;
  UseTranspose_ = src.UseTranspose_;
  constructedWithFilledGraph_ = src.constructedWithFilledGraph_;
  matrixFillCompleteCalled_ = src.matrixFillCompleteCalled_;
  Values_ = 0;
  Values_alloc_lengths_ = 0;
  All_Values_ = 0;
  NormInf_ = -1.0;
  NormOne_ = -1.0;
  NormFrob_ = -1.0;
  NumMyRows_ = src.NumMyRows_;
  ImportVector_ = 0;
  ExportVector_ = 0;

  CV_ = Copy;

  StorageOptimized_ = src.StorageOptimized_;
  if (src.StorageOptimized()) { // Special copy for case where storage is optimized

    int numMyNonzeros = Graph().NumMyEntries();
    if (numMyNonzeros>0) All_Values_ = new double[numMyNonzeros];
    double * srcValues = src.All_Values();
    for (int i=0; i<numMyNonzeros; ++i) All_Values_[i] = srcValues[i];
    Allocated_ = true;
#ifdef Epetra_ENABLE_CASK
    if( matrixFillCompleteCalled_  ) {
      cask = cask_handler_copy(src.cask);
    }
#endif
  }
  else { // copy for non-optimized storage
    
    Allocate();
    for (int i=0; i<NumMyRows_; i++) {
      int NumEntries = src.NumMyEntries(i);
      double * const srcValues = src.Values(i);
      double * targValues = Values(i);
      for (int j=0; j< NumEntries; j++) targValues[j] = srcValues[j];
    }
  }

  return( *this );
}

//==============================================================================
void Epetra_CrsMatrix::InitializeDefaults() { // Initialize all attributes that have trivial default values

  UseTranspose_ = false;
  Values_ = 0;
  Values_alloc_lengths_ = 0;
  All_Values_ = 0;
  NormInf_ = -1.0;
  NormOne_ = -1.0;
  NormFrob_ = -1.0;
  ImportVector_ = 0;
  ExportVector_ = 0;

  return;
}

//==============================================================================
int Epetra_CrsMatrix::Allocate() {

  int i, j;

  // Allocate Values array
  Values_ = NumMyRows_ > 0 ? new double*[NumMyRows_] : NULL;
  Values_alloc_lengths_ = NumMyRows_ > 0 ? new int[NumMyRows_] : NULL;
  if (NumMyRows_ > 0) {
    for(j=0; j<NumMyRows_; ++j) Values_alloc_lengths_[j] = 0;
  }

  // Allocate and initialize entries if we are copying data
  if (CV_==Copy) {
    if (Graph().StaticProfile() || Graph().StorageOptimized()) {
      int numMyNonzeros = Graph().NumMyEntries();
      if (numMyNonzeros>0) All_Values_ = new double[numMyNonzeros];
      if(Graph().StorageOptimized()){
        StorageOptimized_ = true;
      }
    }
    double * all_values = All_Values_;
    for (i=0; i<NumMyRows_; i++) {
      int NumAllocatedEntries = Graph().NumAllocatedMyIndices(i);

      if (NumAllocatedEntries > 0) {
        if (Graph().StaticProfile() || Graph().StorageOptimized()) {
          Values_[i] = all_values;
          all_values += NumAllocatedEntries;
        }
        else {
          Values_[i] = new double[NumAllocatedEntries];
          Values_alloc_lengths_[i] = NumAllocatedEntries;
        }
      }
      else 
        Values_[i] = 0;

      for(j=0; j< NumAllocatedEntries; j++) 
        Values_[i][j] = 0.0; // Fill values with zero
    }
  }	 
  else {
    for (i=0; i<NumMyRows_; i++) {
      Values_[i] = 0;
    }
  }
  SetAllocated(true);
#ifdef Epetra_ENABLE_CASK
  cask=NULL;
#endif
  return(0);
}
//==============================================================================
Epetra_CrsMatrix::~Epetra_CrsMatrix()
{
  DeleteMemory();
}

//==============================================================================
void Epetra_CrsMatrix::DeleteMemory()
{
  int i;

  if (CV_==Copy) {
    if (All_Values_!=0)
      delete [] All_Values_;
    else if (Values_!=0)
      for (i=0; i<NumMyRows_; i++) 
        if (Graph().NumAllocatedMyIndices(i) >0) 
          delete [] Values_[i];
  }

  if (ImportVector_!=0) 
    delete ImportVector_;
  ImportVector_=0;

  if (ExportVector_!=0)
    delete ExportVector_;
  ExportVector_=0;

  delete [] Values_;
  Values_ = 0;

  delete [] Values_alloc_lengths_;
  Values_alloc_lengths_ = 0;

  NumMyRows_ = 0;

#ifdef Epetra_ENABLE_CASK
  if( StorageOptimized_  )
  {
    if( cask != NULL ) 
      cask_handler_destroy(cask);

    cask = NULL;
  }
#endif

  Allocated_ = false;
}

//==============================================================================
int Epetra_CrsMatrix::ReplaceRowMap(const Epetra_BlockMap& newmap)
{
  int err = Graph_.ReplaceRowMap(newmap);
  if (err == 0) {
    //update export vector.

    if (ExportVector_ != 0) {
      delete ExportVector_; 
      ExportVector_= 0;
    }

    ExportVector_ = new Epetra_MultiVector(RowMap(),1);
  }
  return(err);
}

//==============================================================================
int Epetra_CrsMatrix::ReplaceColMap(const Epetra_BlockMap& newmap)
{
  int err = Graph_.ReplaceColMap(newmap);
  if (err == 0) {
    //update import vector.

    if (ImportVector_ != 0) {
      delete ImportVector_;
      ImportVector_= 0;
    }

    ImportVector_ = new Epetra_MultiVector(ColMap(),1);
  }
  return(err);
}

//==============================================================================
int Epetra_CrsMatrix::PutScalar(double ScalarConstant) 
{
  if (StorageOptimized()) {
    int length = NumMyNonzeros();
    for (int i=0; i<length; ++i) All_Values_[i] = ScalarConstant;
  }
  else {
    for(int i=0; i<NumMyRows_; i++) {
      int NumEntries = Graph().NumMyIndices(i);
      double * targValues = Values(i);
      for(int j=0; j< NumEntries; j++) 
	targValues[j] = ScalarConstant;
    }
  }
  return(0);
}
//==============================================================================
int Epetra_CrsMatrix::Scale(double ScalarConstant) 
{
  if (StorageOptimized()) {
    int length = NumMyNonzeros();
    for (int i=0; i<length; ++i) All_Values_[i] *= ScalarConstant;
  }
  else {
    for(int i=0; i<NumMyRows_; i++) {
      int NumEntries = Graph().NumMyIndices(i);
      double * targValues = Values(i);
      for(int j=0; j< NumEntries; j++) 
	targValues[j] *= ScalarConstant;
    }
  }
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::InsertGlobalValues(int Row, int NumEntries,
					 const double* values,
					 const int* Indices)
{
  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // Cannot insert global values into local graph
  if(IndicesAreContiguous()) 
    EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and newed
  Graph_.SetIndicesAreGlobal(true);
  Row = Graph_.LRID(Row); // Find local row number for this global row index

  EPETRA_CHK_ERR( InsertValues(Row, NumEntries, values, Indices) );

  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::InsertGlobalValues(int Row, int NumEntries,
					 double* values,
					 int* Indices)
{
  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // Cannot insert global values into local graph
  if(IndicesAreContiguous()) 
    EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and newed
  Graph_.SetIndicesAreGlobal(true);
  Row = Graph_.LRID(Row); // Find local row number for this global row index

  EPETRA_CHK_ERR( InsertValues(Row, NumEntries, values, Indices) );

  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::InsertMyValues(int Row, int NumEntries,
				     const double* values,
				     const int* Indices)
{
  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // Cannot insert global values into filled graph
  if(IndicesAreContiguous() && CV_==Copy) 
    EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and new
  Graph_.SetIndicesAreLocal(true);

  EPETRA_CHK_ERR( InsertValues(Row, NumEntries, values, Indices) );

  return(0);

}

//==========================================================================
int Epetra_CrsMatrix::InsertMyValues(int Row, int NumEntries,
				     double* values,
				     int* Indices)
{
  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // Cannot insert global values into filled graph
  if(IndicesAreContiguous() && CV_==Copy) 
    EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and new
  Graph_.SetIndicesAreLocal(true);

  EPETRA_CHK_ERR( InsertValues(Row, NumEntries, values, Indices) );

  return(0);

}

//==========================================================================
int Epetra_CrsMatrix::InsertValues(int Row, int NumEntries,
				   const double* values,
				   const int* Indices)
{
  if(CV_ == View){
    //cannot allow View mode with const pointers
    EPETRA_CHK_ERR(-4);
  }
  else{
    //to avoid code duplication I am using a cheap tactic of removing the constness of
    //values and Indices. Since this is only called in copy mode the passed in variables
    //will not be modified.
    return(InsertValues(Row, NumEntries, const_cast<double*>(values), const_cast<int*>(Indices)));
  }
  return 0;
}

//==========================================================================
int Epetra_CrsMatrix::InsertValues(int Row, int NumEntries,
				   double* values,
				   int* Indices)
{
  int j;
  int ierr = 0;

  if(Row < 0 || Row >= NumMyRows_) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  if(CV_ == View) {
    //test indices in static graph
    if(StaticGraph()) {
      int testNumEntries;
      int* testIndices;
      int testRow = Row;
      if(IndicesAreGlobal()) 
        testRow = Graph_.LRID( Row );
      EPETRA_CHK_ERR(Graph_.ExtractMyRowView(testRow, testNumEntries, testIndices));

      bool match = true;
      if(NumEntries != testNumEntries) 
        match = false;
      for(int i = 0; i < NumEntries; ++i)
        match = match && (Indices[i]==testIndices[i]);

      if(!match)
        ierr = -3;
    }

    if(Values_[Row] != 0) 
      ierr = 2; // This row has been defined already.  Issue warning.
    Values_[Row] = values;
  }
  else {    
    if(StaticGraph()) 
      EPETRA_CHK_ERR(-2); // If the matrix graph is fully constructed, we cannot insert new values

    int tmpNumEntries = NumEntries;

    if(Graph_.HaveColMap()) { //must insert only valid indices, values
      const double* tmpValues = values;
      values = new double[NumEntries];
      int loc = 0;
      if(IndicesAreLocal()) {
        for(int i = 0; i < NumEntries; ++i)
          if(Graph_.ColMap().MyLID(Indices[i])) 
            values[loc++] = tmpValues[i];
      }
      else {
        for(int i = 0; i < NumEntries; ++i)
          if(Graph_.ColMap().MyGID(Indices[i])) 
            values[loc++] = tmpValues[i];
      }
      if(NumEntries != loc) 
        ierr = 2; //Some columns excluded
      NumEntries = loc;
    } 

    int start = Graph().NumMyIndices(Row);
    int stop = start + NumEntries;
    int NumAllocatedEntries = Values_alloc_lengths_[Row];
    if(stop > NumAllocatedEntries) {
      if (Graph().StaticProfile() && stop > Graph().NumAllocatedMyIndices(Row)) {
        EPETRA_CHK_ERR(-2); // Cannot expand graph storage if graph created using StaticProfile
      }
      if(NumAllocatedEntries == 0) {
        Values_[Row] = new double[NumEntries]; // Row was never allocated, so do it
        Values_alloc_lengths_[Row] = NumEntries;
      }
      else {
        ierr = 1; // Out of room.  Must delete and allocate more space...
        double* tmp_Values = new double[stop];
        for(j = 0; j < start; j++) 
          tmp_Values[j] = Values_[Row][j]; // Copy existing entries
        delete[] Values_[Row]; // Delete old storage
        Values_[Row] = tmp_Values; // Set pointer to new storage
        Values_alloc_lengths_[Row] = stop;
      }
    }

    for(j = start; j < stop; j++) 
      Values_[Row][j] = values[j-start];


    NumEntries = tmpNumEntries;
    if(Graph_.HaveColMap()) 
      delete[] values;
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  if(!StaticGraph()) {
    EPETRA_CHK_ERR(Graph_.InsertIndices(Row, NumEntries, Indices));
  }

  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::InsertOffsetValues(int Row, int NumEntries,
					 double* values,
					 int* Indices)
{
  return ReplaceOffsetValues(Row, NumEntries, values, Indices);
}

//==========================================================================
int Epetra_CrsMatrix::ReplaceGlobalValues(int Row, int NumEntries, const double * srcValues, const int *Indices) {

  int j;
  int ierr = 0;
  int Loc;

  Row = Graph_.LRID(Row); // Normalize row range
    
  if (Row < 0 || Row >= NumMyRows_) {
    EPETRA_CHK_ERR(-1); // Not in Row range
  }
  double * targValues = Values(Row);
  for (j=0; j<NumEntries; j++) {
    int Index = Indices[j];
    if (Graph_.FindGlobalIndexLoc(Row,Index,j,Loc)) 
      targValues[Loc] = srcValues[j];
    else 
      ierr = 2; // Value Excluded
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::ReplaceMyValues(int Row, int NumEntries, const double * srcValues, const int *Indices) {

  if (!IndicesAreLocal()) 
    EPETRA_CHK_ERR(-4); // Indices must be local.

  int j;
  int ierr = 0;
  int Loc;

  if (Row < 0 || Row >= NumMyRows_) {
    EPETRA_CHK_ERR(-1); // Not in Row range
  }

  double* RowValues = Values(Row); 
  for (j=0; j<NumEntries; j++) {
    int Index = Indices[j];
    if (Graph_.FindMyIndexLoc(Row,Index,j,Loc)) 
      RowValues[Loc] = srcValues[j];
    else 
      ierr = 2; // Value Excluded
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::ReplaceOffsetValues(int Row, int NumEntries,
					  const double * srcValues, const int *Offsets)
{
  int j;
  int ierr = 0;

  Row = Graph_.LRID(Row); // Normalize row range
    
  if (Row < 0 || Row >= NumMyRows_) {
    EPETRA_CHK_ERR(-1); // Not in Row range
  }

  double* RowValues = Values(Row); 
  for(j=0; j<NumEntries; j++) {
    if( Offsets[j] != -1 )
      RowValues[Offsets[j]] = srcValues[j];
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::SumIntoGlobalValues(int Row,
					  int NumEntries,
					  const double * srcValues,
					  const int *Indices)
{
  int j;
  int ierr = 0;
  int Loc = 0;

  Row = Graph_.LRID(Row); // Normalize row range
    
  if (Row < 0 || Row >= NumMyRows_) {
    EPETRA_CHK_ERR(-1); // Not in Row range
  }

  if (StaticGraph() && !Graph_.HaveColMap()) {
    EPETRA_CHK_ERR(-1);
  }

  double * RowValues = Values(Row);

  if (!StaticGraph()) {
    for (j=0; j<NumEntries; j++) {
      int Index = Indices[j];
      if (Graph_.FindGlobalIndexLoc(Row,Index,j,Loc))
        RowValues[Loc] += srcValues[j];
      else
        ierr = 2; // Value Excluded
    }
  }
  else {
    const Epetra_BlockMap& colmap = Graph_.ColMap();
    int NumColIndices = Graph_.NumMyIndices(Row);
    const int* ColIndices = Graph_.Indices(Row);

    if (Graph_.Sorted()) {
      int insertPoint;
      for (j=0; j<NumEntries; j++) {
        int Index = colmap.LID(Indices[j]);

        // Check whether the next added element is the subsequent element in
        // the graph indices, then we can skip the binary search
        if (Loc < NumColIndices && Index == ColIndices[Loc])
          RowValues[Loc] += srcValues[j];
        else {
          Loc = Epetra_Util_binary_search(Index, ColIndices, NumColIndices, insertPoint);
          if (Loc > -1)
            RowValues[Loc] += srcValues[j];
          else 
            ierr = 2; // Value Excluded
        }
        ++Loc;
      }
    }
    else
      for (j=0; j<NumEntries; j++) {
        int Index = colmap.LID(Indices[j]);
        if (Graph_.FindMyIndexLoc(NumColIndices,ColIndices,Index,j,Loc)) 
          RowValues[Loc] += srcValues[j];
        else 
          ierr = 2; // Value Excluded
      }
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);

  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::SumIntoMyValues(int Row, int NumEntries, const double * srcValues, const int *Indices) {

  if (!IndicesAreLocal()) 
    EPETRA_CHK_ERR(-4); // Indices must be local.

  int j;
  int ierr = 0;
  int Loc = 0;
  int insertPoint;

  if (Row < 0 || Row >= NumMyRows_) {
    EPETRA_CHK_ERR(-1); // Not in Row range
  }

  double* RowValues = Values(Row);
  int NumColIndices = Graph_.NumMyIndices(Row);
  const int* ColIndices = Graph_.Indices(Row);
  if (Graph_.Sorted()) {
    for (j=0; j<NumEntries; j++) {
      int Index = Indices[j];

      // Check whether the next added element is the subsequent element in
      // the graph indices.
      if (Loc < NumColIndices && Index == ColIndices[Loc])
        RowValues[Loc] += srcValues[j];
      else {
        Loc = Epetra_Util_binary_search(Index, ColIndices, NumColIndices, insertPoint);
        if (Loc > -1)
          RowValues[Loc] += srcValues[j];
        else 
          ierr = 2; // Value Excluded
      }
      ++Loc;
    }
  }
  else {
    for (j=0; j<NumEntries; j++) {
      int Index = Indices[j];
      if (Graph_.FindMyIndexLoc(Row,Index,j,Loc)) 
        RowValues[Loc] += srcValues[j];
      else 
        ierr = 2; // Value Excluded
    }
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);

  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::SumIntoOffsetValues(int Row, int NumEntries, const double * srcValues, const int *Offsets) {

  int j;
  int ierr = 0;

  Row = Graph_.LRID(Row); // Normalize row range
    
  if (Row < 0 || Row >= NumMyRows_) {
    EPETRA_CHK_ERR(-1); // Not in Row range
  }

  double* RowValues = Values(Row);
  for (j=0; j<NumEntries; j++) {
    if( Offsets[j] != -1 )
      RowValues[Offsets[j]] += srcValues[j];
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);

  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::FillComplete(bool OptimizeDataStorage) {
  squareFillCompleteCalled_ = true;
  EPETRA_CHK_ERR(FillComplete(RowMap(), RowMap(), OptimizeDataStorage));
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::FillComplete(const Epetra_Map& domain_map,
				   const Epetra_Map& range_map, bool OptimizeDataStorage)
{
  int returnValue = 0;

  if (Graph_.Filled()) {
    if (!constructedWithFilledGraph_ && !matrixFillCompleteCalled_) {
      returnValue = 2;
    }
  }

  if (!StaticGraph()) {
    if (Graph_.MakeIndicesLocal(domain_map, range_map) < 0) {
      return(-1);
    }
  }
  SortEntries();  // Sort column entries from smallest to largest
  MergeRedundantEntries(); // Get rid of any redundant index values
  if (!StaticGraph()) {
    if (Graph_.FillComplete(domain_map, range_map) < 0) {
      return(-2);
    }
  }

  matrixFillCompleteCalled_ = true;

  if (squareFillCompleteCalled_) {
    if (DomainMap().NumGlobalElements() != RangeMap().NumGlobalElements()) {
      returnValue = 3;
    }
    squareFillCompleteCalled_ = false;
    EPETRA_CHK_ERR(returnValue);
  }

  if (OptimizeDataStorage) { EPETRA_CHK_ERR(OptimizeStorage()); }

  return(returnValue);
}

//==========================================================================
int Epetra_CrsMatrix::TransformToLocal() {
  EPETRA_CHK_ERR(FillComplete());
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::TransformToLocal(const Epetra_Map* domainMap, const Epetra_Map* rangeMap) {
  EPETRA_CHK_ERR(FillComplete(*domainMap, *rangeMap));
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::SortEntries() {

  if(!IndicesAreLocal()) 
    EPETRA_CHK_ERR(-1);
  if(Sorted())
    return(0);

  // For each row, sort column entries from smallest to largest.
  // Use shell sort. Stable sort so it is fast if indices are already sorted.

  
  for(int i = 0; i < NumMyRows_; i++){

    double* locValues = Values(i);
    int NumEntries = Graph().NumMyIndices(i);
    int* locIndices = Graph().Indices(i);
		
    int n = NumEntries;
    int m = n/2;
    
    while(m > 0) {
      int max = n - m;
      for(int j = 0; j < max; j++) {
	for(int k = j; k >= 0; k-=m) {
	  if(locIndices[k+m] >= locIndices[k])
	    break;
	  double dtemp = locValues[k+m];
	  locValues[k+m] = locValues[k];
	  locValues[k] = dtemp;
	  int itemp = locIndices[k+m];
	  locIndices[k+m] = locIndices[k];
	  locIndices[k] = itemp;
	}
      }
      m = m/2;
    }
  }
  Graph_.SetSorted(true); // This also sorted the graph
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::MergeRedundantEntries() {

  int i;

  if(NoRedundancies()) 
    return(0);
  if(!Sorted()) 
    EPETRA_CHK_ERR(-1);  // Must have sorted entries

  // For each row, remove column indices that are repeated.
  // Also, determine if matrix is upper or lower triangular or has no diagonal (Done in graph)
  // Note:  This function assumes that SortEntries was already called.

  for(i = 0; i<NumMyRows_; i++) {
    int NumEntries = Graph().NumMyIndices(i);
    if(NumEntries > 1) {
      double* const locValues = Values(i);
      int* const locIndices = Graph().Indices(i);		
      int curEntry =0;
      double curValue = locValues[0];
      for(int k = 1; k < NumEntries; k++) {
	if(locIndices[k] == locIndices[k-1]) 
	  curValue += locValues[k];
	else {
	  locValues[curEntry++] = curValue;
	  curValue = locValues[k];
	}
      }
      locValues[curEntry] = curValue;
      
    }
  }
  
  EPETRA_CHK_ERR(Graph_.RemoveRedundantIndices()); // Remove redundant indices and then return
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::OptimizeStorage() {


  if (StorageOptimized()) 
    return(0); // Have we been here before?
  if (!Filled()) EPETRA_CHK_ERR(-1); // Cannot optimize storage before calling FillComplete()


  int ierr = Graph_.OptimizeStorage();
  if (ierr!=0) EPETRA_CHK_ERR(ierr);  // In order for OptimizeStorage to make sense for the matrix, it must work on the graph.

  bool Contiguous = true; // Assume contiguous is true
  for (int i=1; i<NumMyRows_; i++){
    int NumEntries = Graph().NumMyIndices(i-1);
		
    // check if end of beginning of current row starts immediately after end of previous row.
    if (Values_[i]!=Values_[i-1]+NumEntries) {
      Contiguous = false;
      break;
    }
  }

  // NOTE:  At the end of the above loop set, there is a possibility that NumEntries and NumAllocatedEntries
  //        for the last row could be different, but I don't think it matters.


  if ((CV_==View) && !Contiguous) 
    EPETRA_CHK_ERR(-1);  // This is user data, it's not contiguous and we can't make it so.

  
  if(!Contiguous) { // Must pack indices if not already contiguous.  Pack values into All_values_
    
    if (All_Values_==0) {
      // Compute Number of Nonzero entries (Done in FillComplete, but we may not have been there yet.)
      int numMyNonzeros = Graph_.NumMyNonzeros();
    
      // Allocate one big array for all values
      All_Values_ = new double[numMyNonzeros];
      if(All_Values_ == 0) throw ReportError("Error with All_Values_ allocation.", -99);

      const int *const  IndexOffset = Graph().IndexOffset();
      const int numMyRows = NumMyRows_;
      double ** Values_s = Values_;
      double * All_Values_s = All_Values_;
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(Values_s,All_Values_s)
#endif
      for (int i=0; i<numMyRows; i++) {
        int NumEntries = Graph().NumMyIndices(i);
        int curOffset = IndexOffset[i];
        double * values = Values_s[i];
        double * newValues = All_Values_s+curOffset;
	for (int j=0; j<NumEntries; j++) newValues[j] = values[j];
      }
    }
    else { // Static Profile, so just pack into existing storage (can't be threaded)
      double * tmp = All_Values_;
      for (int i=0; i<NumMyRows_; i++) {
        int NumEntries = Graph().NumMyIndices(i);
        double * values = Values_[i];
        if (tmp!=values) // Copy values if not pointing to same location
          for (int j=0; j<NumEntries; j++) tmp[j] = values[j];
        tmp += NumEntries;
      }
    }

    // Free Values_ arrays
    for (int i=0;i<NumMyRows_; ++i) {
      if (Values_alloc_lengths_[i] != 0) delete [] Values_[i];
    }
    
    delete [] Values_alloc_lengths_; Values_alloc_lengths_ = 0;
  } // End of !Contiguous section
  else {
    //if already contiguous, we'll simply set All_Values_ to be
    //a copy of Values_[0].
    All_Values_ = NumMyRows_ > 0 ? Values_[0] : NULL;
  }
  
  // Delete unneeded storage
  delete [] Values_; Values_=0;

#ifdef Epetra_ENABLE_CASK 
  if (cask == NULL  && Graph().StorageOptimized() )  {
     int * Indices = Graph().All_Indices();
     int * IndexOffset = Graph().IndexOffset();
     int NumMyCols_ = NumMyCols();
     cask_handler_initialize(&cask);
     cask_csr_analysis(NumMyRows_, NumMyCols_, IndexOffset, Indices,
                       NumGlobalNonzeros(),cask);
  }
#endif


  StorageOptimized_ = true;

  
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::ExtractGlobalRowCopy(int Row, int Length, int & NumEntries, double * values,
					   int * Indices) const 
{

  int ierr = Graph_.ExtractGlobalRowCopy(Row, Length, NumEntries, Indices);
  if (ierr) 
    EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractGlobalRowCopy(Row, Length, NumEntries, values));
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::ExtractMyRowCopy(int Row, int Length, int & NumEntries, double * values,
				       int * Indices) const 
{

  int ierr = Graph_.ExtractMyRowCopy(Row, Length, NumEntries, Indices);
  if (ierr) 
    EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractMyRowCopy(Row, Length, NumEntries, values));
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::NumMyRowEntries(int Row, int & NumEntries) const 
{

  if (!MyLRID(Row)) 
    EPETRA_CHK_ERR(-1); // Not in the range of local rows
  NumEntries = NumMyEntries(Row);
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::ExtractGlobalRowCopy(int Row, int Length, int & NumEntries, double * values) const 
{

  int Row0 = Graph_.RowMap().LID(Row); // Normalize row range

  EPETRA_CHK_ERR(ExtractMyRowCopy(Row0, Length, NumEntries, values));
  return(0);
}


//==========================================================================
int Epetra_CrsMatrix::ExtractMyRowCopy(int Row, int Length, int & NumEntries, double * targValues) const 
{
  int j;

  if (Row < 0 || Row >= NumMyRows_) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  NumEntries = Graph().NumMyIndices(Row);
  if (Length < NumEntries) 
    EPETRA_CHK_ERR(-2); // Not enough space for copy. Needed size is passed back in NumEntries

  double * srcValues = Values(Row);

  for(j=0; j<NumEntries; j++)
    targValues[j] = srcValues[j];
  
  return(0);
}


//==============================================================================
int Epetra_CrsMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const {
	
  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Can't get diagonal unless matrix is filled (and in local index space)
  if(!RowMap().SameAs(Diagonal.Map())) 
    EPETRA_CHK_ERR(-2); // Maps must be the same

  for(int i = 0; i < NumMyRows_; i++) {
    int ii = GRID(i);
    int NumEntries = Graph().NumMyIndices(i);
    int* Indices = Graph().Indices(i);
    double * srcValues = Values(i);
    
    Diagonal[i] = 0.0;
    for(int j = 0; j < NumEntries; j++) {
      if(ii == GCID(Indices[j])) {
	Diagonal[i] = srcValues[j];
	break;
      }
    }
  }
  return(0);
}
//==============================================================================
int Epetra_CrsMatrix::ReplaceDiagonalValues(const Epetra_Vector & Diagonal) {
	
  if(!Filled())
    EPETRA_CHK_ERR(-1); // Can't replace diagonal unless matrix is filled (and in local index space)
  if(!RowMap().SameAs(Diagonal.Map())) 
    EPETRA_CHK_ERR(-2); // Maps must be the same

  int ierr = 0;
  for(int i = 0; i < NumMyRows_; i++) {
    int ii = GRID(i);
    int NumEntries = Graph().NumMyIndices(i);
    int* Indices = Graph().Indices(i);
    double * targValues = Values(i);
    bool DiagMissing = true;
    for(int j = 0; j < NumEntries; j++) {
      if(ii == GCID(Indices[j])) {
	targValues[j] = Diagonal[i];
	DiagMissing = false;
	break;
      }
    }
    if(DiagMissing) 
      ierr = 1; // flag a warning error
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);

  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::ExtractGlobalRowView(int Row, int & NumEntries, double *& values, int *& Indices) const 
{

  int ierr = Graph_.ExtractGlobalRowView(Row, NumEntries, Indices);
  if (ierr) 
    EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractGlobalRowView(Row, NumEntries, values));
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::ExtractMyRowView(int Row, int & NumEntries, double *& values, int *& Indices) const 
{
  int ierr = Graph_.ExtractMyRowView(Row, NumEntries, Indices);
  if (ierr) 
    EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractMyRowView(Row, NumEntries, values));
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::ExtractGlobalRowView(int Row, int & NumEntries, double *& values) const 
{

  int Row0 = Graph_.RowMap().LID(Row); // Normalize row range

  EPETRA_CHK_ERR(ExtractMyRowView(Row0, NumEntries, values));
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::ExtractMyRowView(int Row, int & NumEntries, double *& targValues) const 
{

  if (Row < 0 || Row >= NumMyRows_) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  NumEntries = Graph().NumMyIndices(Row);

  targValues = Values(Row);
  
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal,
			    const Epetra_Vector& x, Epetra_Vector& y) const
{

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Solve(Upper,Trans,UnitDiag,x,y)");
#endif

  //
  // This function finds y such that Ly = x or Uy = x or the transpose cases.
  //

  // First short-circuit to the pre-5.0 version if no storage optimization was performed
  if (!StorageOptimized() && !Graph().StorageOptimized()) {
    EPETRA_CHK_ERR(Solve1(Upper, Trans, UnitDiagonal, x, y));
    return(0);
  }

  if (!Filled()) {
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  }

  if ((Upper) && (!UpperTriangular())) 
    EPETRA_CHK_ERR(-2);
  if ((!Upper) && (!LowerTriangular())) 
    EPETRA_CHK_ERR(-3);
  if ((!UnitDiagonal) && (NoDiagonal())) 
    EPETRA_CHK_ERR(-4); // If matrix has no diagonal, we must use UnitDiagonal
  if ((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_)) 
    EPETRA_CHK_ERR(-5); // Need each row to have a diagonal

  double *xp = (double*)x.Values();
  double *yp = (double*)y.Values();

      
  GeneralSV(Upper, Trans, UnitDiagonal, xp, yp);

  UpdateFlops(2*NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Solve(Upper,Trans,UnitDiag,X,Y)");
#endif

  //
  // This function find Y such that LY = X or UY = X or the transpose cases.
  //

  // First short-circuit to the pre-5.0 version if no storage optimization was performed
  if (!StorageOptimized() && !Graph().StorageOptimized()) {
    EPETRA_CHK_ERR(Solve1(Upper, Trans, UnitDiagonal, X, Y));
    return(0);
  }
  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  if((Upper) && (!UpperTriangular()))
    EPETRA_CHK_ERR(-2);
  if((!Upper) && (!LowerTriangular()))
    EPETRA_CHK_ERR(-3);
  if((!UnitDiagonal) && (NoDiagonal()))
    EPETRA_CHK_ERR(-4); // If matrix has no diagonal, we must use UnitDiagonal
  if((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_))
    EPETRA_CHK_ERR(-5); // Need each row to have a diagonal

  double** Xp = (double**) X.Pointers();
  double** Yp = (double**) Y.Pointers();
  int LDX = X.ConstantStride() ? X.Stride() : 0;
  int LDY = Y.ConstantStride() ? Y.Stride() : 0;
  int NumVectors = X.NumVectors();


    // Do actual computation
  if (NumVectors==1)
    GeneralSV(Upper, Trans, UnitDiagonal, *Xp, *Yp);
  else
    GeneralSM(Upper, Trans, UnitDiagonal, Xp, LDX, Yp, LDY, NumVectors);

  UpdateFlops(2 * NumVectors * NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::InvRowSums(Epetra_Vector& x) const {
  //
  // Put inverse of the sum of absolute values of the ith row of A in x[i].
  //

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  int ierr = 0;
  int i, j;
  x.PutScalar(0.0); // Make sure we sum into a vector of zeros.
  double * xp = (double*)x.Values();
  if (Graph().RangeMap().SameAs(x.Map()) && Exporter() != 0) {
    Epetra_Vector x_tmp(RowMap());
    x_tmp.PutScalar(0.0);
    double * x_tmp_p = (double*)x_tmp.Values();
    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = NumMyEntries(i);
      double * RowValues  = Values(i);
      for (j=0; j < NumEntries; j++)  x_tmp_p[i] += std::abs(RowValues[j]);
    }
    EPETRA_CHK_ERR(x.Export(x_tmp, *Exporter(), Add)); //Export partial row sums to x.
    int myLength = x.MyLength();
    for (i=0; i<myLength; i++) { 
      if (xp[i]<Epetra_MinDouble) {
        if (xp[i]==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
        else if (ierr!=1) ierr = 2;
        xp[i] = Epetra_MaxDouble;
      }
      else
        xp[i] = 1.0/xp[i];
    }
  }
  else if (Graph().RowMap().SameAs(x.Map())) {
    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = NumMyEntries(i);
      double * RowValues  = Values(i);
      double scale = 0.0;
      for (j=0; j < NumEntries; j++) scale += std::abs(RowValues[j]);
      if (scale<Epetra_MinDouble) {
        if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
        else if (ierr!=1) ierr = 2;
        xp[i] = Epetra_MaxDouble;
      }
      else
        xp[i] = 1.0/scale;
    }
  }
  else { // x.Map different than both Graph().RowMap() and Graph().RangeMap()
    EPETRA_CHK_ERR(-2); // The map of x must be the RowMap or RangeMap of A.
  }
  UpdateFlops(NumGlobalNonzeros());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::InvRowMaxs(Epetra_Vector& x) const {
  //
  // Put inverse of the max of absolute values of the ith row of A in x[i].
  //

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  int ierr = 0;
  int i, j;
  bool needExport = false;
  double * xp = (double*)x.Values();
  Epetra_Vector* x_tmp = 0;
  if (Graph().RangeMap().SameAs(x.Map())) {
    if (Exporter() != 0) {
      needExport = true; //Having this information later avoids a .SameAs
      x_tmp = new Epetra_Vector(RowMap()); // Create import vector if needed
      xp = (double*)x_tmp->Values();
    }
  }
  else if (!Graph().RowMap().SameAs(x.Map())) {
    EPETRA_CHK_ERR(-2); // The map of x must be the RowMap or RangeMap of A.
  }
  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = NumMyEntries(i);
    double * RowValues  = Values(i);
    double scale = 0.0;
    for (j=0; j < NumEntries; j++) scale = EPETRA_MAX(std::abs(RowValues[j]),scale);
    if (scale<Epetra_MinDouble) {
      if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowmax found (supercedes ierr = 2)
      else if (ierr!=1) ierr = 2;
      xp[i] = Epetra_MaxDouble;
    }
    else
      xp[i] = 1.0/scale;
  }
  if(needExport) {
    x.PutScalar(0.0);
    EPETRA_CHK_ERR(x.Export(*x_tmp, *Exporter(), Insert)); // Fill x with values from temp vector
    delete x_tmp;
  }
  UpdateFlops(NumGlobalNonzeros());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::InvColSums(Epetra_Vector& x) const {
  //
  // Put inverse of the sum of absolute values of the jth column of A in x[j].
  //

  if(!Filled())  EPETRA_CHK_ERR(-1); // Matrix must be filled.
  int ierr = 0;
  int i, j;
  int MapNumMyElements = x.Map().NumMyElements();
  x.PutScalar(0.0); // Make sure we sum into a vector of zeros.
  double* xp = (double*)x.Values();
  if(Graph().DomainMap().SameAs(x.Map()) && Importer() != 0) {
    Epetra_Vector x_tmp(ColMap());
    x_tmp.PutScalar(0.0);
    double * x_tmp_p = (double*)x_tmp.Values();
    for(i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    ColIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for(j = 0; j < NumEntries; j++) 
        x_tmp_p[ColIndices[j]] += std::abs(RowValues[j]);
    }
    EPETRA_CHK_ERR(x.Export(x_tmp, *Importer(), Add)); // Fill x with partial column sums
  }
  else if(Graph().ColMap().SameAs(x.Map())) {
    for(i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    ColIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for(j = 0; j < NumEntries; j++) 
        xp[ColIndices[j]] += std::abs(RowValues[j]);
    }
  }
  else { //x.Map different than both Graph().ColMap() and Graph().DomainMap()
    EPETRA_CHK_ERR(-2); // x must have the same distribution as the domain of A
  }

  // Invert values, don't allow them to get too large
  for(i = 0; i < MapNumMyElements; i++) {
    double scale = xp[i];
    if(scale < Epetra_MinDouble) {
      if(scale == 0.0) 
	ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if(ierr != 1) 
	ierr = 2;
      xp[i] = Epetra_MaxDouble;
    }
    else
      xp[i] = 1.0 / scale;
  }

  UpdateFlops(NumGlobalNonzeros());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::InvColMaxs(Epetra_Vector& x) const {
  //
  // Put inverse of the max of absolute values of the jth column of A in x[j].
  //

  if(!Filled())  EPETRA_CHK_ERR(-1); // Matrix must be filled.
  int ierr = 0;
  int i, j;
  int MapNumMyElements = x.Map().NumMyElements();
  x.PutScalar(0.0); // Make sure we sum into a vector of zeros.
  double* xp = (double*)x.Values();
  if(Graph().DomainMap().SameAs(x.Map()) && Importer() != 0) {
    Epetra_Vector x_tmp(ColMap());
    x_tmp.PutScalar(0.0);
    double * x_tmp_p = (double*)x_tmp.Values();
    for(i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    ColIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for(j = 0; j < NumEntries; j++) 
        x_tmp_p[ColIndices[j]] = EPETRA_MAX(std::abs(RowValues[j]),x_tmp_p[ColIndices[j]]);
    }
    EPETRA_CHK_ERR(x.Export(x_tmp, *Importer(), AbsMax)); // Fill x with partial column sums
  }
  else if(Graph().ColMap().SameAs(x.Map())) {
    for(i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    ColIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for(j = 0; j < NumEntries; j++) 
        xp[ColIndices[j]] = EPETRA_MAX(std::abs(RowValues[j]),xp[ColIndices[j]]);
    }
  }
  else { //x.Map different than both Graph().ColMap() and Graph().DomainMap()
    EPETRA_CHK_ERR(-2); // x must have the same distribution as the domain of A
  }

  // Invert values, don't allow them to get too large
  for(i = 0; i < MapNumMyElements; i++) {
    double scale = xp[i];
    if(scale < Epetra_MinDouble) {
      if(scale == 0.0) 
	ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if(ierr != 1) 
	ierr = 2;
      xp[i] = Epetra_MaxDouble;
    }
    else
      xp[i] = 1.0 / scale;
  }

  UpdateFlops(NumGlobalNonzeros());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::LeftScale(const Epetra_Vector& x) {
  //
  // This function scales the ith row of A by x[i].
  //

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  double* xp = 0;
  if(Graph().RangeMap().SameAs(x.Map()))  
    // If we have a non-trivial exporter, we must import elements that are 
    // permuted or are on other processors.  (We will use the exporter to
    // perform the import.)
    if(Exporter() != 0) {
      UpdateExportVector(1);
      EPETRA_CHK_ERR(ExportVector_->Import(x,*Exporter(), Insert));
      xp = (double*) ExportVector_->Values();
    }
    else
      xp = (double*)x.Values();
  else if (Graph().RowMap().SameAs(x.Map()))
    xp = (double*)x.Values();
  else {
    EPETRA_CHK_ERR(-2); // The Map of x must be the RowMap or RangeMap of A.
  }
  int i, j;

  for(i = 0; i < NumMyRows_; i++) {
    int      NumEntries = NumMyEntries(i);
    double* RowValues  = Values(i);
    double scale = xp[i];
    for(j = 0; j < NumEntries; j++)  
      RowValues[j] *= scale;
  }
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  UpdateFlops(NumGlobalNonzeros());

  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::RightScale(const Epetra_Vector& x) {
  //
  // This function scales the jth column of A by x[j].
  //

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  double* xp = 0;
  if(Graph().DomainMap().SameAs(x.Map())) 
    // If we have a non-trivial exporter, we must import elements that are 
    // permuted or are on other processors.
    if(Importer() != 0) {
      UpdateImportVector(1);
      EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer(), Insert));
      xp = (double*) ImportVector_->Values();
    }
    else
      xp = (double*)x.Values();
  else if(Graph().ColMap().SameAs(x.Map()))
    xp = (double*)x.Values(); 
  else
    EPETRA_CHK_ERR(-2); // The Map of x must be the RowMap or RangeMap of A.
  int i, j;

  for(i = 0; i < NumMyRows_; i++) {
    int     NumEntries = NumMyEntries(i);
    int*    ColIndices = Graph().Indices(i);
    double* RowValues  = Values(i);
    for(j = 0; j < NumEntries; j++)  
      RowValues[j] *=  xp[ColIndices[j]];
  }
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  UpdateFlops(NumGlobalNonzeros());
  return(0);
}

//=============================================================================
double Epetra_CrsMatrix::NormInf() const {

#if 0
  //
  //  Commenting this section out disables caching, ie. 
  //  causes the norm to be computed each time NormInf is called.
  //  See bug #1151 for a full discussion.  
  //
  double MinNorm ; 
  Comm().MinAll( &NormInf_, &MinNorm, 1 ) ; 

  if( MinNorm >= 0.0) 
    return(NormInf_);
#endif

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  Epetra_Vector x(RangeMap()); // Need temp vector for row sums
  double* xp = (double*)x.Values();
  Epetra_MultiVector* x_tmp = 0;

  // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
  if(Exporter() != 0) {
    x_tmp = new Epetra_Vector(RowMap()); // Create temporary import vector if needed
    xp = (double*)x_tmp->Values();
  }
  int i, j;

  for(i = 0; i < NumMyRows_; i++) {
    xp[i] = 0.0;
    int     NumEntries = NumMyEntries(i);
    double* RowValues  = Values(i);
    for(j = 0; j < NumEntries; j++) 
      xp[i] += std::abs(RowValues[j]);
  }
  if(Exporter() != 0) {
    x.PutScalar(0.0);
    EPETRA_CHK_ERR(x.Export(*x_tmp, *Exporter(), Add)); // Fill x with Values from temp vector
  }
  x.MaxValue(&NormInf_); // Find max
  if(x_tmp != 0) 
    delete x_tmp;
  UpdateFlops(NumGlobalNonzeros());
  return(NormInf_);
}
//=============================================================================
double Epetra_CrsMatrix::NormOne() const {

#if 0
  //
  //  Commenting this section out disables caching, ie. 
  //  causes the norm to be computed each time NormOne is called.  
  //  See bug #1151 for a full discussion.  
  //
  double MinNorm ; 
  Comm().MinAll( &NormOne_, &MinNorm, 1 ) ; 

  if( MinNorm >= 0.0) 
    return(NormOne_);
#endif

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  Epetra_Vector x(DomainMap()); // Need temp vector for column sums
  
  double* xp = (double*)x.Values();
  Epetra_MultiVector* x_tmp = 0;
  int NumCols = NumMyCols();
  

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if(Importer() != 0) {
    x_tmp = new Epetra_Vector(ColMap()); // Create temporary import vector if needed
    xp = (double*)x_tmp->Values();
  }
  int i, j;

  for(i = 0; i < NumCols; i++) 
    xp[i] = 0.0;

  for(i = 0; i < NumMyRows_; i++) {
    int     NumEntries = NumMyEntries(i);
    int*    ColIndices = Graph().Indices(i);
    double* RowValues  = Values(i);
    for(j = 0; j < NumEntries; j++) 
      xp[ColIndices[j]] += std::abs(RowValues[j]);
  }
  if(Importer() != 0) {
    x.PutScalar(0.0);
    EPETRA_CHK_ERR(x.Export(*x_tmp, *Importer(), Add)); // Fill x with Values from temp vector
  }
  x.MaxValue(&NormOne_); // Find max
  if(x_tmp != 0) 
    delete x_tmp;
  UpdateFlops(NumGlobalNonzeros());
  return(NormOne_);
}
//=============================================================================
double Epetra_CrsMatrix::NormFrobenius() const {

#if 0
  //
  //  Commenting this section out disables caching, ie. 
  //  causes the norm to be computed each time NormFrobenius is called.  
  //  See bug #1151 for a full discussion.  
  //
  double MinNorm ; 
  Comm().MinAll( &NormFrob_, &MinNorm, 1 ) ; 

  if( MinNorm >= 0.0) 
    return(NormFrob_);
#endif

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  double local_sum = 0.0;

  for(int i = 0; i < NumMyRows_; i++) {
    int     NumEntries = NumMyEntries(i);
    double* RowValues  = Values(i);
    for(int j = 0; j < NumEntries; j++) {
      local_sum += RowValues[j]*RowValues[j];
    }
  }

  double global_sum = 0.0;
  Comm().SumAll(&local_sum, &global_sum, 1);

  NormFrob_ = std::sqrt(global_sum);

  UpdateFlops(NumGlobalNonzeros());

  return(NormFrob_);
}
//=========================================================================
int Epetra_CrsMatrix::CheckSizes(const Epetra_SrcDistObject & Source) {
  try {
    const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Source);
    if (!A.Graph().GlobalConstantsComputed()) EPETRA_CHK_ERR(-1); // Must have global constants to proceed
  }
  catch (...) {
    return(0); // No error at this point, object could be a RowMatrix
  }
  return(0);
}
//=========================================================================
int Epetra_CrsMatrix::CopyAndPermute(const Epetra_SrcDistObject & Source,
				     int NumSameIDs, 
				     int NumPermuteIDs,
                                     int * PermuteToLIDs,
				     int *PermuteFromLIDs,
                                     const Epetra_OffsetIndex * Indexor ) {
 
  try {
    const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Source);
    EPETRA_CHK_ERR(CopyAndPermuteCrsMatrix(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs,
					   PermuteFromLIDs,Indexor));
  }
  catch (...) {
    try {
      const Epetra_RowMatrix & A = dynamic_cast<const Epetra_RowMatrix &>(Source);
      EPETRA_CHK_ERR(CopyAndPermuteRowMatrix(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs,
					     PermuteFromLIDs,Indexor));
    }
    catch (...) {
      EPETRA_CHK_ERR(-1); // Incompatible SrcDistObject
    }
  }
  
  return(0);
}

//=========================================================================
int Epetra_CrsMatrix::CopyAndPermuteCrsMatrix(const Epetra_CrsMatrix & A,
                                              int NumSameIDs, 
					      int NumPermuteIDs,
                                              int * PermuteToLIDs,
					      int *PermuteFromLIDs,
                                              const Epetra_OffsetIndex * Indexor) {
  
  int i, ierr;
  
  int Row, NumEntries;
  int maxNumEntries = A.MaxNumEntries();
  int * Indices = 0;
  double * values = 0;

  if (maxNumEntries>0 && A.IndicesAreLocal() ) { //Need Temp Space
    Indices = new int[maxNumEntries];
    values = new double[maxNumEntries];
  }
  
  // Do copy first
  if (NumSameIDs>0) {
    if (A.IndicesAreLocal()) {
      if (StaticGraph() || IndicesAreLocal()) {
        if(Indexor) {
          for (i=0; i<NumSameIDs; i++) {
	    Row = GRID(i);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(Row, maxNumEntries, NumEntries, values, Indices)); // Set pointers
            ierr = ReplaceOffsetValues(Row, NumEntries, values, Indexor->SameOffsets()[i]);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumSameIDs; i++) {
	    Row = GRID(i);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(Row, maxNumEntries, NumEntries, values, Indices)); // Set pointers
            ierr = ReplaceGlobalValues(Row, NumEntries, values, Indices);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
      }
      else {
        if(Indexor) {
          for (i=0; i<NumSameIDs; i++) {
	    Row = GRID(i);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(Row, maxNumEntries, NumEntries, values, Indices)); // Set pointers
            ierr = InsertOffsetValues(Row, NumEntries, values, Indexor->SameOffsets()[i]);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumSameIDs; i++) {
	    Row = GRID(i);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(Row, maxNumEntries, NumEntries, values, Indices)); // Set pointers
            ierr = InsertGlobalValues(Row, NumEntries, values, Indices);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
      } 
    }
    else { // A.IndicesAreGlobal()
      if (StaticGraph() || IndicesAreLocal()) {
        if(Indexor) {
          for (i=0; i<NumSameIDs; i++) {
	    Row = GRID(i);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowView(Row, NumEntries, values, Indices)); // Set pointers
            ierr = ReplaceOffsetValues(Row, NumEntries, values, Indexor->SameOffsets()[i]);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumSameIDs; i++) {
	    Row = GRID(i);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowView(Row, NumEntries, values, Indices)); // Set pointers
            ierr = ReplaceGlobalValues(Row, NumEntries, values, Indices);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
      }
      else {
        if(Indexor) {
          for (i=0; i<NumSameIDs; i++) {
	    Row = GRID(i);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowView(Row, NumEntries, values, Indices)); // Set pointers
            ierr = InsertOffsetValues(Row, NumEntries, values, Indexor->SameOffsets()[i]);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumSameIDs; i++) {
	    Row = GRID(i);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowView(Row, NumEntries, values, Indices)); // Set pointers
            ierr = InsertGlobalValues(Row, NumEntries, values, Indices);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
      } 
    }
  }
  
  // Do local permutation next
  int FromRow, ToRow;
  if (NumPermuteIDs>0) {
    if (A.IndicesAreLocal()) {
      if (StaticGraph() || IndicesAreLocal()) {
        if(Indexor) {
          for (i=0; i<NumPermuteIDs; i++) {
	    FromRow = A.GRID(PermuteFromLIDs[i]);
	    ToRow = GRID(PermuteToLIDs[i]);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, maxNumEntries, NumEntries, values, Indices)); // Set pointers
	    ierr = ReplaceOffsetValues(ToRow, NumEntries, values, Indexor->PermuteOffsets()[i]);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumPermuteIDs; i++) {
	    FromRow = A.GRID(PermuteFromLIDs[i]);
	    ToRow = GRID(PermuteToLIDs[i]);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, maxNumEntries, NumEntries, values, Indices)); // Set pointers
	    ierr = ReplaceGlobalValues(ToRow, NumEntries, values, Indices);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
      }
      else {
        if(Indexor) {
          for (i=0; i<NumPermuteIDs; i++) {
	    FromRow = A.GRID(PermuteFromLIDs[i]);
	    ToRow = GRID(PermuteToLIDs[i]);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, maxNumEntries, NumEntries, values, Indices)); // Set pointers
	    ierr = InsertOffsetValues(ToRow, NumEntries, values, Indexor->PermuteOffsets()[i]);
	    if (ierr<0) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumPermuteIDs; i++) {
	    FromRow = A.GRID(PermuteFromLIDs[i]);
	    ToRow = GRID(PermuteToLIDs[i]);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, maxNumEntries, NumEntries, values, Indices)); // Set pointers
	    ierr = InsertGlobalValues(ToRow, NumEntries, values, Indices);
	    if (ierr<0) EPETRA_CHK_ERR(ierr);
          }
        }
      }
    }
    else { // A.IndicesAreGlobal()
      if (StaticGraph() || IndicesAreLocal()) {
        if(Indexor) {
          for (i=0; i<NumPermuteIDs; i++) {
	    FromRow = A.GRID(PermuteFromLIDs[i]);
	    ToRow = GRID(PermuteToLIDs[i]);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowView(FromRow, NumEntries, values, Indices)); // Set pointers
	    ierr = ReplaceOffsetValues(ToRow, NumEntries, values, Indexor->PermuteOffsets()[i]);
	    if (ierr<0) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumPermuteIDs; i++) {
	    FromRow = A.GRID(PermuteFromLIDs[i]);
	    ToRow = GRID(PermuteToLIDs[i]);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowView(FromRow, NumEntries, values, Indices)); // Set pointers
	    ierr = ReplaceGlobalValues(ToRow, NumEntries, values, Indices);
	    if (ierr<0) EPETRA_CHK_ERR(ierr);
          }
        }
      }
      else {
        if(Indexor) {
          for (i=0; i<NumPermuteIDs; i++) {
	    FromRow = A.GRID(PermuteFromLIDs[i]);
	    ToRow = GRID(PermuteToLIDs[i]);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowView(FromRow, NumEntries, values, Indices)); // Set pointers
	    ierr = InsertOffsetValues(ToRow, NumEntries, values, Indexor->PermuteOffsets()[i]);
	    if (ierr<0) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumPermuteIDs; i++) {
	    FromRow = A.GRID(PermuteFromLIDs[i]);
	    ToRow = GRID(PermuteToLIDs[i]);
	    EPETRA_CHK_ERR(A.ExtractGlobalRowView(FromRow, NumEntries, values, Indices)); // Set pointers
	    ierr = InsertGlobalValues(ToRow, NumEntries, values, Indices);
	    if (ierr<0) EPETRA_CHK_ERR(ierr);
          }
        }
      }
    }
  }

  if (maxNumEntries>0 && A.IndicesAreLocal() ) { // Delete Temp Space
    delete [] values;
    delete [] Indices;
  }
	
  return(0);
}

//=========================================================================
int Epetra_CrsMatrix::CopyAndPermuteRowMatrix(const Epetra_RowMatrix & A,
                                              int NumSameIDs, 
					      int NumPermuteIDs,
                                              int * PermuteToLIDs,
					      int *PermuteFromLIDs,
                                              const Epetra_OffsetIndex * Indexor ) {
  
  int i, j, ierr;
  
  int Row, NumEntries;
  int FromRow, ToRow;
  int maxNumEntries = A.MaxNumEntries();
  int * Indices = 0;
  double * values = 0;

  if (maxNumEntries>0) {
    Indices = new int[maxNumEntries];
    values = new double[maxNumEntries]; // Must extract values even though we discard them
  }
  
  const Epetra_Map & rowMap = A.RowMatrixRowMap();
  const Epetra_Map & colMap = A.RowMatrixColMap();

  // Do copy first
  if (NumSameIDs>0) {
    if (StaticGraph() || IndicesAreLocal()) {
      if( Indexor ) {
        for (i=0; i<NumSameIDs; i++) {
          Row = GRID(i);
          int AlocalRow = rowMap.LID(Row);
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(AlocalRow, maxNumEntries, NumEntries, values, Indices));
	  ierr = ReplaceOffsetValues(Row, NumEntries, values, Indexor->SameOffsets()[i]);
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
      else {
        for (i=0; i<NumSameIDs; i++) {
          Row = GRID(i);
          int AlocalRow = rowMap.LID(Row);
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(AlocalRow, maxNumEntries, NumEntries, values, Indices));
          for(j=0; j<NumEntries; ++j) {
            Indices[j] = LCID(colMap.GID(Indices[j]));
          }
	  ierr = ReplaceMyValues(i, NumEntries, values, Indices);
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
    }
    else {
      if( Indexor ) {
        for (i=0; i<NumSameIDs; i++) {
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(i, maxNumEntries, NumEntries, values, Indices));
          Row = GRID(i);
	  ierr = InsertOffsetValues(Row, NumEntries, values, Indexor->SameOffsets()[i]); 
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
      else {
        for (i=0; i<NumSameIDs; i++) {
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(i, maxNumEntries, NumEntries, values, Indices));
          Row = GRID(i);
          for( j=0; j<NumEntries; ++j ) Indices[j] = colMap.GID(Indices[j]); //convert to GIDs
	  ierr = InsertGlobalValues(Row, NumEntries, values, Indices); 
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
    }
  }
  
  // Do local permutation next
  if (NumPermuteIDs>0) {
    if (StaticGraph() || IndicesAreLocal()) {
      if( Indexor ) {
        for (i=0; i<NumPermuteIDs; i++) {
          FromRow = PermuteFromLIDs[i];
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(FromRow, maxNumEntries, NumEntries, values, Indices));
          ToRow = GRID(PermuteToLIDs[i]);
	  ierr = ReplaceOffsetValues(ToRow, NumEntries, values, Indexor->PermuteOffsets()[i]);
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
      else {
        for (i=0; i<NumPermuteIDs; i++) {
          FromRow = PermuteFromLIDs[i];
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(FromRow, maxNumEntries, NumEntries, values, Indices));
          ToRow = GRID(PermuteToLIDs[i]);
          for(j=0; j<NumEntries; ++j) {
            Indices[j] = LCID(colMap.GID(Indices[j]));
          }
	  ierr = ReplaceMyValues(ToRow, NumEntries, values, Indices);
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
    }
    else {
      if( Indexor ) {
        for (i=0; i<NumPermuteIDs; i++) {
          FromRow = PermuteFromLIDs[i];
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(FromRow, maxNumEntries, NumEntries, values, Indices));
          ToRow = GRID(PermuteToLIDs[i]);
	  ierr = InsertOffsetValues(ToRow, NumEntries, values, Indexor->PermuteOffsets()[i]); 
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
      else {
        for (i=0; i<NumPermuteIDs; i++) {
          FromRow = PermuteFromLIDs[i];
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(FromRow, maxNumEntries, NumEntries, values, Indices));
	  for (j=0; j<NumEntries; j++) Indices[j] = colMap.GID(Indices[j]); // convert to GIDs
          ToRow = GRID(PermuteToLIDs[i]);
	  ierr = InsertGlobalValues(ToRow, NumEntries, values, Indices); 
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
    }
  }	

  if (maxNumEntries>0) {
    delete [] values;
    delete [] Indices;
  }
  
  return(0);
}

//=========================================================================
int Epetra_CrsMatrix::PackAndPrepare(const Epetra_SrcDistObject & Source, 
				     int NumExportIDs,
                                     int * ExportLIDs,
				     int & LenExports,
                                     char *& Exports,
				     int & SizeOfPacket,
                                     int * Sizes,
                                     bool & VarSizes,
                                     Epetra_Distributor & Distor)
{
  (void)Distor;	
  // Rest of work can be done using RowMatrix only  
  const Epetra_RowMatrix & A = dynamic_cast<const Epetra_RowMatrix &>(Source);

  VarSizes = true; //enable variable block size data comm

  int TotalSendLength = 0;
  int * IntSizes = 0; 
  if( NumExportIDs>0 ) IntSizes = new int[NumExportIDs];

  for( int i = 0; i < NumExportIDs; ++i )
  {    
    int NumEntries;
    A.NumMyRowEntries( ExportLIDs[i], NumEntries );
    // Will have NumEntries doubles, NumEntries +2 ints, pack them interleaved     Sizes[i] = NumEntries;
    Sizes[i] = NumEntries;
    IntSizes[i] = 1 + (((NumEntries+2)*(int)sizeof(int))/(int)sizeof(double));
    TotalSendLength += (Sizes[i]+IntSizes[i]);
  }    
         
  double * DoubleExports = 0; 
  SizeOfPacket = (int)sizeof(double);
       
  //setup buffer locally for memory management by this object
  if( TotalSendLength*SizeOfPacket > LenExports )
  {
    if( LenExports > 0 ) delete [] Exports;
    LenExports = TotalSendLength*SizeOfPacket;
    DoubleExports = new double[TotalSendLength];
    for( int i = 0; i < TotalSendLength; ++i ) DoubleExports[i] = 0.0;
    Exports = (char *) DoubleExports;
  } 
 
  int NumEntries;
  int * Indices;
  double * values;
  int FromRow; 
  double * valptr, * dintptr; 
  int * intptr;                         
 
  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  NumEntries, Number of indices in row.
  // next NumEntries: The actual indices for the row.
 
  const Epetra_Map & rowMap = A.RowMatrixRowMap();
  const Epetra_Map & colMap = A.RowMatrixColMap();
 
  if( NumExportIDs > 0 )
  {
    int maxNumEntries = A.MaxNumEntries();
    dintptr = (double *) Exports;
    valptr = dintptr + IntSizes[0];
    intptr = (int *) dintptr;
    for (int i=0; i<NumExportIDs; i++)
    {
      FromRow = rowMap.GID(ExportLIDs[i]);
      intptr[0] = FromRow;
      values = valptr;
      Indices = intptr + 2;
      EPETRA_CHK_ERR(A.ExtractMyRowCopy(ExportLIDs[i], maxNumEntries, NumEntries, values, Indices));
      for (int j=0; j<NumEntries; j++) Indices[j] = colMap.GID(Indices[j]); // convert to GIDs
      intptr[1] = NumEntries; // Load second slot of segment
      if( i < (NumExportIDs-1) )
      {
        dintptr += (IntSizes[i]+Sizes[i]);
        valptr = dintptr + IntSizes[i+1];
        intptr = (int *) dintptr;
      }
    }
 
    for( int i = 0; i < NumExportIDs; ++i )
      Sizes[i] += IntSizes[i];
  }
 
  if( IntSizes ) delete [] IntSizes;
 
  return(0);
}

//=========================================================================
int Epetra_CrsMatrix::UnpackAndCombine(const Epetra_SrcDistObject & Source, 
				       int NumImportIDs,
                                       int * ImportLIDs, 
                                       int LenImports,
				       char * Imports,
                                       int & SizeOfPacket, 
				       Epetra_Distributor & Distor, 
				       Epetra_CombineMode CombineMode,
                                       const Epetra_OffsetIndex * Indexor )
{
  (void)Source;
  (void)LenImports;
  (void)SizeOfPacket;
  (void)Distor;
  if (NumImportIDs<=0) return(0);
	
  if (   CombineMode != Add
	 && CombineMode != Insert
	 && CombineMode != Zero )
    EPETRA_CHK_ERR(-1); //Unsupported CombineMode, defaults to Zero

  int NumEntries;
  int * Indices;
  double * values;
  int ToRow;
  int i, ierr;
  int IntSize;
  
  double * valptr, *dintptr;
  int * intptr;

  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  NumEntries, Number of indices in row.
  // next NumEntries: The actual indices for the row.

  dintptr = (double *) Imports;
  intptr = (int *) dintptr;
  NumEntries = intptr[1];
  IntSize = 1 + (((NumEntries+2)*(int)sizeof(int))/(int)sizeof(double));
  valptr = dintptr + IntSize;
 
  for (i=0; i<NumImportIDs; i++)
  {
    ToRow = GRID(ImportLIDs[i]);
    assert((intptr[0])==ToRow); // Sanity check
    values = valptr;
    Indices = intptr + 2;
 
    if (CombineMode==Add) {
      if (StaticGraph() || IndicesAreLocal()) {
        if( Indexor )
          ierr = SumIntoOffsetValues(ToRow, NumEntries, values, Indexor->RemoteOffsets()[i]);
        else
          ierr = SumIntoGlobalValues(ToRow, NumEntries, values, Indices);
      }
      else {
        if( Indexor )
          ierr = InsertOffsetValues(ToRow, NumEntries, values, Indexor->RemoteOffsets()[i]);
        else
          ierr = InsertGlobalValues(ToRow, NumEntries, values, Indices);
      }
      if (ierr<0) EPETRA_CHK_ERR(ierr);
    }
    else if (CombineMode==Insert) {
      if (StaticGraph() || IndicesAreLocal()) {
        if( Indexor )
          ierr = ReplaceOffsetValues(ToRow, NumEntries, values, Indexor->RemoteOffsets()[i]);
        else
          ierr = ReplaceGlobalValues(ToRow, NumEntries, values, Indices);
      }
      else {
        if( Indexor )
          ierr = InsertOffsetValues(ToRow, NumEntries, values, Indexor->RemoteOffsets()[i]);
        else
          ierr = InsertGlobalValues(ToRow, NumEntries, values, Indices);
      }
      if (ierr<0) EPETRA_CHK_ERR(ierr);
    }
 
    if( i < (NumImportIDs-1) )
    {
      dintptr += IntSize + NumEntries;
      intptr = (int *) dintptr;
      NumEntries = intptr[1];
      IntSize = 1 + (((NumEntries+2)*(int)sizeof(int))/(int)sizeof(double));
      valptr = dintptr + IntSize;
    }
  }

  return(0);
}

//=========================================================================

void Epetra_CrsMatrix::Print(ostream& os) const {
  int MyPID = RowMap().Comm().MyPID();
  int NumProc = RowMap().Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      /*      const Epetra_fmtflags olda = os.setf(ios::right,ios::adjustfield);
	      const Epetra_fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
	      const int             oldp = os.precision(12); */
      if (MyPID==0) {
	os <<  "\nNumber of Global Rows        = "; os << NumGlobalRows(); os << endl;
	os <<    "Number of Global Cols        = "; os << NumGlobalCols(); os << endl;
	os <<    "Number of Global Diagonals   = "; os << NumGlobalDiagonals(); os << endl;
	os <<    "Number of Global Nonzeros    = "; os << NumGlobalNonzeros(); os << endl;
	os <<    "Global Maximum Num Entries   = "; os << GlobalMaxNumEntries(); os << endl;
	if (LowerTriangular()) os <<    " ** Matrix is Lower Triangular **"; os << endl;
	if (UpperTriangular()) os <<    " ** Matrix is Upper Triangular **"; os << endl;
	if (NoDiagonal())      os <<    " ** Matrix has no diagonal     **"; os << endl; os << endl;
      }
			
      os <<  "\nNumber of My Rows        = "; os << NumMyRows(); os << endl;
      os <<    "Number of My Cols        = "; os << NumMyCols(); os << endl;
      os <<    "Number of My Diagonals   = "; os << NumMyDiagonals(); os << endl;
      os <<    "Number of My Nonzeros    = "; os << NumMyNonzeros(); os << endl;
      os <<    "My Maximum Num Entries   = "; os << MaxNumEntries(); os << endl; os << endl;

      os << flush;
      
      // Reset os flags
      
      /*      os.setf(olda,ios::adjustfield);
	      os.setf(oldf,ios::floatfield);
	      os.precision(oldp); */
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }
	
  {for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumMyRows1 = NumMyRows();
      int MaxNumIndices = MaxNumEntries();
      int * Indices  = new int[MaxNumIndices];
      double * values  = new double[MaxNumIndices];
      int NumIndices;
      int i, j;
			
      if (MyPID==0) {
	os.width(8);
	os <<  "   Processor ";
	os.width(10);
	os <<  "   Row Index ";
	os.width(10);
	os <<  "   Col Index ";
	os.width(20);
	os <<  "   Value     ";
	os << endl;
      }
      for (i=0; i<NumMyRows1; i++) {
	int Row = GRID(i); // Get global row number
	ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, values, Indices);
				
	for (j = 0; j < NumIndices ; j++) {   
	  os.width(8);
	  os <<  MyPID ; os << "    ";	
	  os.width(10);
	  os <<  Row ; os << "    ";	
	  os.width(10);
	  os <<  Indices[j]; os << "    ";
	  os.width(20);
	  os <<  values[j]; os << "    ";
	  os << endl;
	}
      }
			
      delete [] Indices;
      delete [] values;
      
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    RowMap().Comm().Barrier();
    RowMap().Comm().Barrier();
    RowMap().Comm().Barrier();
  }}
	
  return;
}
//=============================================================================
int Epetra_CrsMatrix::Multiply(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const {

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Multiply(TransA,x,y)");
#endif

  //
  // This function forms the product y = A * x or y = A' * x
  //

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  double* xp = (double*) x.Values();
  double* yp = (double*) y.Values();

  Epetra_Vector * xcopy = 0;
  if (&x==&y && Importer()==0 && Exporter()==0) {
    xcopy = new Epetra_Vector(x);
    xp = (double *) xcopy->Values();
  }
  UpdateImportVector(1); // Refresh import and output vectors if needed
  UpdateExportVector(1);

  if(!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if(Importer() != 0) {
      EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer(), Insert));
      xp = (double*) ImportVector_->Values();
    }
		
    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if(Exporter() != 0)  yp = (double*) ExportVector_->Values();
		
    // Do actual computation
    GeneralMV(xp, yp);

    if(Exporter() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*ExportVector_, *Exporter(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
  }
	
  else { // Transpose operation

    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    if(Exporter() != 0) {
      EPETRA_CHK_ERR(ExportVector_->Import(x, *Exporter(), Insert));
      xp = (double*) ExportVector_->Values();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if(Importer() != 0) yp = (double*) ImportVector_->Values();

    // Do actual computation
    GeneralMTV(xp, yp);

    if(Importer() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*ImportVector_, *Importer(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
  }


  UpdateFlops(2 * NumGlobalNonzeros());
  if (xcopy!=0) {
    delete xcopy;
    EPETRA_CHK_ERR(1); // Return positive code to alert the user about needing extra copy of x
    return(1);
  }
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Multiply(TransA,X,Y)");
#endif

#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  if(Epetra_CrsMatrixTraceDumpMultiply) {
    *out << std::boolalpha;
    *out << "\nEntering Epetra_CrsMatrix::Multipy("<<TransA<<",X,Y) ...\n";
    if(!TransA) {
      *out << "\nDomainMap =\n";
      this->DomainMap().Print(Teuchos::OSTab(out).o());
    }
    else {
      *out << "\nRangeMap =\n";
      this->RangeMap().Print(Teuchos::OSTab(out).o());
    }
    *out << "\nInitial input X with " << ( TransA ? "RangeMap" : "DomainMap" ) << " =\n\n";
    X.Print(Teuchos::OSTab(out).o());
  }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY

  //
  // This function forms the product Y = A * Y or Y = A' * X
  //
  if(!Filled()) {
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  }

  int NumVectors = X.NumVectors();
  if (NumVectors!=Y.NumVectors()) {
    EPETRA_CHK_ERR(-2); // Need same number of vectors in each MV
  }

  double** Xp = (double**) X.Pointers();
  double** Yp = (double**) Y.Pointers();

  int LDX = X.ConstantStride() ? X.Stride() : 0;
  int LDY = Y.ConstantStride() ? Y.Stride() : 0;

  Epetra_MultiVector * Xcopy = 0;
  if (&X==&Y && Importer()==0 && Exporter()==0) {
    Xcopy = new Epetra_MultiVector(X);
    Xp = (double **) Xcopy->Pointers();
    LDX = Xcopy->ConstantStride() ? Xcopy->Stride() : 0;
  }
  UpdateImportVector(NumVectors); // Make sure Import and Export Vectors are compatible
  UpdateExportVector(NumVectors);

  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer(), Insert));
      Xp = (double**)ImportVector_->Pointers();
      LDX = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
      if(Epetra_CrsMatrixTraceDumpMultiply) {
        *out << "\nColMap =\n";
        this->ColMap().Print(Teuchos::OSTab(out).o());
        *out << "\nX after import from DomainMap to ColMap =\n\n";
        ImportVector_->Print(Teuchos::OSTab(out).o());
      }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      Yp = (double**)ExportVector_->Pointers();
      LDY = ExportVector_->ConstantStride() ? ExportVector_->Stride() : 0;
    }

    // Do actual computation
    if (NumVectors==1)
      GeneralMV(*Xp, *Yp);
    else
      GeneralMM(Xp, LDX, Yp, LDY, NumVectors);
#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    if(Epetra_CrsMatrixTraceDumpMultiply) {
      *out << "\nRowMap =\n";
      this->RowMap().Print(Teuchos::OSTab(out).o());
      *out << "\nY after local mat-vec where Y has RowMap =\n\n";
      if(Exporter()!=0)
        ExportVector_->Print(Teuchos::OSTab(out).o());
      else
        Y.Print(Teuchos::OSTab(out).o());
    }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    if (Exporter()!=0) {
      Y.PutScalar(0.0);  // Make sure target is zero
      Y.Export(*ExportVector_, *Exporter(), Add); // Fill Y with Values from export vector
#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
      if(Epetra_CrsMatrixTraceDumpMultiply) {
        *out << "\nRangeMap =\n";
        this->RangeMap().Print(Teuchos::OSTab(out).o());
        *out << "\nY after export from RowMap to RangeMap = \n\n";
        Y.Print(Teuchos::OSTab(out).o());
      }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Y.Reduce());
  }
  else { // Transpose operation

    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors

    if (Exporter()!=0) {
      EPETRA_CHK_ERR(ExportVector_->Import(X, *Exporter(), Insert));
      Xp = (double**)ExportVector_->Pointers();
      LDX = ExportVector_->ConstantStride() ? ExportVector_->Stride() : 0;
#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
      if(Epetra_CrsMatrixTraceDumpMultiply) {
        *out << "\nRowMap =\n";
        this->RowMap().Print(Teuchos::OSTab(out).o());
        *out << "\nX after import from RangeMap to RowMap =\n\n";
        ExportVector_->Print(Teuchos::OSTab(out).o());
      }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      Yp = (double**)ImportVector_->Pointers();
      LDY = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
    }

    // Do actual computation
    if (NumVectors==1)
      GeneralMTV(*Xp, *Yp);
    else
      GeneralMTM(Xp, LDX, Yp, LDY, NumVectors);
#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    if(Epetra_CrsMatrixTraceDumpMultiply) {
      *out << "\nColMap =\n";
      this->ColMap().Print(Teuchos::OSTab(out).o());
      *out << "\nY after local transpose mat-vec where Y has ColMap =\n\n";
      if(Importer()!=0)
        ImportVector_->Print(Teuchos::OSTab(out).o());
      else
        Y.Print(Teuchos::OSTab(out).o());
    }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    if (Importer()!=0) {
      Y.PutScalar(0.0);  // Make sure target is zero
      EPETRA_CHK_ERR(Y.Export(*ImportVector_, *Importer(), Add)); // Fill Y with Values from export vector
#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
      if(Epetra_CrsMatrixTraceDumpMultiply) {
        *out << "\nDomainMap =\n";
        this->DomainMap().Print(Teuchos::OSTab(out).o());
        *out << "\nY after export from ColMap to DomainMap =\n\n";
        Y.Print(Teuchos::OSTab(out).o());
      }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1)  EPETRA_CHK_ERR(Y.Reduce());
  }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  if (Xcopy!=0) {
    delete Xcopy;
    EPETRA_CHK_ERR(1); // Return positive code to alert the user about needing extra copy of X
    return(1);
  }

#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
  if(Epetra_CrsMatrixTraceDumpMultiply) {
    *out << "\nFinal output Y is the last Y printed above!\n";
    *out << "\nLeaving Epetra_CrsMatrix::Multipy("<<TransA<<",X,Y) ...\n";
  }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY

  return(0);
}
//=======================================================================================================
void Epetra_CrsMatrix::UpdateImportVector(int NumVectors) const {    
  if(Importer() != 0) {
    if(ImportVector_ != 0) {
      if(ImportVector_->NumVectors() != NumVectors) { 
	delete ImportVector_; 
	ImportVector_= 0;
      }
    }
    if(ImportVector_ == 0) 
      ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors); // Create import vector if needed
  }
  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::UpdateExportVector(int NumVectors) const {    
  if(Exporter() != 0) {
    if(ExportVector_ != 0) {
      if(ExportVector_->NumVectors() != NumVectors) { 
	delete ExportVector_; 
	ExportVector_= 0;
      }
    }
    if(ExportVector_ == 0) 
      ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors); // Create Export vector if needed
  }
  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::GeneralMV(double * x, double * y)  const {
  
if (StorageOptimized() && Graph().StorageOptimized()) {

  double * values = All_Values();
  int * Indices = Graph().All_Indices();
  int * IndexOffset = Graph().IndexOffset();
#ifdef EPETRA_HAVE_OMP
  const int numMyRows = NumMyRows_;
#pragma omp parallel for default(none) shared(IndexOffset,values,Indices,y,x)
     for (int row=0; row<numMyRows; ++row)
        {
 	  const int curOffset = IndexOffset[row];
          const double *val_ptr    = values+curOffset;
          const int    *colnum_ptr = Indices+curOffset;
 	  double s = 0.;
 	  const double *const val_end_of_row = &values[IndexOffset[row+1]];
 	  while (val_ptr != val_end_of_row)
 	    s += *val_ptr++ * x[*colnum_ptr++];
 	  y[row] = s;
	}
#else
#ifndef Epetra_ENABLE_CASK
       const double *val_ptr    = values;
       const int    *colnum_ptr = Indices;
       double       * dst_ptr = y;
       for (int row=0; row<NumMyRows_; ++row)
 	{
 	  double s = 0.;
 	  const double *const val_end_of_row = &values[IndexOffset[row+1]];
 	  while (val_ptr != val_end_of_row)
 	    s += *val_ptr++ * x[*colnum_ptr++];
 	  *dst_ptr++ = s;
 	}
#else
       cask_csr_dax_new(NumMyRows_, IndexOffset, Indices,
                        values, x, y, cask);
#endif // Epetra_ENABLE_CASK
#endif // EPETRA_HAVE_OMP

    return;
  }
  else if (!StorageOptimized() && !Graph().StorageOptimized()) {


    int* NumEntriesPerRow = Graph().NumIndicesPerRow();
    int** Indices = Graph().Indices();
    double** srcValues = Values();
        const int numMyRows = NumMyRows_;

    // Do actual computation
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(NumEntriesPerRow,Indices,srcValues,y,x)
#endif
    for(int i = 0; i < numMyRows; i++) {
      int     NumEntries = NumEntriesPerRow[i];
      int*    RowIndices = Indices[i];
      double* RowValues  = srcValues[i];
      double sum = 0.0;
      for(int j = 0; j < NumEntries; j++) 
	sum += *RowValues++ * x[*RowIndices++];
      
      y[i] = sum;
      
    }
  } 
  else { // Case where StorageOptimized is incompatible:  Use general accessors.

    const int numMyRows = NumMyRows_;

    // Do actual computation
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(x,y)
#endif
    for(int i = 0; i < numMyRows; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    RowIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      double sum = 0.0;
      for(int j = 0; j < NumEntries; j++) 
	sum += *RowValues++ * x[*RowIndices++];
      
      y[i] = sum;
      
    }
  } 
  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::GeneralMTV(double * x, double * y) const {

  int NumCols = NumMyCols();
#if !defined(FORTRAN_DISABLED) || defined(Epetra_ENABLE_CASK)
  if (StorageOptimized() && Graph().StorageOptimized()) {
    double * values = All_Values_;
    int * Indices = Graph().All_Indices();
    int * IndexOffset = Graph().IndexOffset();
#ifndef Epetra_ENABLE_CASK
   int ione = 1;
   EPETRA_DCRSMV_F77(&ione, &NumMyRows_, &NumCols, values, Indices, IndexOffset, x, y);
#else
   cask_csr_datx( NumMyRows_, NumCols, IndexOffset,  Indices,  values,x ,y );
#endif

    return;
  }
#endif // FORTRAN_DISABLED
  for(int i = 0; i < NumCols; i++) 
    y[i] = 0.0; // Initialize y for transpose multiply

  if (StorageOptimized() && Graph().StorageOptimized()) {
    double * values = All_Values_;
    int * Indices = Graph().All_Indices();
    int * IndexOffset = Graph().IndexOffset();
    for(int i = 0; i < NumMyRows_; ++i) {
      int prevOffset = *IndexOffset++;
      int NumEntries = *IndexOffset - prevOffset;
      double xi = x[i];
      for(int j = 0; j < NumEntries; j++) 
	y[*Indices++] += *values++ * xi;
    }
  }
  else if (!StorageOptimized() && !Graph().StorageOptimized()) {

    int* NumEntriesPerRow = Graph().NumIndicesPerRow();
    int** Indices = Graph().Indices();
    double** srcValues = Values();
    
    for(int i = 0; i < NumMyRows_; i++) {
      int     NumEntries = *NumEntriesPerRow++;
      int*    RowIndices = *Indices++;
      double* RowValues  = *srcValues++;
      double xi = x[i];
      for(int j = 0; j < NumEntries; j++) 
	y[*RowIndices++] += *RowValues++ * xi;
    }
  }
  else { // Case where StorageOptimized is incompatible:  Use general accessors.
  
    for(int i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    RowIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      double xi = x[i];
      for(int j = 0; j < NumEntries; j++) 
	y[*RowIndices++] += *RowValues++ * xi;
    }
  }

  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::GeneralMM(double ** X, int LDX, double ** Y, int LDY, int NumVectors) const {

#if !defined(FORTRAN_DISABLED) || defined(Epetra_ENABLE_CASK)
  if (StorageOptimized() && Graph().StorageOptimized()) {
    double * values = All_Values_;
    int * Indices = Graph().All_Indices();
    int * IndexOffset = Graph().IndexOffset();

    if (LDX!=0 && LDY!=0) {
#ifndef Epetra_ENABLE_CASK
    int izero = 0;
    EPETRA_DCRSMM_F77(&izero, &NumMyRows_, &NumMyRows_, values, Indices, IndexOffset, *X, &LDX, *Y, &LDY, &NumVectors);
#else
    cask_csr_dgesmm_new(0, 1.0, NumMyRows_, NumMyRows_,  NumVectors,
                    IndexOffset, Indices, values, *X, LDX, 0.0,  *Y, LDY,cask);
#endif
    return;
    }

    double ** const xp = X;
    double ** const yp = Y;
    const int numMyRows = NumMyRows_;
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(IndexOffset,Indices,values,NumVectors)
#endif
    for (int i=0; i < numMyRows; i++) {
      int prevOffset = IndexOffset[i];
      int NumEntries = IndexOffset[i+1] - prevOffset;
      int *    RowIndices = Indices+prevOffset;
      double * RowValues  = values+prevOffset;
      for (int k=0; k<NumVectors; k++) {
	double sum = 0.0;
	const double * const x = xp[k];
	double * const y = yp[k];
	for (int j=0; j < NumEntries; j++) sum += RowValues[j] * x[RowIndices[j]];
	y[i] = sum;
      }
    }
  }
  else if (!StorageOptimized() && !Graph().StorageOptimized()) {
#else
  if (!StorageOptimized() && !Graph().StorageOptimized()) {
#endif

    int* NumEntriesPerRow = Graph().NumIndicesPerRow();
    int** Indices = Graph().Indices();
    double** srcValues = Values();
    double ** const xp = X;
    double ** const yp = Y;
    const int numMyRows = NumMyRows_;

#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(NumEntriesPerRow,Indices,srcValues,NumVectors)
#endif
    for (int i=0; i < numMyRows; i++) {
      int      NumEntries = NumEntriesPerRow[i];
      int *    RowIndices = Indices[i];
      double * RowValues  = srcValues[i];
      for (int k=0; k<NumVectors; k++) {
	double sum = 0.0;
	const double * const x = xp[k];
	double * const y = yp[k];
	for (int j=0; j < NumEntries; j++) sum += RowValues[j] * x[RowIndices[j]];
	y[i] = sum;
      }
    }
  }
  else {

    double ** const xp = X;
    double ** const yp = Y;
    const int numMyRows = NumMyRows_;
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(NumVectors)
#endif
    for (int i=0; i < numMyRows; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    RowIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for (int k=0; k<NumVectors; k++) {
	double sum = 0.0;
	double * x = xp[k];
	for (int j=0; j < NumEntries; j++) sum += RowValues[j] * x[RowIndices[j]];
	yp[k][i] = sum;
      }
    }
  }
  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::GeneralMTM(double ** X, int LDX, double ** Y, int LDY, int NumVectors)  const{

  int NumCols = NumMyCols();
#if !defined(FORTRAN_DISABLED) || defined(Epetra_ENABLE_CASK)
  if (StorageOptimized() && Graph().StorageOptimized()) {
    if (LDX!=0 && LDY!=0) {
      double * values = All_Values_;
      int * Indices = Graph().All_Indices();
      int * IndexOffset = Graph().IndexOffset();

#ifndef Epetra_ENABLE_CASK
      int ione = 1;
      EPETRA_DCRSMM_F77(&ione, &NumMyRows_, &NumCols, values, Indices, IndexOffset, *X, &LDX, *Y, &LDY, &NumVectors);
#else
      cask_csr_dgesmm_new(1, 1.0, NumMyRows_, NumCols,  NumVectors,
                          IndexOffset, Indices, values, *X, LDX, 0.0,
                          *Y, LDY, cask);
#endif
      return;
    }
  }
#endif
  for (int k=0; k<NumVectors; k++) 
    for (int i=0; i < NumCols; i++) 
      Y[k][i] = 0.0; // Initialize y for transpose multiply
  
  if (StorageOptimized() && Graph().StorageOptimized()) {
    double * values = All_Values_;
    int * Indices = Graph().All_Indices();
    int * IndexOffset = Graph().IndexOffset();
    for (int i=0; i < NumMyRows_; i++) {
      int prevOffset = *IndexOffset++;
      int NumEntries = *IndexOffset - prevOffset;
      int *    RowIndices = Indices+prevOffset;
      double * RowValues  = values+prevOffset;
      
      for (int k=0; k<NumVectors; k++) {
	double * y = Y[k];
	double * x = X[k];
	for (int j=0; j < NumEntries; j++) 
	  y[RowIndices[j]] += RowValues[j] * x[i];
      }
    }
  }
  else if (!StorageOptimized() && !Graph().StorageOptimized()) {
    
    int* NumEntriesPerRow = Graph().NumIndicesPerRow();
    int** Indices = Graph().Indices();
    double** srcValues = Values();

    for (int i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *srcValues++;
      for (int k=0; k<NumVectors; k++) {
	double * y = Y[k];
	double * x = X[k];
	for (int j=0; j < NumEntries; j++) 
	  y[RowIndices[j]] += RowValues[j] * x[i];
      }
    }
  }
  else { // Case where StorageOptimized is incompatible:  Use general accessors.
    
    for (int i=0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    RowIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for (int k=0; k<NumVectors; k++) {
	double * y = Y[k];
	double * x = X[k];
	for (int j=0; j < NumEntries; j++) 
	  y[RowIndices[j]] += RowValues[j] * x[i];
      }
    }
  }
  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::GeneralSV(bool Upper, bool Trans, bool UnitDiagonal, double * xp, double * yp)  const {


  int i, j, j0;

#if !defined(FORTRAN_DISABLED) || defined(Epetra_ENABLE_CASK)
  if (StorageOptimized() && Graph().StorageOptimized() && ((UnitDiagonal && NoDiagonal())|| (!UnitDiagonal && !NoDiagonal()))) {
    double * values = All_Values();
    int * Indices = Graph().All_Indices();
    int * IndexOffset = Graph().IndexOffset();

    int iupper = Upper ? 1:0;
    int itrans = Trans ? 1:0;
    int udiag =  UnitDiagonal ? 1:0;
    int nodiag = NoDiagonal() ? 1:0;
    int xysame = (xp==yp) ? 1:0;

#ifndef Epetra_ENABLE_CASK
    EPETRA_DCRSSV_F77( &iupper, &itrans, &udiag, &nodiag, &NumMyRows_, &NumMyRows_, values, Indices, IndexOffset, xp, yp, &xysame);
#else
    cask_csr_dtrsv_new( iupper, itrans, udiag, nodiag, 0, xysame, NumMyRows_,
                    NumMyRows_, IndexOffset, Indices, values, xp, yp, cask);
#endif
    return;
  }
  //=================================================================
  else { // !StorageOptimized()
  //=================================================================
#endif    
    if (!Trans) {
      
      if (Upper) {
	
	j0 = 1;
	if (NoDiagonal()) 
	  j0--; // Include first term if no diagonal
	for (i=NumMyRows_-1; i >=0; i--) {
	  int      NumEntries = NumMyEntries(i);
	  int *    RowIndices = Graph().Indices(i);
	  double * RowValues  = Values(i);
	  double sum = 0.0;
	  for (j=j0; j < NumEntries; j++) 
	    sum += RowValues[j] * yp[RowIndices[j]];
	  
	  if (UnitDiagonal) 
	    yp[i] = xp[i] - sum;
	  else 
	    yp[i] = (xp[i] - sum)/RowValues[0];
	  
	}
      }
      else {
	j0 = 1;
	if (NoDiagonal())
	  j0--; // Include first term if no diagonal
	for (i=0; i < NumMyRows_; i++) {
	  int      NumEntries = NumMyEntries(i) - j0;
	  int *    RowIndices = Graph().Indices(i);
	  double * RowValues  = Values(i);
	  double sum = 0.0;
	  for (j=0; j < NumEntries; j++) 
	    sum += RowValues[j] * yp[RowIndices[j]];
	  
	  if (UnitDiagonal) 
	    yp[i] = xp[i] - sum;
	  else 
	    yp[i] = (xp[i] - sum)/RowValues[NumEntries];
	  
	}
      }
    }
    
    // ***********  Transpose case *******************************
    
    else {
      
      if (xp!=yp) 
	for (i=0; i < NumMyRows_; i++) 
	  yp[i] = xp[i]; // Initialize y for transpose solve
      
      if (Upper) {
	
	j0 = 1;
	if (NoDiagonal()) 
	  j0--; // Include first term if no diagonal
	
	for (i=0; i < NumMyRows_; i++) {
	  int      NumEntries = NumMyEntries(i);
	  int *    RowIndices = Graph().Indices(i);
	  double * RowValues  = Values(i);
	  if (!UnitDiagonal) 
	    yp[i] = yp[i]/RowValues[0];
	  double ytmp = yp[i];
	  for (j=j0; j < NumEntries; j++) 
	    yp[RowIndices[j]] -= RowValues[j] * ytmp;
	}
      }
      else {
	
	j0 = 1;
	if (NoDiagonal()) 
	  j0--; // Include first term if no diagonal
	
	for (i=NumMyRows_-1; i >= 0; i--) {
	  int      NumEntries = NumMyEntries(i) - j0;
	  int *    RowIndices = Graph().Indices(i);
	  double * RowValues  = Values(i);
	  if (!UnitDiagonal) 
	    yp[i] = yp[i]/RowValues[NumEntries];
	  double ytmp = yp[i];
	  for (j=0; j < NumEntries; j++) 
	    yp[RowIndices[j]] -= RowValues[j] * ytmp;
	}
      }
      
    }
#if !defined(FORTRAN_DISABLED) || defined(Epetra_ENABLE_CASK)
  }
#endif
  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::GeneralSM(bool Upper, bool Trans, bool UnitDiagonal, double ** Xp, int LDX, double ** Yp, int LDY, int NumVectors)  const{

  int i, j, j0, k;
  double diag = 0.0;

  if (StorageOptimized() && Graph().StorageOptimized()) {
    double * values = All_Values();
    int * Indices = Graph().All_Indices();
    int * IndexOffset = Graph().IndexOffset();
#if !defined(FORTRAN_DISABLED) || defined(Epetra_ENABLE_CASK)
    if (LDX!=0 && LDY!=0 && ((UnitDiagonal && NoDiagonal()) || (!UnitDiagonal && !NoDiagonal()))) {
      int iupper = Upper ? 1:0;
      int itrans = Trans ? 1:0;
      int udiag =  UnitDiagonal ? 1:0;
      int nodiag = NoDiagonal() ? 1:0;
      int xysame = (Xp==Yp) ? 1:0;

#ifndef Epetra_ENABLE_CASK
      EPETRA_DCRSSM_F77( &iupper, &itrans, &udiag, &nodiag, &NumMyRows_, &NumMyRows_, values, Indices, IndexOffset, 
			 *Xp, &LDX, *Yp, &LDY, &xysame, &NumVectors);
#else
      cask_csr_dtrsm( iupper, itrans, udiag, nodiag, 0, xysame,  NumMyRows_, 
                      NumMyRows_, NumVectors, IndexOffset, Indices, values, 
                      *Xp, LDX, *Yp, LDY);
#endif
      return;
    }
#endif
    if(!Trans) {   
      if(Upper) {   
	j0 = 1;
	if(NoDiagonal()) 
	  j0--; // Include first term if no diagonal
	for(i = NumMyRows_ - 1; i >= 0; i--) {
	  int Offset = IndexOffset[i];
	  int      NumEntries = IndexOffset[i+1]-Offset;
	  int *    RowIndices = Indices+Offset;
	  double * RowValues  = values+Offset;
	  if(!UnitDiagonal) 
	    diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
	  for(k = 0; k < NumVectors; k++) {
	    double sum = 0.0;
	    for(j = j0; j < NumEntries; j++) 
	      sum += RowValues[j] * Yp[k][RowIndices[j]];
					
	    if(UnitDiagonal) 
	      Yp[k][i] = Xp[k][i] - sum;
	    else 
	      Yp[k][i] = (Xp[k][i] - sum) * diag;
	  }
	}
      }
      else {
	j0 = 1;
	if(NoDiagonal()) 
	  j0--; // Include first term if no diagonal
	for(i = 0; i < NumMyRows_; i++) {
	  int Offset = IndexOffset[i];
	  int      NumEntries = IndexOffset[i+1]-Offset - j0;
	  int *    RowIndices = Indices+Offset;
	  double * RowValues  = values+Offset;
	  if(!UnitDiagonal)
	    diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
	  for(k = 0; k < NumVectors; k++) {
	    double sum = 0.0;
	    for(j = 0; j < NumEntries; j++) 
	      sum += RowValues[j] * Yp[k][RowIndices[j]];
					
	    if(UnitDiagonal) 
	      Yp[k][i] = Xp[k][i] - sum;
	    else 
	      Yp[k][i] = (Xp[k][i] - sum)*diag;
	  }
	}
      }
    }
    // ***********  Transpose case *******************************

    else {
      for(k = 0; k < NumVectors; k++) 
	if(Yp[k] != Xp[k]) 
	  for(i = 0; i < NumMyRows_; i++)
	    Yp[k][i] = Xp[k][i]; // Initialize y for transpose multiply
    
      if(Upper) {
	j0 = 1;
	if(NoDiagonal()) 
	  j0--; // Include first term if no diagonal
      
	for(i = 0; i < NumMyRows_; i++) {
	  int Offset = IndexOffset[i];
	  int      NumEntries = IndexOffset[i+1]-Offset;
	  int *    RowIndices = Indices+Offset;
	  double * RowValues  = values+Offset;
	  if(!UnitDiagonal) 
	    diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
	  for(k = 0; k < NumVectors; k++) {
	    if(!UnitDiagonal) 
	      Yp[k][i] = Yp[k][i]*diag;
	    double ytmp = Yp[k][i];
	    for(j = j0; j < NumEntries; j++) 
	      Yp[k][RowIndices[j]] -= RowValues[j] * ytmp;
	  }
	}
      }
      else {
	j0 = 1;
	if(NoDiagonal()) 
	  j0--; // Include first term if no diagonal  
	for(i = NumMyRows_ - 1; i >= 0; i--) {
	  int Offset = IndexOffset[i];
	  int      NumEntries = IndexOffset[i+1]-Offset - j0;
	  int *    RowIndices = Indices+Offset;
	  double * RowValues  = values+Offset;
	  if(!UnitDiagonal) 
	    diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
	  for(k = 0; k < NumVectors; k++) {
	    if(!UnitDiagonal)  
	      Yp[k][i] = Yp[k][i]*diag;
	    double ytmp = Yp[k][i];
	    for(j = 0; j < NumEntries; j++)
	      Yp[k][RowIndices[j]] -= RowValues[j] * ytmp;
	  }
	}
      }
    }
  }
    // ========================================================
  else { // !StorageOptimized()
    // ========================================================

    if(!Trans) {   
      if(Upper) {   
	j0 = 1;
	if(NoDiagonal()) 
	  j0--; // Include first term if no diagonal
	for(i = NumMyRows_ - 1; i >= 0; i--) {
	  int     NumEntries = NumMyEntries(i);
	  int*    RowIndices = Graph().Indices(i);
	  double* RowValues  = Values(i);
	  if(!UnitDiagonal) 
	    diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
	  for(k = 0; k < NumVectors; k++) {
	    double sum = 0.0;
	    for(j = j0; j < NumEntries; j++) 
	      sum += RowValues[j] * Yp[k][RowIndices[j]];
					
	    if(UnitDiagonal) 
	      Yp[k][i] = Xp[k][i] - sum;
	    else 
	      Yp[k][i] = (Xp[k][i] - sum) * diag;
	  }
	}
      }
      else {
	j0 = 1;
	if(NoDiagonal()) 
	  j0--; // Include first term if no diagonal
	for(i = 0; i < NumMyRows_; i++) {
	  int     NumEntries = NumMyEntries(i) - j0;
	  int*    RowIndices = Graph().Indices(i);
	  double* RowValues  = Values(i);
	  if(!UnitDiagonal)
	    diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
	  for(k = 0; k < NumVectors; k++) {
	    double sum = 0.0;
	    for(j = 0; j < NumEntries; j++) 
	      sum += RowValues[j] * Yp[k][RowIndices[j]];
					
	    if(UnitDiagonal) 
	      Yp[k][i] = Xp[k][i] - sum;
	    else 
	      Yp[k][i] = (Xp[k][i] - sum)*diag;
	  }
	}
      }
    }
    // ***********  Transpose case *******************************

    else {
      for(k = 0; k < NumVectors; k++) 
	if(Yp[k] != Xp[k]) 
	  for(i = 0; i < NumMyRows_; i++)
	    Yp[k][i] = Xp[k][i]; // Initialize y for transpose multiply
    
      if(Upper) {
	j0 = 1;
	if(NoDiagonal()) 
	  j0--; // Include first term if no diagonal
      
	for(i = 0; i < NumMyRows_; i++) {
	  int     NumEntries = NumMyEntries(i);
	  int*    RowIndices = Graph().Indices(i);
	  double* RowValues  = Values(i);
	  if(!UnitDiagonal) 
	    diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
	  for(k = 0; k < NumVectors; k++) {
	    if(!UnitDiagonal) 
	      Yp[k][i] = Yp[k][i]*diag;
	    double ytmp = Yp[k][i];
	    for(j = j0; j < NumEntries; j++) 
	      Yp[k][RowIndices[j]] -= RowValues[j] * ytmp;
	  }
	}
      }
      else {
	j0 = 1;
	if(NoDiagonal()) 
	  j0--; // Include first term if no diagonal  
	for(i = NumMyRows_ - 1; i >= 0; i--) {
	  int     NumEntries = NumMyEntries(i) - j0;
	  int*    RowIndices = Graph().Indices(i);
	  double* RowValues  = Values(i);
	  if(!UnitDiagonal) 
	    diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
	  for(k = 0; k < NumVectors; k++) {
	    if(!UnitDiagonal)  
	      Yp[k][i] = Yp[k][i]*diag;
	    double ytmp = Yp[k][i];
	    for(j = 0; j < NumEntries; j++)
	      Yp[k][RowIndices[j]] -= RowValues[j] * ytmp;
	  }
	}
      }
    }
  }
  
  return;
}
//=============================================================================
int Epetra_CrsMatrix::Multiply1(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const {

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Multiply1(TransA,x,y)");
#endif

  //
  // This function forms the product y = A * x or y = A' * x
  //

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  int i, j;
  double* xp = (double*) x.Values();
  double* yp = (double*) y.Values();
  int NumMyCols_ = NumMyCols();

  if(!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if(Importer() != 0) {
      if(ImportVector_ != 0) {
	if(ImportVector_->NumVectors() != 1) { 
	  delete ImportVector_; 
	  ImportVector_= 0;
	}
      }
      if(ImportVector_ == 0) 
	ImportVector_ = new Epetra_MultiVector(ColMap(),1); // Create import vector if needed
      EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer(), Insert));
      xp = (double*) ImportVector_->Values();
    }
		
    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if(Exporter() != 0) {
      if(ExportVector_ != 0) {
	if(ExportVector_->NumVectors() != 1) { 
	  delete ExportVector_; 
	  ExportVector_= 0;
	}
      }
      if(ExportVector_ == 0) 
	ExportVector_ = new Epetra_MultiVector(RowMap(),1); // Create Export vector if needed
      yp = (double*) ExportVector_->Values();
    }
		
    // Do actual computation
    for(i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    RowIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      double sum = 0.0;
      for(j = 0; j < NumEntries; j++) 
	sum += RowValues[j] * xp[RowIndices[j]];
			
      yp[i] = sum;
			
    }
    if(Exporter() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*ExportVector_, *Exporter(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
  }
	
  else { // Transpose operation

    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    if(Exporter() != 0) {
      if(ExportVector_ != 0) {
	if(ExportVector_->NumVectors() != 1) { 
	  delete ExportVector_; 
	  ExportVector_= 0;
	}
      }
      if(ExportVector_ == 0) 
	ExportVector_ = new Epetra_MultiVector(RowMap(),1); // Create Export vector if needed
      EPETRA_CHK_ERR(ExportVector_->Import(x, *Exporter(), Insert));
      xp = (double*) ExportVector_->Values();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if(Importer() != 0) {
      if(ImportVector_ != 0) {
	if(ImportVector_->NumVectors() != 1) { 
	  delete ImportVector_; 
	  ImportVector_= 0;
	}
      }
      if(ImportVector_ == 0) 
	ImportVector_ = new Epetra_MultiVector(ColMap(),1); // Create import vector if needed
      yp = (double*) ImportVector_->Values();
    }

    // Do actual computation
    for(i = 0; i < NumMyCols_; i++) 
      yp[i] = 0.0; // Initialize y for transpose multiply
        
    for(i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    RowIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for(j = 0; j < NumEntries; j++) 
	yp[RowIndices[j]] += RowValues[j] * xp[i];
    }
    if(Importer() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*ImportVector_, *Importer(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
  }

  UpdateFlops(2 * NumGlobalNonzeros());
  return(0);
}
//=============================================================================
int Epetra_CrsMatrix::Multiply1(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Multiply1(TransA,X,Y)");
#endif

  //
  // This function forms the product Y = A * Y or Y = A' * X
  //
  if((X.NumVectors() == 1) && (Y.NumVectors() == 1)) {
    double* xp = (double*) X[0];
    double* yp = (double*) Y[0];
    Epetra_Vector x(View, X.Map(), xp);
    Epetra_Vector y(View, Y.Map(), yp);
    EPETRA_CHK_ERR(Multiply1(TransA, x, y));
    return(0);
  }
  if(!Filled()) {
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  }

  int i, j, k;

  double** Xp = (double**) X.Pointers();
  double** Yp = (double**) Y.Pointers();

  int NumVectors = X.NumVectors();
  int NumMyCols_ = NumMyCols();


  // Need to better manage the Import and Export vectors:
  // - Need accessor functions
  // - Need to make the NumVector match (use a View to do this)
  // - Need to look at RightScale and ColSum routines too.

  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()!=NumVectors) { 
	  delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) 
	ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors); // Create import vector if needed
      EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer(), Insert));
      Xp = (double**)ImportVector_->Pointers();
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()!=NumVectors) { 
	  delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) 
	ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors); // Create Export vector if needed
      Yp = (double**)ExportVector_->Pointers();
    }

    // Do actual computation

    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = NumMyEntries(i);
      int *    RowIndices = Graph().Indices(i);
      double * RowValues  = Values(i);
      for (k=0; k<NumVectors; k++) {
	double sum = 0.0;
	for (j=0; j < NumEntries; j++) sum += RowValues[j] * Xp[k][RowIndices[j]];
	Yp[k][i] = sum;
      }
    }
    if (Exporter()!=0) {
      Y.PutScalar(0.0); // Make sure target is zero
      Y.Export(*ExportVector_, *Exporter(), Add); // Fill Y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Y.Reduce());
  }
  else { // Transpose operation
		

    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors

    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()!=NumVectors) { 
	  delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) 
	ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors); // Create Export vector if needed
      EPETRA_CHK_ERR(ExportVector_->Import(X, *Exporter(), Insert));
      Xp = (double**)ExportVector_->Pointers();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()!=NumVectors) { 
	  delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) 
	ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors); // Create import vector if needed
      Yp = (double**)ImportVector_->Pointers();
    }

    // Do actual computation



    for (k=0; k<NumVectors; k++) 
      for (i=0; i < NumMyCols_; i++) 
	Yp[k][i] = 0.0; // Initialize y for transpose multiply
    
    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = NumMyEntries(i);
      int *    RowIndices = Graph().Indices(i);
      double * RowValues  = Values(i);
      for (k=0; k<NumVectors; k++) {
	for (j=0; j < NumEntries; j++) 
	  Yp[k][RowIndices[j]] += RowValues[j] * Xp[k][i];
      }
    }
    if (Importer()!=0) {
      Y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(Y.Export(*ImportVector_, *Importer(), Add)); // Fill Y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1)  EPETRA_CHK_ERR(Y.Reduce());
  }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Solve1(bool Upper, bool Trans, bool UnitDiagonal,
			    const Epetra_Vector& x, Epetra_Vector& y) const
{

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Solve1(Upper,Trans,UnitDiag,x,y)");
#endif

  //
  // This function finds y such that Ly = x or Uy = x or the transpose cases.
  //

  if (!Filled()) {
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  }

  if ((Upper) && (!UpperTriangular())) 
    EPETRA_CHK_ERR(-2);
  if ((!Upper) && (!LowerTriangular())) 
    EPETRA_CHK_ERR(-3);
  if ((!UnitDiagonal) && (NoDiagonal())) 
    EPETRA_CHK_ERR(-4); // If matrix has no diagonal, we must use UnitDiagonal
  if ((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_)) 
    EPETRA_CHK_ERR(-5); // Need each row to have a diagonal
      

  int i, j, j0;
  int * NumEntriesPerRow = Graph().NumIndicesPerRow();
  int ** Indices = Graph().Indices();
  double ** Vals = Values();
  int NumMyCols_ = NumMyCols();

  // If upper, point to last row
  if ((Upper && !Trans) || (!Upper && Trans)) {
    NumEntriesPerRow += NumMyRows_-1;
    Indices += NumMyRows_-1;
    Vals += NumMyRows_-1;
  }
    
  double *xp = (double*)x.Values();
  double *yp = (double*)y.Values();

  if (!Trans) {

    if (Upper) {

      j0 = 1;
      if (NoDiagonal()) 
	j0--; // Include first term if no diagonal
      for (i=NumMyRows_-1; i >=0; i--) {
	int      NumEntries = *NumEntriesPerRow--;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Vals--;
	double sum = 0.0;
	for (j=j0; j < NumEntries; j++) 
	  sum += RowValues[j] * yp[RowIndices[j]];
				
	if (UnitDiagonal) 
	  yp[i] = xp[i] - sum;
	else 
	  yp[i] = (xp[i] - sum)/RowValues[0];
				
      }
    }
    else {
      j0 = 1;
      if (NoDiagonal())
	j0--; // Include first term if no diagonal
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++ - j0;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Vals++;
	double sum = 0.0;
	for (j=0; j < NumEntries; j++) 
	  sum += RowValues[j] * yp[RowIndices[j]];
				
	if (UnitDiagonal) 
	  yp[i] = xp[i] - sum;
	else 
	  yp[i] = (xp[i] - sum)/RowValues[NumEntries];
				
      }
    }
  }
	
  // ***********  Transpose case *******************************
	
  else {

    if (xp!=yp) 
      for (i=0; i < NumMyCols_; i++) 
	yp[i] = xp[i]; // Initialize y for transpose solve
    
    if (Upper) {

      j0 = 1;
      if (NoDiagonal()) 
	j0--; // Include first term if no diagonal
    
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Vals++;
	if (!UnitDiagonal) 
	  yp[i] = yp[i]/RowValues[0];
     double ytmp = yp[i];
	for (j=j0; j < NumEntries; j++) 
	  yp[RowIndices[j]] -= RowValues[j] * ytmp;
      }
    }
    else {
			
      j0 = 1;
      if (NoDiagonal()) 
	j0--; // Include first term if no diagonal
    
      for (i=NumMyRows_-1; i >= 0; i--) {
	int      NumEntries = *NumEntriesPerRow-- - j0;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Vals--;
	if (!UnitDiagonal) 
	  yp[i] = yp[i]/RowValues[NumEntries];
     double ytmp = yp[i];
	for (j=0; j < NumEntries; j++) 
	  yp[RowIndices[j]] -= RowValues[j] * ytmp;
      }
    }
		
  }
  UpdateFlops(2*NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Solve1(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Solve(Upper,Trans,UnitDiag,X,Y)");
#endif

  //
  // This function find Y such that LY = X or UY = X or the transpose cases.
  //
  if((X.NumVectors() == 1) && (Y.NumVectors() == 1)) {
    double* xp = (double*) X[0];
    double* yp = (double*) Y[0];
    Epetra_Vector x(View, X.Map(), xp);
    Epetra_Vector y(View, Y.Map(), yp);
    EPETRA_CHK_ERR(Solve1(Upper, Trans, UnitDiagonal, x, y));
    return(0);
  }
  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  if((Upper) && (!UpperTriangular()))
    EPETRA_CHK_ERR(-2);
  if((!Upper) && (!LowerTriangular()))
    EPETRA_CHK_ERR(-3);
  if((!UnitDiagonal) && (NoDiagonal()))
    EPETRA_CHK_ERR(-4); // If matrix has no diagonal, we must use UnitDiagonal
  if((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_))
    EPETRA_CHK_ERR(-5); // Need each row to have a diagonal

  int i, j, j0, k;
  int* NumEntriesPerRow = Graph().NumIndicesPerRow();
  int** Indices = Graph().Indices();
  double** Vals = Values();
  double diag = 0.0;

  // If upper, point to last row
  if((Upper && !Trans) || (!Upper && Trans)) {
    NumEntriesPerRow += NumMyRows_-1;
    Indices += NumMyRows_-1;
    Vals += NumMyRows_-1;
  }

  double** Xp = (double**) X.Pointers();
  double** Yp = (double**) Y.Pointers();

  int NumVectors = X.NumVectors();

  if(!Trans) {   
    if(Upper) {   
      j0 = 1;
      if(NoDiagonal()) 
	j0--; // Include first term if no diagonal
      for(i = NumMyRows_ - 1; i >= 0; i--) {
	int     NumEntries = *NumEntriesPerRow--;
	int*    RowIndices = *Indices--;
	double* RowValues  = *Vals--;
	if(!UnitDiagonal) 
	  diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
	for(k = 0; k < NumVectors; k++) {
	  double sum = 0.0;
	  for(j = j0; j < NumEntries; j++) 
	    sum += RowValues[j] * Yp[k][RowIndices[j]];
					
	  if(UnitDiagonal) 
	    Yp[k][i] = Xp[k][i] - sum;
	  else 
	    Yp[k][i] = (Xp[k][i] - sum) * diag;
	}
      }
    }
    else {
      j0 = 1;
      if(NoDiagonal()) 
	j0--; // Include first term if no diagonal
      for(i = 0; i < NumMyRows_; i++) {
	int     NumEntries = *NumEntriesPerRow++ - j0;
	int*    RowIndices = *Indices++;
	double* RowValues  = *Vals++;
	if(!UnitDiagonal)
	  diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
	for(k = 0; k < NumVectors; k++) {
	  double sum = 0.0;
	  for(j = 0; j < NumEntries; j++) 
	    sum += RowValues[j] * Yp[k][RowIndices[j]];
					
	  if(UnitDiagonal) 
	    Yp[k][i] = Xp[k][i] - sum;
	  else 
	    Yp[k][i] = (Xp[k][i] - sum)*diag;
	}
      }
    }
  }
  // ***********  Transpose case *******************************

  else {
    for(k = 0; k < NumVectors; k++) 
      if(Yp[k] != Xp[k]) 
	for(i = 0; i < NumMyRows_; i++)
	  Yp[k][i] = Xp[k][i]; // Initialize y for transpose multiply
    
    if(Upper) {
      j0 = 1;
      if(NoDiagonal()) 
	j0--; // Include first term if no diagonal
      
      for(i = 0; i < NumMyRows_; i++) {
	int     NumEntries = *NumEntriesPerRow++;
	int*    RowIndices = *Indices++;
	double* RowValues  = *Vals++;
	if(!UnitDiagonal) 
	  diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
	for(k = 0; k < NumVectors; k++) {
	  if(!UnitDiagonal) 
	    Yp[k][i] = Yp[k][i]*diag;
       double ytmp = Yp[k][i];
	  for(j = j0; j < NumEntries; j++) 
	    Yp[k][RowIndices[j]] -= RowValues[j] * ytmp;
	}
      }
    }
    else {
      j0 = 1;
      if(NoDiagonal()) 
	j0--; // Include first term if no diagonal  
      for(i = NumMyRows_ - 1; i >= 0; i--) {
	int     NumEntries = *NumEntriesPerRow-- - j0;
	int*    RowIndices = *Indices--;
	double* RowValues  = *Vals--;
     if (!UnitDiagonal)
       diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
	for(k = 0; k < NumVectors; k++) {
	  if(!UnitDiagonal)  
	    Yp[k][i] = Yp[k][i]*diag;
       double ytmp = Yp[k][i];
	  for(j = 0; j < NumEntries; j++)
	    Yp[k][RowIndices[j]] -= RowValues[j] * ytmp;
        }
      }
    }
  }
  
  UpdateFlops(2 * NumVectors * NumGlobalNonzeros());
  return(0);
}

