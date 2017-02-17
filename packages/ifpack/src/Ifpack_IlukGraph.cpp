/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_IlukGraph.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_Import.h"

#include <Teuchos_ParameterList.hpp>
#include <ifp_parameters.h>

//==============================================================================
Ifpack_IlukGraph::Ifpack_IlukGraph(const Epetra_CrsGraph & Graph_in, int LevelFill_in, int LevelOverlap_in)
  : Graph_(Graph_in),
    DomainMap_(Graph_in.DomainMap()),
    RangeMap_(Graph_in.RangeMap()),
    Comm_(Graph_in.Comm()),
    LevelFill_(LevelFill_in),
    LevelOverlap_(LevelOverlap_in),
    IndexBase_(Graph_in.IndexBase64()),
    NumGlobalRows_(Graph_in.NumGlobalRows64()),
    NumGlobalCols_(Graph_in.NumGlobalCols64()),
    NumGlobalBlockRows_(Graph_in.NumGlobalBlockRows64()),
    NumGlobalBlockCols_(Graph_in.NumGlobalBlockCols64()),
    NumGlobalBlockDiagonals_(0),
    NumGlobalNonzeros_(0),
    NumGlobalEntries_(0),
    NumMyBlockRows_(Graph_in.NumMyBlockRows()),
    NumMyBlockCols_(Graph_in.NumMyBlockCols()),
    NumMyRows_(Graph_in.NumMyRows()),
    NumMyCols_(Graph_in.NumMyCols()),
    NumMyBlockDiagonals_(0),
    NumMyNonzeros_(0),
    NumMyEntries_(0)
{
}

//==============================================================================
Ifpack_IlukGraph::Ifpack_IlukGraph(const Ifpack_IlukGraph & Graph_in)
  : Graph_(Graph_in.Graph_),
    DomainMap_(Graph_in.DomainMap()),
    RangeMap_(Graph_in.RangeMap()),
    Comm_(Graph_in.Comm()),
    OverlapGraph_(Graph_in.OverlapGraph_),
    OverlapRowMap_(Graph_in.OverlapRowMap_),
    OverlapImporter_(Graph_in.OverlapImporter_),
    LevelFill_(Graph_in.LevelFill_),
    LevelOverlap_(Graph_in.LevelOverlap_),
    IndexBase_(Graph_in.IndexBase_),
    NumGlobalRows_(Graph_in.NumGlobalRows_),
    NumGlobalCols_(Graph_in.NumGlobalCols_),
    NumGlobalBlockRows_(Graph_in.NumGlobalBlockRows_),
    NumGlobalBlockCols_(Graph_in.NumGlobalBlockCols_),
    NumGlobalBlockDiagonals_(Graph_in.NumGlobalBlockDiagonals_),
    NumGlobalNonzeros_(Graph_in.NumGlobalNonzeros_),
    NumGlobalEntries_(Graph_in.NumGlobalEntries_),
    NumMyBlockRows_(Graph_in.NumMyBlockRows_),
    NumMyBlockCols_(Graph_in.NumMyBlockCols_),
    NumMyRows_(Graph_in.NumMyRows_),
    NumMyCols_(Graph_in.NumMyCols_),
    NumMyBlockDiagonals_(Graph_in.NumMyBlockDiagonals_),
    NumMyNonzeros_(Graph_in.NumMyNonzeros_),
    NumMyEntries_(Graph_in.NumMyEntries_)
{
  Epetra_CrsGraph & L_Graph_In = Graph_in.L_Graph();
  Epetra_CrsGraph & U_Graph_In = Graph_in.U_Graph();
  L_Graph_ = Teuchos::rcp( new Epetra_CrsGraph(L_Graph_In) );
  U_Graph_ = Teuchos::rcp( new Epetra_CrsGraph(U_Graph_In) );
}

//==============================================================================
Ifpack_IlukGraph::~Ifpack_IlukGraph()
{
}

//==============================================================================
int Ifpack_IlukGraph::SetParameters(const Teuchos::ParameterList& parameterlist,
                                    bool cerr_warning_if_unused)
{
  Ifpack::param_struct params;
  params.int_params[Ifpack::level_fill-FIRST_INT_PARAM] = LevelFill_;
  params.int_params[Ifpack::level_overlap-FIRST_INT_PARAM] = LevelOverlap_;

  Ifpack::set_parameters(parameterlist, params, cerr_warning_if_unused);

  LevelFill_ = params.int_params[Ifpack::level_fill-FIRST_INT_PARAM];
  LevelOverlap_ = params.int_params[Ifpack::level_overlap-FIRST_INT_PARAM];
  return(0);
}

//==============================================================================
int Ifpack_IlukGraph::ConstructOverlapGraph() {

  OverlapGraph_ = Teuchos::rcp( (Epetra_CrsGraph *) &Graph_, false );
  OverlapRowMap_ = Teuchos::rcp( (Epetra_BlockMap *) &Graph_.RowMap(), false );

  if (LevelOverlap_==0 || !Graph_.DomainMap().DistributedGlobal()) return(0); // Nothing to do

  Teuchos::RefCountPtr<Epetra_CrsGraph> OldGraph;
  Teuchos::RefCountPtr<Epetra_BlockMap> OldRowMap;
  Epetra_BlockMap * DomainMap_tmp = (Epetra_BlockMap *) &Graph_.DomainMap();
  Epetra_BlockMap * RangeMap_tmp = (Epetra_BlockMap *) &Graph_.RangeMap();
  for (int level=1; level <= LevelOverlap_; level++) {
    OldGraph = OverlapGraph_;
    OldRowMap = OverlapRowMap_;

    OverlapImporter_ = Teuchos::rcp( (Epetra_Import *) OldGraph->Importer(), false );
    OverlapRowMap_ = Teuchos::rcp( new Epetra_BlockMap(OverlapImporter_->TargetMap()) );


    if (level<LevelOverlap_)
      OverlapGraph_ = Teuchos::rcp( new Epetra_CrsGraph(Copy, *OverlapRowMap_, 0) );
    else
      // On last iteration, we want to filter out all columns except those that correspond
      // to rows in the graph.  This assures that our matrix is square
      OverlapGraph_ = Teuchos::rcp( new Epetra_CrsGraph(Copy, *OverlapRowMap_, *OverlapRowMap_, 0) );

    EPETRA_CHK_ERR(OverlapGraph_->Import( Graph_, *OverlapImporter_, Insert));
    if (level<LevelOverlap_) {
      EPETRA_CHK_ERR(OverlapGraph_->FillComplete(*DomainMap_tmp, *RangeMap_tmp));
    }
    else {
      // Copy last OverlapImporter because we will use it later
      OverlapImporter_ = Teuchos::rcp( new Epetra_Import(*OverlapRowMap_, *DomainMap_tmp) );
      EPETRA_CHK_ERR(OverlapGraph_->FillComplete(*DomainMap_tmp, *RangeMap_tmp));
    }
  }

    NumMyBlockRows_ = OverlapGraph_->NumMyBlockRows();
    NumMyBlockCols_ = OverlapGraph_->NumMyBlockCols();
    NumMyRows_ = OverlapGraph_->NumMyRows();
    NumMyCols_ = OverlapGraph_->NumMyCols();

  return(0);
}

//==============================================================================
int Ifpack_IlukGraph::ConstructFilledGraph() {
  using std::cout;
  using std::endl;

  int ierr = 0;
  int i, j;
  int * In=0;
  int NumIn, NumL, NumU;
  bool DiagFound;


  EPETRA_CHK_ERR(ConstructOverlapGraph());

  L_Graph_ = Teuchos::rcp( new Epetra_CrsGraph(Copy, OverlapGraph_->RowMap(), OverlapGraph_->RowMap(),  0) );
  U_Graph_ = Teuchos::rcp( new Epetra_CrsGraph(Copy, OverlapGraph_->RowMap(), OverlapGraph_->RowMap(),  0));


  // Get Maximun Row length
  int MaxNumIndices = OverlapGraph_->MaxNumIndices();

  std::vector<int> L(MaxNumIndices);
  std::vector<int> U(MaxNumIndices);

  // First we copy the user's graph into L and U, regardless of fill level

  for (i=0; i< NumMyBlockRows_; i++) {


    OverlapGraph_->ExtractMyRowView(i, NumIn, In); // Get Indices


    // Split into L and U (we don't assume that indices are ordered).

    NumL = 0;
    NumU = 0;
    DiagFound = false;

    for (j=0; j< NumIn; j++) {
      int k = In[j];

      if (k<NumMyBlockRows_) { // Ignore column elements that are not in the square matrix

        if (k==i) DiagFound = true;

        else if (k < i) {
          L[NumL] = k;
          NumL++;
        }
        else {
          U[NumU] = k;
          NumU++;
        }
      }
    }

    // Check in things for this row of L and U

    if (DiagFound) NumMyBlockDiagonals_++;
    if (NumL) L_Graph_->InsertMyIndices(i, NumL, &L[0]);
    if (NumU) U_Graph_->InsertMyIndices(i, NumU, &U[0]);

  }

  if (LevelFill_ > 0) {

    // Complete Fill steps
    Epetra_BlockMap * L_DomainMap = (Epetra_BlockMap *) &OverlapGraph_->RowMap();
    Epetra_BlockMap * L_RangeMap = (Epetra_BlockMap *) &Graph_.RangeMap();
    Epetra_BlockMap * U_DomainMap = (Epetra_BlockMap *) &Graph_.DomainMap();
    Epetra_BlockMap * U_RangeMap = (Epetra_BlockMap *) &OverlapGraph_->RowMap();
    EPETRA_CHK_ERR(L_Graph_->FillComplete(*L_DomainMap, *L_RangeMap));
    EPETRA_CHK_ERR(U_Graph_->FillComplete(*U_DomainMap, *U_RangeMap));

    // At this point L_Graph and U_Graph are filled with the pattern of input graph,
    // sorted and have redundant indices (if any) removed.  Indices are zero based.
    // LevelFill is greater than zero, so continue...

    int MaxRC = NumMyBlockRows_;
    std::vector<std::vector<int> > Levels(MaxRC);
    std::vector<int> LinkList(MaxRC);
    std::vector<int> CurrentLevel(MaxRC);
    std::vector<int> CurrentRow(MaxRC);
    std::vector<int> LevelsRowU(MaxRC);

    for (i=0; i<NumMyBlockRows_; i++)
    {
      int First, Next;

      // copy column indices of row into workspace and sort them

      int LenL = L_Graph_->NumMyIndices(i);
      int LenU = U_Graph_->NumMyIndices(i);
      int Len = LenL + LenU + 1;

      EPETRA_CHK_ERR(L_Graph_->ExtractMyRowCopy(i, LenL, LenL, &CurrentRow[0]));      // Get L Indices
      CurrentRow[LenL] = i;                                     // Put in Diagonal
      //EPETRA_CHK_ERR(U_Graph_->ExtractMyRowCopy(i, LenU, LenU, CurrentRow+LenL+1)); // Get U Indices
      int ierr1 = 0;
      if (LenU) {
        // Get U Indices
        ierr1 = U_Graph_->ExtractMyRowCopy(i, LenU, LenU, &CurrentRow[LenL+1]);
      }
      if (ierr1!=0) {
        cout << "ierr1 = "<< ierr1 << endl;
        cout << "i = " << i << endl;
        cout << "NumMyBlockRows_ = " << U_Graph_->NumMyBlockRows() << endl;
      }

      // Construct linked list for current row

      for (j=0; j<Len-1; j++) {
        LinkList[CurrentRow[j]] = CurrentRow[j+1];
        CurrentLevel[CurrentRow[j]] = 0;
      }

      LinkList[CurrentRow[Len-1]] = NumMyBlockRows_;
      CurrentLevel[CurrentRow[Len-1]] = 0;

      // Merge List with rows in U

      First = CurrentRow[0];
      Next = First;
      while (Next < i)
        {
          int PrevInList = Next;
          int NextInList = LinkList[Next];
          int RowU = Next;
          int LengthRowU;
          int * IndicesU;
          // Get Indices for this row of U
          EPETRA_CHK_ERR(U_Graph_->ExtractMyRowView(RowU, LengthRowU, IndicesU));

          int ii;

          // Scan RowU

          for (ii=0; ii<LengthRowU; /*nop*/)
            {
              int CurInList = IndicesU[ii];
              if (CurInList < NextInList)
                {
                  // new fill-in
                  int NewLevel = CurrentLevel[RowU] + Levels[RowU][ii+1] + 1;
                  if (NewLevel <= LevelFill_)
                    {
                      LinkList[PrevInList]  = CurInList;
                      LinkList[CurInList] = NextInList;
                      PrevInList = CurInList;
                      CurrentLevel[CurInList] = NewLevel;
                    }
                  ii++;
                }
              else if (CurInList == NextInList)
                {
                  PrevInList = NextInList;
                  NextInList = LinkList[PrevInList];
                  int NewLevel = CurrentLevel[RowU] + Levels[RowU][ii+1] + 1;
                  CurrentLevel[CurInList] = EPETRA_MIN(CurrentLevel[CurInList], NewLevel);
                  ii++;
                }
              else // (CurInList > NextInList)
                {
                  PrevInList = NextInList;
                  NextInList = LinkList[PrevInList];
                }
            }
          Next = LinkList[Next];
        }

      // Put pattern into L and U

      LenL = 0;

      Next = First;

      // Lower

      while (Next < i) {
        CurrentRow[LenL++] = Next;
        Next = LinkList[Next];
      }

      EPETRA_CHK_ERR(L_Graph_->RemoveMyIndices(i)); // Delete current set of Indices
      int ierr11 = L_Graph_->InsertMyIndices(i, LenL, &CurrentRow[0]);
      if (ierr11 < 0) EPETRA_CHK_ERR(ierr1);

      // Diagonal

      if (Next != i) return(-2); // Fatal:  U has zero diagonal.
      else {
        LevelsRowU[0] = CurrentLevel[Next];
        Next = LinkList[Next];
      }

      // Upper

      LenU = 0;

      while (Next < NumMyBlockRows_) // Should be "Next < NumMyBlockRows_"?
        {
          LevelsRowU[LenU+1] = CurrentLevel[Next];
          CurrentRow[LenU++] = Next;
          Next = LinkList[Next];
        }

      EPETRA_CHK_ERR(U_Graph_->RemoveMyIndices(i)); // Delete current set of Indices
      int ierr2 = U_Graph_->InsertMyIndices(i, LenU, &CurrentRow[0]);
      if (ierr2<0) EPETRA_CHK_ERR(ierr2);

      // Allocate and fill Level info for this row
      Levels[i] = std::vector<int>(LenU+1);
      for (int jj=0; jj<LenU+1; jj++) Levels[i][jj] = LevelsRowU[jj];

    }
  }

  // Complete Fill steps
  Epetra_BlockMap L_DomainMap = (Epetra_BlockMap) OverlapGraph_->RowMap();
  Epetra_BlockMap L_RangeMap = (Epetra_BlockMap) Graph_.RangeMap();
  Epetra_BlockMap U_DomainMap = (Epetra_BlockMap) Graph_.DomainMap();
  Epetra_BlockMap U_RangeMap = (Epetra_BlockMap) OverlapGraph_->RowMap();
  EPETRA_CHK_ERR(L_Graph_->FillComplete(L_DomainMap, L_RangeMap));
  EPETRA_CHK_ERR(U_Graph_->FillComplete(U_DomainMap, U_RangeMap));

  // Optimize graph storage

  EPETRA_CHK_ERR(L_Graph_->OptimizeStorage());
  EPETRA_CHK_ERR(U_Graph_->OptimizeStorage());

  // Compute global quantities

  NumGlobalBlockDiagonals_ = 0;
  long long NumMyBlockDiagonals_LL = NumMyBlockDiagonals_;
  EPETRA_CHK_ERR(L_Graph_->Comm().SumAll(&NumMyBlockDiagonals_LL, &NumGlobalBlockDiagonals_, 1));

  NumGlobalNonzeros_ = L_Graph_->NumGlobalNonzeros64()+U_Graph_->NumGlobalNonzeros64();
  NumMyNonzeros_ = L_Graph_->NumMyNonzeros()+U_Graph_->NumMyNonzeros();
  NumGlobalEntries_ = L_Graph_->NumGlobalEntries64()+U_Graph_->NumGlobalEntries64();
  NumMyEntries_ = L_Graph_->NumMyEntries()+U_Graph_->NumMyEntries();
  return(ierr);
}
//==========================================================================

// Non-member functions

std::ostream& operator << (std::ostream& os, const Ifpack_IlukGraph& A)
{
  using std::endl;

/*  Epetra_fmtflags olda = os.setf(ios::right,ios::adjustfield);
  Epetra_fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
  int oldp = os.precision(12); */
  int LevelFill = A.LevelFill();
  Epetra_CrsGraph & L = (Epetra_CrsGraph &) A.L_Graph();
  Epetra_CrsGraph & U = (Epetra_CrsGraph &) A.U_Graph();
  os.width(14);
  os <<  "     Level of Fill = "; os << LevelFill;
  os << endl;

  os.width(14);
  os <<  "     Graph of L = ";
  os << endl;
  os << L; // Let Epetra_CrsGraph handle the rest.

  os.width(14);
  os <<  "     Graph of U = ";
  os << endl;
  os << U; // Let Epetra_CrsGraph handle the rest.

  // Reset os flags

/*  os.setf(olda,ios::adjustfield);
  os.setf(oldf,ios::floatfield);
  os.precision(oldp); */

  return os;
}
