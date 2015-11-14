// @HEADER
// ***********************************************************************
//
//                      Didasko Tutorial Package
//                 Copyright (2005) Sandia Corporation
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
// Questions about Didasko? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ***********************************************************************
// @HEADER

// Output a CRS matrix in MATLAB format

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

// ============================================================
// define a class, derived from Epetra_CrsMatrix, which
// initializes the matrix entires. User has to provide
// a valid Epetra_Map in the constructor, plus the diagonal
// value, and the sub- and super-diagonal values.
// ============================================================

class TridiagonalCrsMatrix : public Epetra_CrsMatrix {

  public:
    TridiagonalCrsMatrix(const Epetra_Map & Map,
        double a,
        double diag, double c) :
      Epetra_CrsMatrix(Copy,Map,3)
  {
    // global number of rows
    int NumGlobalElements = Map.NumGlobalElements();
    // local number of rows
    int NumMyElements = Map.NumMyElements();
    // get update list
    int * MyGlobalElements = new int [NumMyElements];
    Map.MyGlobalElements( MyGlobalElements );

    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1

    double *Values = new double[2];
    Values[0] = a; Values[1] = c;
    int *Indices = new int[2];
    int NumEntries;

    for( int i=0 ; i<NumMyElements; ++i ) {
      if (MyGlobalElements[i]==0) {
        Indices[0] = 1;
        NumEntries = 1;
      } else if (MyGlobalElements[i] == NumGlobalElements-1) {
        Indices[0] = NumGlobalElements-2;
        NumEntries = 1;
      } else {
        Indices[0] = MyGlobalElements[i]-1;
        Indices[1] = MyGlobalElements[i]+1;
        NumEntries = 2;
      }
      InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
      // Put in the diagonal entry
      InsertGlobalValues(MyGlobalElements[i], 1, &diag, MyGlobalElements+i);
    }

    FillComplete();
    delete[] MyGlobalElements;
    delete[] Values;
    delete[] Indices;
  }

};

/* ======== ================ *
 * function CrsMatrix2MATLAB *
 * ======== ================ *
 *
 * Print out a CrsMatrix in a MATLAB format. Each processor prints out
 * its part, starting from proc 0 to proc NumProc-1. The first line of
 * each processor's output states the number of local rows and of
 * local nonzero elements. Output is finished by "End of Matrix Output".
 *
 *
 * Return code:        true if matrix has been printed out
 * -----------         false otherwise
 *
 * Parameters:
 * ----------
 *
 * - Epetra_CrsMatrix  reference to the ditributed CrsMatrix to
 *                     print out
 */

bool CrsMatrix2MATLAB( const Epetra_CrsMatrix & A )

{

  int MyPID = A.Comm().MyPID();
  int NumProc = A.Comm().NumProc();

  // work only on transformed matrices;
  if( A.IndicesAreLocal() == false ) {
    if( MyPID == 0 ) {
      cerr << "ERROR in "<< __FILE__ << ", line " << __LINE__ << endl;
      cerr << "Function CrsMatrix2MATLAB accepts\n";
      cerr << "transformed matrices ONLY. Please call A.FillComplete()\n";
      cerr << "on your matrix A to that purpose.\n";
      cerr << "Now returning...\n";
    }
    return false;
  }

  int NumMyRows = A.NumMyRows(); // number of rows on this process
  int NumNzRow;   // number of nonzero elements for each row
  int NumEntries; // number of extracted elements for each row
  int NumGlobalRows; // global dimensio of the problem
  int GlobalRow;  // row in global ordering
  int NumGlobalNonzeros; // global number of nonzero elements

  NumGlobalRows = A.NumGlobalRows();
  NumGlobalNonzeros = A.NumGlobalNonzeros();

  // print out on cout if no filename is provided

  int IndexBase = A.IndexBase(); // MATLAB start from 0
  if( IndexBase == 0 ) IndexBase = 1;

  // write on file the dimension of the matrix

  if( MyPID==0 ) {
    cout << "A = spalloc(";
    cout << NumGlobalRows << ',' << NumGlobalRows;
    cout << ',' << NumGlobalNonzeros << ");\n";
  }

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {

    if( MyPID == Proc ) {

      cout << "% On proc " << Proc << ": ";
      cout << NumMyRows << " rows and ";
      cout << A.NumMyNonzeros() << " nonzeros\n";

      // cycle over all local rows to find out nonzero elements
      for( int MyRow=0 ; MyRow<NumMyRows ; ++MyRow ) {

        GlobalRow = A.GRID(MyRow);

        NumNzRow = A.NumMyEntries(MyRow);
        double *Values = new double[NumNzRow];
        int *Indices = new int[NumNzRow];

        A.ExtractMyRowCopy(MyRow, NumNzRow,
            NumEntries, Values, Indices);
        // print out the elements with MATLAB syntax
        for( int j=0 ; j<NumEntries ; ++j ) {
          cout << "A(" << GlobalRow  + IndexBase
            << "," << A.GCID(Indices[j]) + IndexBase
            << ") = " << Values[j] << ";\n";
        }

        delete[] Values;
        delete[] Indices;
      }

    }
    A.Comm().Barrier();
    if( MyPID == 0 ) {
      cout << " %End of Matrix Output\n";
    }
  }

  return true;

}

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // set global dimension to 5, could be any number
  int NumGlobalElements = 5;

  // define a linear map
  Epetra_Map Map(NumGlobalElements,0,Comm);

  // create the matrix
  TridiagonalCrsMatrix A( Map, -1.0, 2.0, -1.0);

  // output informationto stdout
  CrsMatrix2MATLAB( A );

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure Didasko with:\n"
      "--enable-epetra");

  return 0;
}
#endif

