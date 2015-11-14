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

// Solve a 2D Laplacian problem
// This example builds the matrix and solves it with AztecOO.

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_AZTECOO)

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
#include "AztecOO.h"

// external function
void  get_neighbours( const int i, const int nx, const int ny,
    int & left, int & right,
    int & lower, int & upper);

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // number of nodes along the x- and y-axis
  int nx = 5;
  int ny = 6;
  int NumGlobalElements = nx * ny;

  // create a linear map
  Epetra_Map Map(NumGlobalElements,0,Comm);

  // local number of rows
  int NumMyElements = Map.NumMyElements();
  // get update list
  int * MyGlobalElements = new int [NumMyElements];
  Map.MyGlobalElements( MyGlobalElements );

  // Create a Epetra_Matrix with 5 nonzero per rows

  Epetra_CrsMatrix A(Copy,Map,5);

  // Add  rows one-at-a-time
  // Need some vectors to help

  double Values[4];
  int Indices[4];
  int left, right, lower, upper;
  double diag = 4.0;

  for( int i=0 ; i<NumMyElements; ++i ) {
    int NumEntries=0;
    get_neighbours(  MyGlobalElements[i], nx, ny,
        left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    // put the off-diagonal entries
    A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
    // Put in the diagonal entry
    A.InsertGlobalValues(MyGlobalElements[i], 1, &diag, MyGlobalElements+i);
  }

  // Finish up
  A.FillComplete();

  // create x and b vectors
  Epetra_Vector x(Map);
  Epetra_Vector b(Map);
  b.PutScalar(1.0);

  // ==================== AZTECOO INTERFACE ======================

  // create linear problem
  Epetra_LinearProblem Problem(&A,&x,&b);
  // create AztecOO instance
  AztecOO Solver(Problem);

  Solver.SetAztecOption( AZ_precond, AZ_Jacobi );

  Solver.Iterate(1000,1E-9);

  // ==================== END OF AZTECOO INTERFACE ================

  if( Comm.MyPID() == 0 ) {
    cout << "Solver performed " << Solver.NumIters()
      << "iterations.\n";
    cout << "Norm of the true residual = " << Solver.TrueResidual() << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);

}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void  get_neighbours( const int i, const int nx, const int ny,
    int & left, int & right,
    int & lower, int & upper)
{

  int ix, iy;
  ix = i%nx;
  iy = (i - ix)/nx;

  if( ix == 0 )
    left = -1;
  else
    left = i-1;
  if( ix == nx-1 )
    right = -1;
  else
    right = i+1;
  if( iy == 0 )
    lower = -1;
  else
    lower = i-nx;
  if( iy == ny-1 )
    upper = -1;
  else
    upper = i+nx;

  return;

}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure Didasko with:\n"
      "--enable-epetra\n"
      "--enable-aztecoo");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

#endif

