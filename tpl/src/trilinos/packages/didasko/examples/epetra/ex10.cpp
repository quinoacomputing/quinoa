
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

// Define dense serial matrices
// This code should be run with one process

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_SerialDenseMatrix.h"

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // declare two dense matrix, whose dimensions are still not specified
  Epetra_SerialDenseMatrix A, B;

  // Total number of rows ans columns for dense matrix A
  int NumRowsA = 2, NumColsA = 2;

  // shape A
  A.Shape( NumRowsA, NumColsA );

  // set the element of A using the () operator.
  // Note that i is the row-index, and j the column-index
  for( int i=0 ; i<NumRowsA ; ++i )
    for( int j=0 ; j<NumColsA ; ++j )
      A(i,j) = i+100*j;

  // Epetra_SerialDenseMatrix overloads the << operator
  cout << A;

  // get matrix norms
  cout << "Inf norm of A = " << A.OneNorm() << endl;
  cout << "One norm of A = " << A.InfNorm() << endl;

  // now define an other matrix, B, for matrix multiplication
  int NumRowsB = 2, NumColsB=1;
  B.Shape(NumRowsB, NumColsB);

  // enter the values of B
  for( int i=0 ; i<NumRowsB ; ++i )
    for( int j=0 ; j<NumColsB ; ++j )
      B(i,j) = 11.0+i+100*j;

  cout << B;

  // define the matrix which will hold A * B
  Epetra_SerialDenseMatrix AtimesB;

  // same number of rows than A, same columns than B
  AtimesB.Shape(NumRowsA,NumColsB);

  // A * B
  AtimesB.Multiply('N','N',1.0, A, B, 0.0);
  cout << AtimesB;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

} /* main */

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
