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

// Show the usage of RowMap, ColumnMap, RangeMap and DomainMap
// This code must be run with two processes

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if( NumProc != 2 ) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return(EXIT_SUCCESS);
  }

  // I define two maps:
  // - MapA has 2 nodes on proc 0, 2 nodes on proc 1
  // - MapB has 4 nodes on proc 0, 0 nodes on proc 1

  int NumElementsA = 0, NumElementsB = 0;

  switch( MyPID ) {
    case 0:
      NumElementsA = 2;
      NumElementsB = 4;
      break;

    case 1:
      NumElementsA = 2;
      NumElementsB = 0;
      break;
  }

  Epetra_Map MapA(-1,NumElementsA,0,Comm);
  Epetra_Map MapB(-1,NumElementsB,0,Comm);

  // I create a diagonal matrix Q with RowMap,ColMap distibutions = MapA
  Epetra_CrsMatrix Q(Copy,MapA,MapA,1);

  int * MyGlobalElementsA = MapA.MyGlobalElements();

  // Now I put in Q the diagonal elements, using MapA
  for( int i=0 ; i<NumElementsA ; ++i ) {
    double one = 2.0;
    int indices = MyGlobalElementsA[i];
    Q.InsertGlobalValues(MyGlobalElementsA[i], 1, &one, &indices );
  }

  // Now, I would like to apply Q to a vector defined on
  // MapB, giving as a result a vector on MapB
  Q.FillComplete(MapB,MapB);

  // the instruction above causes the code to crash if
  // I build up the matrix as `Epetra_CrsMatrix Q(Copy,MapA,1)'
  // (without specifing the col map)

  // create few vectors on the two maps
  Epetra_Vector VecA(MapA);   Epetra_Vector VecA2(MapA);
  Epetra_Vector VecB(MapB);   Epetra_Vector VecB2(MapB);

  VecA.PutScalar(1.0);        VecA2.PutScalar(1.0);
  VecB.PutScalar(1.0);        VecB2.PutScalar(1.0);

  Q.Multiply(false,VecB,VecB2);

  cout << VecB2;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);

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

