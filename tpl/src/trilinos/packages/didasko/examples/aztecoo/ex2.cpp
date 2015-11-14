
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

//  Solve a linear system with Aztec, using Aztec as a preconditioner
// (recursive way)
//
// NOTE: This example first builds the matrix, then solves it with AztecOO
//
// NOTE2: this example implemenets minor modifications to one of the
// examples included in the AztecOO package. Please give a look  to file
// ${TRILINOS_HOME}/packages/aztecoo/examples/AztecOO_RecursiveCall/cxx_main.cpp
// for more details.

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_TRIUTILS)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "AztecOO_Operator.h"

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  /* this example creates a tridiagonal matrix of type
   *
   *     |  2  -1            |
   *     | -1   2   -1       |
   * A = |      ...  ... ... |
   *     |            -1  2  |
   */

  // set global dimension to 5, could be any number
  int NumGlobalElements = 5;

  // create a linear map
  Epetra_Map Map(NumGlobalElements,0,Comm);
  // local number of rows
  int NumMyElements = Map.NumMyElements();
  // get update list
  int * MyGlobalElements = new int [NumMyElements];
  Map.MyGlobalElements( MyGlobalElements );

  // Create a Epetra_Matrix
  // create a CSR matrix

  Epetra_CrsMatrix A(Copy,Map,3);

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1

  double *Values = new double[2];
  Values[0] = -1.0; Values[1] = -1.0;
  int *Indices = new int[2];
  double two = 2.0;
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
    A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
    // Put in the diagonal entry
    A.InsertGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i);
  }

  // Finish up
  A.FillComplete();

  // E N D   O F   M A T R I X   C O N S T R U C T I O N

  Epetra_Vector x(Map);
  Epetra_Vector b(Map);

  x.Random();

  // this is the linear problem to be solved: set the linear operator,
  // the solution and the right-hand side
  Epetra_LinearProblem A_Problem(&A, &x, &b);

  // and this is the AztecOO solver
  AztecOO A_Solver(A_Problem);

  // --- Here we define the precondioner ---
  // create P as a copy of A, in principle can be different
  Epetra_CrsMatrix P(A);

  // Here we create the linear problem which will be used as
  // preconditioner. This requires sereval steps.
  // (Note that all the P_ prefix indentify preconditioner'
  // objects)

  // 1. we create the linear system solve at each prec step
  Epetra_LinearProblem P_Problem;
  // and we assign the linear operator (in this case, the
  // matrix A itself)
  P_Problem.SetOperator(&P);

  // as we wish to use AztecOO to solve the prec step
  // (in a recursive way), we have to define an AztecOO
  // object.
  AztecOO P_Solver(P_Problem);

  // now, we customize certain parameters
  P_Solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  P_Solver.SetAztecOption(AZ_output, AZ_none);
  P_Solver.SetAztecOption(AZ_solver, AZ_cg);
  // The last step is to create an AztecOO_Operator, so that
  // we can set the Aztec's preconditioner with.
  AztecOO_Operator P_Operator(&P_Solver, 10);

  // Here we set the user's defined preconditioners
  A_Solver.SetPrecOperator(&P_Operator);

  // --- up to here ---

  // Finally, we solve the linear system:

  int Niters=100;
  A_Solver.SetAztecOption(AZ_kspace, Niters);
  A_Solver.SetAztecOption(AZ_solver, AZ_GMRESR);

  A_Solver.Iterate(Niters, 1.0E-12);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return 0;

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
      "--enable-triutils\n"
      "--enable-aztecoo\n");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

#endif
