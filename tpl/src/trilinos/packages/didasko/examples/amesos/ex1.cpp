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

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_AMESOS) && defined(HAVE_DIDASKO_TEUCHOS) && defined(HAVE_DIDASKO_TRIUTILS)

#include "Epetra_ConfigDefs.h"
#include "Amesos_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Amesos.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Trilinos_Util;

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // initialize an Gallery object
  CrsMatrixGallery Gallery("laplace_2d", Comm, false); // CJ TODO FIXME: change for Epetra64
  Gallery.Set("problem_size", 100); //must be a square number

  // get pointers to the linear problem, containing matrix, LHS and RHS.
  // if you need to access them, you can for example uncomment the following
  // code:
  // Epetra_CrsMatrix* Matrix = Gallery.GetMatrix();
  // Epetra_MultiVector* LHS = Gallery.GetStartingSolution();
  // Epetra_MultiVector* RHS = Gallery.GetRHS();
  //
  // NOTE: StartingSolution and RHS are pointers to Gallery's internally stored
  // vectors. Using StartingSolution and RHS, we can verify the residual
  // after the solution of the linear system. However, users may define as well
  // their own vectors for solution and RHS.

  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  // initialize Amesos solver:
  // `Solver' is the pointer to the Amesos solver
  // (note the use of the base class Amesos_BaseSolver)
  Amesos_BaseSolver* Solver;
  // Amesos_Factory is the function class used to create the solver.
  // This class contains no data.
  Amesos Amesos_Factory;

  // empty parameter list
  Teuchos::ParameterList List;

  // may also try: "Amesos_Umfpack", "Amesos_Lapack", ...
  string SolverType = "Amesos_Klu";

  Solver = Amesos_Factory.Create(SolverType, *Problem);
  // Amesos_Factory returns 0 is the selected solver is not
  // available
  assert (Solver);

  // start solving
  Solver->SymbolicFactorization();
  Solver->NumericFactorization();
  Solver->Solve();

  // verify that residual is really small
  double residual, diff;

  Gallery.ComputeResidual(&residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);

  if( Comm.MyPID() == 0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }

  // delete Solver
  delete Solver;

  if (residual > 1e-5)
    exit(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  exit(EXIT_SUCCESS);
}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure Didasko with:\n"
      "--enable-epetra\n"
      "--enable-teuchos\n"
      "--enable-triutils\n"
      "--enable-amesos");

  return 0;
}

#endif

