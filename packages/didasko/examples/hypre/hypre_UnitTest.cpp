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

#include "Ifpack.h"
#include "AztecOO.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "EpetraExt_HypreIJMatrix.h"
#include "Epetra_MultiVector.h"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "hypre_Helpers.hpp"
#include "Teuchos_Array.hpp"
#include <string>
#include <stdio.h>
#include <map>

using Teuchos::RCP;
using Teuchos::rcp;

const double tol = 1E-7;
/*
TEUCHOS_UNIT_TEST( Ifpack_Hypre, AztecOO ){
  Epetra_CrsMatrix* Crs_Matrix = newCrsMatrix(3);
  Ifpack_Hypre preconditioner(Crs_Matrix);
  preconditioner.Initialize();
  int NumProc = Crs_Matrix->Comm().NumProc();
  HYPRE_IJMatrix hypre_mat = preconditioner.HypreMatrix();
  EpetraExt_HypreIJMatrix Hyp_Matrix(hypre_mat);
  TEST_EQUALITY(EquivalentMatrices(Hyp_Matrix, *Crs_Matrix, tol), true);

  int numVec = 1;
  Epetra_MultiVector X(Hyp_Matrix.OperatorDomainMap(), numVec);
  Epetra_MultiVector KnownX(Hyp_Matrix.OperatorDomainMap(), numVec);
  KnownX.Random();
  Epetra_MultiVector B(Hyp_Matrix.OperatorRangeMap(), numVec);
  Hyp_Matrix.Multiply(false, KnownX, B);

  Epetra_LinearProblem Problem(&Hyp_Matrix, &X, &B);
  AztecOO solver(Problem);
  Hyp_Matrix.SetParameter(Solver, &HYPRE_ParCSRPCGSetMaxIter, 1000);
  Hyp_Matrix.SetParameter(Solver, &HYPRE_ParCSRPCGSetTol, 1e-9);
  solver.SetPrecOperator(&Hyp_Matrix);
  solver.Iterate(1000, tol);
  TEST_EQUALITY(EquivalentVectors(X, KnownX, tol*10*pow(10.0,NumProc)), true);
  
  X.PutScalar(0.0);
  Epetra_LinearProblem Problem2(Crs_Matrix, &X, &B);
  AztecOO solver2(Problem2);
  preconditioner.SetParameter(Solver, PCG); // Use a PCG Solver
  preconditioner.SetParameter(Preconditioner, ParaSails); // Use a Euclid Preconditioner (but not really used)
  preconditioner.SetParameter(Solver, &HYPRE_ParCSRPCGSetMaxIter, 1000); // Set maximum iterations
  preconditioner.SetParameter(Solver, &HYPRE_ParCSRPCGSetTol, 1e-9); // Set a tolerance
  preconditioner.SetParameter(Solver); // Solve the problem
  preconditioner.SetParameter(true); // Don't use the preconditioner
  solver2.SetPrecOperator(&preconditioner);
  preconditioner.Compute();
  preconditioner.Condest(Ifpack_Cheap, 1550, 1e-9, Crs_Matrix);
  cout << endl << preconditioner << endl;
  solver2.Iterate(1000, tol);
  TEST_EQUALITY(EquivalentVectors(X, KnownX, tol*10*pow(10.0,NumProc)), true);
}*/

TEUCHOS_UNIT_TEST( Ifpack_Hypre, Ifpack ){

  Epetra_CrsMatrix* Crs_Matrix = newCrsMatrix(4);
  Ifpack Factory;
  RCP<Ifpack_Preconditioner> preconditioner= rcp(Factory.Create("Hypre", Crs_Matrix));
  int NumProc = Crs_Matrix->Comm().NumProc();
  int MyPID = Crs_Matrix->Comm().MyPID();

  int numVec = 2;
  Epetra_MultiVector X(preconditioner->OperatorRangeMap(), numVec);
  Epetra_MultiVector KnownX(preconditioner->OperatorRangeMap(), numVec);
  TEST_EQUALITY(KnownX.Random(), 0);
  Epetra_MultiVector B(preconditioner->OperatorDomainMap(), numVec);
  TEST_EQUALITY(preconditioner->Initialize(), 0);
  TEST_EQUALITY(preconditioner->Apply(KnownX, B), 0);

  Teuchos::ParameterList list("New List");
  RCP<FunctionParameter> functs[5];
  functs[0] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetMaxIter, 1000));
  functs[1] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTol, 1e-9));
  functs[2] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetLogging, 1));
  functs[3] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetPrintLevel, 0));
  functs[4] = rcp(new FunctionParameter(Preconditioner, &HYPRE_ParaSailsSetLogging, 0));
  list.set("NumFunctions", 5);
  list.set<RCP<FunctionParameter>*>("Functions", functs);
  list.set("SolveOrPrecondition", Solver);
  list.set("Solver", PCG);
  list.set("Preconditioner", ParaSails);
  list.set("SetPreconditioner", true);
  TEST_EQUALITY(preconditioner->SetParameters(list), 0);
  //(dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Solver, PCG); // Use a PCG Solver
  //(dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Preconditioner, ParaSails); // Use a ParaSails Preconditioner
  //(dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Solver, &HYPRE_ParCSRPCGSetMaxIter, 1000); // Set maximum iterations
  //(dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Solver, &HYPRE_ParCSRPCGSetTol, 1e-9); // Set a tolerance
  //(dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Solver); // Solve the problem
  //(dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(true); // Use the preconditioner
  TEST_EQUALITY(preconditioner->Compute(),0);
  preconditioner->Condest(Ifpack_Cheap, 1000, 1e-9, Crs_Matrix);
  TEST_EQUALITY(preconditioner->ApplyInverse(B, X),0);
  TEST_EQUALITY(EquivalentVectors(X, KnownX, tol*10*pow(10.0,NumProc)), true);
  cout << *preconditioner;
  int numIters;
  double residual;
  (dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Solver, &HYPRE_ParCSRPCGGetNumIterations, &numIters); 
  (dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->SetParameter(Solver, &HYPRE_ParCSRPCGGetFinalRelativeResidualNorm, &residual); 
  (dynamic_cast<Ifpack_Hypre*> (preconditioner.get()))->CallFunctions();
  if(MyPID == 0) printf("It took %d iterations, and achieved %e residual.\n", numIters, residual);
  delete Crs_Matrix;
}

TEUCHOS_UNIT_TEST( Ifpack_Hypre, EpetraExt ){
  
  Epetra_CrsMatrix* Crs_Matrix = newCrsMatrix(3);
  Ifpack_Hypre preconditioner(Crs_Matrix);
  int NumProc = Crs_Matrix->Comm().NumProc();
  int MyPID = Crs_Matrix->Comm().MyPID();
  HYPRE_IJMatrix hypre_mat = preconditioner.HypreMatrix();
  EpetraExt_HypreIJMatrix Hyp_Matrix(hypre_mat);
  EquivalentMatrices(Hyp_Matrix, *Crs_Matrix, tol);

  int numVec = 2;
  Epetra_MultiVector X(Hyp_Matrix.RowMatrixRowMap(), numVec);
  Epetra_MultiVector KnownX(Hyp_Matrix.RowMatrixRowMap(), numVec);
  KnownX.Random();
  Epetra_MultiVector B(Hyp_Matrix.RowMatrixRowMap(), numVec);
  Hyp_Matrix.Multiply(false, KnownX, B);
  Hyp_Matrix.SetParameter(Solver, HYPRE_ParCSRPCGSetMaxIter, 1000);
  Hyp_Matrix.SetParameter(Solver, HYPRE_ParCSRPCGSetTol, 1e-9);
  Hyp_Matrix.SetParameter(false);
  Hyp_Matrix.ApplyInverse(B, X);
  TEST_EQUALITY(EquivalentVectors(X, KnownX, tol*10*pow(10.0,NumProc)), true);

  int numIters;
  double residual;
  Hyp_Matrix.SetParameter(Solver, &HYPRE_ParCSRPCGGetFinalRelativeResidualNorm, &residual);
  Hyp_Matrix.SetParameter(Solver, &HYPRE_ParCSRPCGGetNumIterations, &numIters);
  if(MyPID == 0) printf("It took %d iterations, and achieved %e residual.\n", numIters, residual);
  /* This proves the problem is in the maps of the vectors
  X.PutScalar(0.0);
  Teuchos::ParameterList list("New List");
  RCP<FunctionParameter> functs[2];
  functs[0] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetMaxIter, 1000));
  functs[1] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTol, 1e-9)); 
  list.set("NumFunctions", 2);
  list.set<RCP<FunctionParameter>*>("Functions", functs);
  list.set("SolveOrPrecondition", Solver);
  list.set("SetPreconditioner", true);
  preconditioner.SetParameters(list);
  preconditioner.Initialize();
  preconditioner.Compute();
  preconditioner.ApplyInverse(B, X);
  TEST_EQUALITY(EquivalentVectors(X, KnownX, tol*10*pow(10.0,NumProc)), true);
  */
  delete Crs_Matrix;
}
