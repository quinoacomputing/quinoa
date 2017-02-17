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

#include "Ifpack_ConfigDefs.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Ifpack_AdditiveSchwarz.h"
#include "AztecOO.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_IC.h"
#include "Ifpack_ILU.h"
#include "Ifpack_Amesos.h"

bool verbose = false;

bool CompareWithAztecOO(Epetra_LinearProblem& Problem, const std::string what,
                       int Overlap, int ival)
{
  using std::cout;
  using std::endl;

  AztecOO AztecOOSolver(Problem);
  AztecOOSolver.SetAztecOption(AZ_solver,AZ_gmres);
  AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
  AztecOOSolver.SetAztecOption(AZ_overlap,Overlap);
  AztecOOSolver.SetAztecOption(AZ_graph_fill,ival);
  AztecOOSolver.SetAztecOption(AZ_poly_ord, ival);
  AztecOOSolver.SetAztecParam(AZ_drop, 0.0);
  AztecOOSolver.SetAztecParam(AZ_athresh, 0.0);
  AztecOOSolver.SetAztecParam(AZ_rthresh, 0.0);

  Epetra_MultiVector& RHS = *(Problem.GetRHS());
  Epetra_MultiVector& LHS = *(Problem.GetLHS());
  Teuchos::RefCountPtr<Epetra_RowMatrix> A = Teuchos::rcp(Problem.GetMatrix(), false);

  LHS.Random();
  A->Multiply(false,LHS,RHS);

  Teuchos::ParameterList List;
  List.set("fact: level-of-fill", ival);
  List.set("relaxation: sweeps", ival);
  List.set("relaxation: damping factor", 1.0);
  List.set("relaxation: zero starting solution", true);

  //default combine mode is as for AztecOO
  List.set("schwarz: combine mode", Zero);

  Epetra_Time Time(A->Comm());

  Teuchos::RefCountPtr<Ifpack_Preconditioner> Prec;

  if (what == "Jacobi") {
    Prec = Teuchos::rcp( new Ifpack_PointRelaxation(&*A) );
    List.set("relaxation: type", "Jacobi");
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_Jacobi);
    AztecOOSolver.SetAztecOption(AZ_reorder,0);
  }
  else if (what == "IC no reord") {
    Prec = Teuchos::rcp( new Ifpack_AdditiveSchwarz<Ifpack_IC>(&*A,Overlap) );
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    AztecOOSolver.SetAztecOption(AZ_subdomain_solve,AZ_icc);
    AztecOOSolver.SetAztecOption(AZ_reorder,0);
  }
  else if (what == "IC reord") {
    Prec = Teuchos::rcp( new Ifpack_AdditiveSchwarz<Ifpack_IC>(&*A,Overlap) );
    List.set("schwarz: use reordering", true);
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    AztecOOSolver.SetAztecOption(AZ_subdomain_solve,AZ_icc);
    AztecOOSolver.SetAztecOption(AZ_reorder,1);
  }
  else if (what == "ILU no reord") {
    Prec = Teuchos::rcp( new Ifpack_AdditiveSchwarz<Ifpack_ILU>(&*A,Overlap) );
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    AztecOOSolver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
    AztecOOSolver.SetAztecOption(AZ_reorder,0);
  }
  else if (what == "ILU reord") {
    Prec = Teuchos::rcp( new Ifpack_AdditiveSchwarz<Ifpack_ILU>(&*A,Overlap) );
    List.set("schwarz: use reordering", true);
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    AztecOOSolver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
    AztecOOSolver.SetAztecOption(AZ_reorder,1);
  }
#ifdef HAVE_IFPACK_AMESOS
  else if (what == "LU") {
    Prec = Teuchos::rcp( new Ifpack_AdditiveSchwarz<Ifpack_Amesos>(&*A,Overlap) );
    List.set("amesos: solver type", "Klu");
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    AztecOOSolver.SetAztecOption(AZ_subdomain_solve,AZ_lu);
  }
#endif
  else {
    cerr << "Option not recognized" << endl;
    exit(EXIT_FAILURE);
  }

  // ==================================== //
  // Solve with AztecOO's preconditioners //
  // ==================================== //

  LHS.PutScalar(0.0);

  Time.ResetStartTime();
  AztecOOSolver.Iterate(150,1e-5);

  if (verbose) {
    cout << endl;
    cout << "==================================================" << endl;
    cout << "Testing `" << what << "', Overlap = "
         << Overlap << ", ival = " << ival << endl;
    cout << endl;
    cout << "[AztecOO] Total time = " << Time.ElapsedTime() << " (s)" << endl;
    cout << "[AztecOO] Residual   = " << AztecOOSolver.TrueResidual() << " (s)" << endl;
    cout << "[AztecOO] Iterations = " << AztecOOSolver.NumIters() << endl;
    cout << endl;
  }

  int AztecOOPrecIters = AztecOOSolver.NumIters();

  // =========================================== //
  // Create the IFPACK preconditioner and solver //
  // =========================================== //

  Epetra_Time Time2(A->Comm());
  assert(Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->SetParameters(List));

  Time.ResetStartTime();
  IFPACK_CHK_ERR(Prec->Initialize());
  if (verbose)
    cout << "[IFPACK] Time for Initialize() = "
         << Time.ElapsedTime() << " (s)" << endl;

  Time.ResetStartTime();
  IFPACK_CHK_ERR(Prec->Compute());
  if (verbose)
    cout << "[IFPACK] Time for Compute() = "
         << Time.ElapsedTime() << " (s)" << endl;


  AztecOOSolver.SetPrecOperator(&*Prec);

  LHS.PutScalar(0.0);

  Time.ResetStartTime();
  AztecOOSolver.Iterate(150,1e-5);

  if (verbose) {
    cout << "[IFPACK] Total time = " << Time2.ElapsedTime() << " (s)" << endl;
    cout << "[IFPACK] Residual   = " << AztecOOSolver.TrueResidual() << " (s)" << endl;
    cout << "[IFPACK] Iterations = " << AztecOOSolver.NumIters() << endl;
    cout << endl;
  }

  int IFPACKPrecIters = AztecOOSolver.NumIters();

  if (IFPACK_ABS(AztecOOPrecIters - IFPACKPrecIters) > 3) {
    cerr << "TEST FAILED (" << AztecOOPrecIters << " != "
         << IFPACKPrecIters << ")" << endl;
    return(false);
  }
  else
    return(true);

}

// ======================================================================
int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int nx = 30;
  Teuchos::ParameterList GaleriList;
  GaleriList.set("n", nx * nx);
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx);

  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap("Linear", Comm, GaleriList) );
  Teuchos::RefCountPtr<Epetra_RowMatrix> A = Teuchos::rcp( Galeri::CreateCrsMatrix("Laplace2D", &*Map, GaleriList) );
  Epetra_Vector LHS(*Map);
  Epetra_Vector RHS(*Map);
  Epetra_LinearProblem Problem(&*A, &LHS, &RHS);

  int TestPassed = true;

  // Jacobi as in AztecOO (no overlap)
  for (int ival = 1 ; ival < 10 ; ival += 3) {
    TestPassed = TestPassed &&
      CompareWithAztecOO(Problem,"Jacobi",0,ival);
  }

#if 0
  // AztecOO with IC and overlap complains, also with
  // large fill-ins (in parallel)
  TestPassed = TestPassed &&
    CompareWithAztecOO(Problem,"IC no reord",0,0);
  TestPassed = TestPassed &&
    CompareWithAztecOO(Problem,"IC reord",0,0);

  vector<std::string> Tests;
  // now test solvers that accept overlap
  Tests.push_back("ILU no reord");
  Tests.push_back("ILU reord");
  // following requires --enable-aztecoo-azlu
#ifdef HAVE_IFPACK_AMESOS
  //Tests.push_back("LU");
#endif

  for (unsigned int i = 0 ; i < Tests.size() ; ++i) {
    for (int overlap = 0 ; overlap < 1 ; overlap += 2) {
      for (int ival = 0 ; ival < 10 ; ival += 4)
        TestPassed = TestPassed &&
          CompareWithAztecOO(Problem,Tests[i],overlap,ival);
    }
  }
#endif

  if (!TestPassed) {
    cerr << "Test `CompareWithAztecOO.exe' FAILED!" << endl;
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  cout << "Test `CompareWithAztecOO.exe' passed!" << endl;

  exit(EXIT_SUCCESS);
}
