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


#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialSpdDenseSolver.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Version.h"

// prototypes

int check(Epetra_SerialSpdDenseSolver & solver, double * A1, int LDA,
	  int N1, int NRHS1, double OneNorm1,
	  double * B1, int LDB1,
	  double * X1, int LDX1,
	  bool Upper, bool verbose);

void GenerateHilbert(double *A, int LDA, int N);

bool Residual( int N, int NRHS, double * A, int LDA,
	       double * X, int LDX, double * B, int LDB, double * resid);


int main(int argc, char *argv[])
{
  int ierr = 0, i, j, k;

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if(verbose && Comm.MyPID()==0)
    std::cout << Epetra_Version() << std::endl << std::endl;

  int rank = Comm.MyPID();
  //  char tmp;
  //  if (rank==0) std::cout << "Press any key to continue..."<< std::endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  if (verbose) std::cout << Comm << std::endl;

  //  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;

  int N = 20;
  int NRHS = 4;
  double * A = new double[N*N];
  double * A1 = new double[N*N];
  double * X = new double[(N+1)*NRHS];
  double * X1 = new double[(N+1)*NRHS];
  int LDX = N+1;
  int LDX1 = N+1;
  double * B = new double[N*NRHS];
  double * B1 = new double[N*NRHS];
  int LDB = N;
  int LDB1 = N;

  int LDA = N;
  int LDA1 = LDA;
  double OneNorm1;
  bool Upper = false;

  Epetra_SerialSpdDenseSolver solver;
  Epetra_SerialSymDenseMatrix * Matrix;
  for (int kk=0; kk<2; kk++) {
    for (i=1; i<=N; i++) {
      GenerateHilbert(A, LDA, i);
      OneNorm1 = 0.0;
      for (j=1; j<=i; j++) OneNorm1 += 1.0/((double) j); // 1-Norm = 1 + 1/2 + ...+1/n

      if (kk==0) {
	Matrix = new Epetra_SerialSymDenseMatrix(View, A, LDA, i);
	LDA1 = LDA;
      }
      else {
	Matrix = new Epetra_SerialSymDenseMatrix(Copy, A, LDA, i);
	LDA1 = i;
      }
      GenerateHilbert(A1, LDA1, i);
	
      if (kk==1) {
	solver.FactorWithEquilibration(true);
	Matrix->SetUpper();
	Upper = true;
	solver.SolveToRefinedSolution(false);
      }

      for (k=0; k<NRHS; k++)
	for (j=0; j<i; j++) {
	  B[j+k*LDB] = 1.0/((double) (k+3)*(j+3));
	  B1[j+k*LDB1] = B[j+k*LDB1];
	}
      Epetra_SerialDenseMatrix Epetra_B(View, B, LDB, i, NRHS);
      Epetra_SerialDenseMatrix Epetra_X(View, X, LDX, i, NRHS);
      solver.SetMatrix(*Matrix);
      solver.SetVectors(Epetra_X, Epetra_B);

      ierr = check(solver, A1, LDA1,  i, NRHS, OneNorm1, B1, LDB1,  X1, LDX1, Upper, verbose);
      assert (ierr>-1);
      delete Matrix;
      if (ierr!=0) {
	if (verbose) std::cout << "Factorization failed due to bad conditioning.  This is normal if SCOND is small."
			  << std::endl;
	break;
      }
    }
  }

  delete [] A;
  delete [] A1;
  delete [] X;
  delete [] X1;
  delete [] B;
  delete [] B1;

  /////////////////////////////////////////////////////////////////////
  // Now test norms and scaling functions
  /////////////////////////////////////////////////////////////////////

  Epetra_SerialSymDenseMatrix D;
  double ScalarA = 2.0;

  int DM = 10;
  int DN = 10;
  D.Shape(DM);
  for (j=0; j<DN; j++)
    for (i=0; i<DM; i++) D[j][i] = (double) (1+i+j*DM) ;

  //std::cout << D << std::endl;

  double NormInfD_ref = (double)(DM*(DN*(DN+1))/2);
  double NormOneD_ref = NormInfD_ref;

  double NormInfD = D.NormInf();
  double NormOneD = D.NormOne();

  if (verbose) {
    std::cout << " *** Before scaling *** " << std::endl
	 << " Computed one-norm of test matrix = " << NormOneD << std::endl
	 << " Expected one-norm                = " << NormOneD_ref << std::endl
	 << " Computed inf-norm of test matrix = " << NormInfD << std::endl
	 << " Expected inf-norm                = " << NormInfD_ref << std::endl;
  }
  D.Scale(ScalarA); // Scale entire D matrix by this value

  //std::cout << D << std::endl;

  NormInfD = D.NormInf();
  NormOneD = D.NormOne();
  if (verbose) {
    std::cout << " *** After scaling *** " << std::endl
	 << " Computed one-norm of test matrix = " << NormOneD << std::endl
	 << " Expected one-norm                = " << NormOneD_ref*ScalarA << std::endl
	 << " Computed inf-norm of test matrix = " << NormInfD << std::endl
	 << " Expected inf-norm                = " << NormInfD_ref*ScalarA << std::endl;
  }



  /////////////////////////////////////////////////////////////////////
  // Now test for larger system, both correctness and performance.
  /////////////////////////////////////////////////////////////////////


  N = 2000;
  NRHS = 5;
  LDA = N;
  LDB = N;
  LDX = N;

  if (verbose) std::cout << "\n\nComputing factor of an " << N << " x " << N << " SPD matrix...Please wait.\n\n" << std::endl;

  // Define A and X

  A = new double[LDA*N];
  X = new double[LDB*NRHS];

  for (j=0; j<N; j++) {
    for (k=0; k<NRHS; k++) X[j+k*LDX] = 1.0/((double) (j+5+k));
    for (i=0; i<N; i++) {
      if (i==j) A[i+j*LDA] = 100.0 + i;
      else A[i+j*LDA] = -1.0/((double) (i+10)*(j+10));
    }
  }

  // Define Epetra_SerialDenseMatrix object

  Epetra_SerialSymDenseMatrix BigMatrix(Copy, A, LDA, N);
  Epetra_SerialSymDenseMatrix OrigBigMatrix(View, A, LDA, N);

  Epetra_SerialSpdDenseSolver BigSolver;
  BigSolver.FactorWithEquilibration(true);
  BigSolver.SetMatrix(BigMatrix);

  // Time factorization

  Epetra_Flops counter;
  BigSolver.SetFlopCounter(counter);
  Epetra_Time Timer(Comm);
  double tstart = Timer.ElapsedTime();
  ierr = BigSolver.Factor();
  if (ierr!=0 && verbose) std::cout << "Error in factorization = "<<ierr<< std::endl;
  assert(ierr==0);
  double time = Timer.ElapsedTime() - tstart;

  double FLOPS = counter.Flops();
  double MFLOPS = FLOPS/time/1000000.0;
  if (verbose) std::cout << "MFLOPS for Factorization = " << MFLOPS << std::endl;

  // Define Left hand side and right hand side
  Epetra_SerialDenseMatrix LHS(View, X, LDX, N, NRHS);
  Epetra_SerialDenseMatrix RHS;
  RHS.Shape(N,NRHS); // Allocate RHS

  // Compute RHS from A and X

  Epetra_Flops RHS_counter;
  RHS.SetFlopCounter(RHS_counter);
  tstart = Timer.ElapsedTime();
  RHS.Multiply('L', 1.0, OrigBigMatrix, LHS, 0.0); // Symmetric Matrix-multiply
  time = Timer.ElapsedTime() - tstart;

  Epetra_SerialDenseMatrix OrigRHS = RHS;

  FLOPS = RHS_counter.Flops();
  MFLOPS = FLOPS/time/1000000.0;
  if (verbose) std::cout << "MFLOPS to build RHS (NRHS = " << NRHS <<") = " << MFLOPS << std::endl;

  // Set LHS and RHS and solve
  BigSolver.SetVectors(LHS, RHS);

  tstart = Timer.ElapsedTime();
  ierr = BigSolver.Solve();
  if (ierr==1 && verbose) std::cout << "LAPACK guidelines suggest this matrix might benefit from equilibration." << std::endl;
  else if (ierr!=0 && verbose) std::cout << "Error in solve = "<<ierr<< std::endl;
  assert(ierr>=0);
  time = Timer.ElapsedTime() - tstart;

  FLOPS = BigSolver.Flops();
  MFLOPS = FLOPS/time/1000000.0;
  if (verbose) std::cout << "MFLOPS for Solve (NRHS = " << NRHS <<") = " << MFLOPS << std::endl;

  double * resid = new double[NRHS];
  bool OK = Residual(N, NRHS, A, LDA, BigSolver.X(), BigSolver.LDX(),
		     OrigRHS.A(), OrigRHS.LDA(), resid);

  if (verbose) {
    if (!OK) std::cout << "************* Residual do not meet tolerance *************" << std::endl;
    for (i=0; i<NRHS; i++)
      std::cout << "Residual[" << i <<"] = "<< resid[i] << std::endl;
    std::cout  << std::endl;
  }

  // Solve again using the Epetra_SerialDenseVector class for LHS and RHS

  Epetra_SerialDenseVector X2;
  Epetra_SerialDenseVector B2;
  X2.Size(BigMatrix.N());
  B2.Size(BigMatrix.M());
  int length = BigMatrix.N();
  {for (int kk=0; kk<length; kk++) X2[kk] = ((double ) kk)/ ((double) length);} // Define entries of X2

  RHS_counter.ResetFlops();
  B2.SetFlopCounter(RHS_counter);
  tstart = Timer.ElapsedTime();
  B2.Multiply('N', 'N', 1.0, OrigBigMatrix, X2, 0.0); // Define B2 = A*X2
  time = Timer.ElapsedTime() - tstart;

  Epetra_SerialDenseVector OrigB2 = B2;

  FLOPS = RHS_counter.Flops();
  MFLOPS = FLOPS/time/1000000.0;
  if (verbose) std::cout << "MFLOPS to build single RHS = " << MFLOPS << std::endl;

  // Set LHS and RHS and solve
  BigSolver.SetVectors(X2, B2);

  tstart = Timer.ElapsedTime();
  ierr = BigSolver.Solve();
  time = Timer.ElapsedTime() - tstart;
  if (ierr==1 && verbose) std::cout << "LAPACK guidelines suggest this matrix might benefit from equilibration." << std::endl;
  else if (ierr!=0 && verbose) std::cout << "Error in solve = "<<ierr<< std::endl;
  assert(ierr>=0);

  FLOPS = counter.Flops();
  MFLOPS = FLOPS/time/1000000.0;
  if (verbose) std::cout << "MFLOPS to solve single RHS = " << MFLOPS << std::endl;

  OK = Residual(N, 1, A, LDA, BigSolver.X(), BigSolver.LDX(), OrigB2.A(),
		OrigB2.LDA(), resid);

  if (verbose) {
    if (!OK) std::cout << "************* Residual do not meet tolerance *************" << std::endl;
      std::cout << "Residual = "<< resid[0] << std::endl;
  }
  delete [] resid;
  delete [] A;
  delete [] X;

  ///////////////////////////////////////////////////
  // Now test default constructor and index operators
  ///////////////////////////////////////////////////

  N = 5;
  Epetra_SerialSymDenseMatrix C; // Implicit call to default constructor, should not need to call destructor
  C.Shape(5); // Make it 5 by 5
  double * C1 = new double[N*N];
  GenerateHilbert(C1, N, N); // Generate Hilber matrix

  C1[1+2*N] = 1000.0;  // Make matrix nonsymmetric

  // Fill values of C with Hilbert values
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      C(i,j) = C1[i+j*N];

  // Test if values are correctly written and read
  for (i=0; i<N; i++)
    for (j=0; j<N; j++) {
      assert(C(i,j) == C1[i+j*N]);
      assert(C(i,j) == C[j][i]);
    }

  if (verbose)
    std::cout << "Default constructor and index operator check OK.  Values of Hilbert matrix = "
	 << std::endl << C << std::endl
	 << "Values should be 1/(i+j+1), except value (1,2) should be 1000" << std::endl;

  delete [] C1;


#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}

int check(Epetra_SerialSpdDenseSolver &solver, double * A1, int LDA1,
	  int N1, int NRHS1, double OneNorm1,
	  double * B1, int LDB1,
	  double * X1, int LDX1,
	  bool Upper, bool verbose) {
  (void)OneNorm1;
  int i;
  bool OK;
  // Test query functions

  int M= solver.M();
  if (verbose) std::cout << "\n\nNumber of Rows = " << M << std::endl<< std::endl;
  assert(M==N1);

  int N= solver.N();
  if (verbose) std::cout << "\n\nNumber of Equations = " << N << std::endl<< std::endl;
  assert(N==N1);

  int LDA = solver.LDA();
  if (verbose) std::cout << "\n\nLDA = " << LDA << std::endl<< std::endl;
  assert(LDA==LDA1);

  int LDB = solver.LDB();
  if (verbose) std::cout << "\n\nLDB = " << LDB << std::endl<< std::endl;
  assert(LDB==LDB1);

  int LDX = solver.LDX();
  if (verbose) std::cout << "\n\nLDX = " << LDX << std::endl<< std::endl;
  assert(LDX==LDX1);

  int NRHS = solver.NRHS();
  if (verbose) std::cout << "\n\nNRHS = " << NRHS << std::endl<< std::endl;
  assert(NRHS==NRHS1);

  assert(solver.ANORM()==-1.0);
  assert(solver.RCOND()==-1.0);
  if (!solver.A_Equilibrated() && !solver.B_Equilibrated()) {
    assert(solver.SCOND()==-1.0);
    assert(solver.AMAX()==-1.0);
  }


  // Other binary tests

  assert(!solver.Factored());
  assert(solver.SymMatrix()->Upper()==Upper);
  assert(!solver.SolutionErrorsEstimated());
  assert(!solver.Inverted());
  assert(!solver.ReciprocalConditionEstimated());
  assert(!solver.Solved());

  assert(!solver.SolutionRefined());

  //std::cout << "Matrix before factorization " << std::endl << *solver.SymMatrix() << std::endl << std::endl;

  int ierr = solver.Factor();
  //std::cout << "Matrix after factorization " << std::endl << *solver.SymMatrix() << std::endl << std::endl;
  //std::cout << "Factor after factorization " << std::endl << *solver.SymFactoredMatrix() << std::endl << std::endl;
  assert(ierr>-1);
  if (ierr!=0) return(ierr); // Factorization failed due to poor conditioning.
  double rcond;
  ierr = solver.ReciprocalConditionEstimate(rcond);
  assert(ierr==0);
  if (verbose) {

    double rcond1 = 1.0/std::exp(3.5*((double)N));
    if (N==1) rcond1 = 1.0;
    std::cout << "\n\nSCOND = "<< rcond << " should be approx = "
		    << rcond1 << std::endl << std::endl;
  }

  ierr = solver.Solve();
  assert(ierr>-1);
  if (ierr!=0 && verbose) std::cout << "LAPACK rules suggest system should be equilibrated." << std::endl;

  assert(solver.Factored());
  assert(solver.SymMatrix()->Upper()==Upper);
  assert(solver.ReciprocalConditionEstimated());
  assert(solver.Solved());

  if (solver.SolutionErrorsEstimated()) {
    if (verbose) {
      std::cout << "\n\nFERR[0] = "<< solver.FERR()[0] << std::endl;
      std::cout << "\n\nBERR[0] = "<< solver.BERR()[0] << std::endl<< std::endl;
    }
  }

  double * resid = new double[NRHS];
  OK = Residual(N, NRHS, A1, LDA1, solver.X(), solver.LDX(), B1, LDB1, resid);
  if (verbose) {
    if (!OK) std::cout << "************* Residual do not meet tolerance *************" << std::endl;
    std::cout << "\n\nResiduals using factorization to solve" << std::endl;
    for (i=0; i<NRHS; i++)
      std::cout << "Residual[" << i <<"] = "<< resid[i] << std::endl;
    std::cout  << std::endl;
  }


  ierr = solver.Invert();
  assert(ierr>-1);

  assert(solver.Inverted());
  assert(!solver.Factored());

  Epetra_SerialDenseMatrix RHS1(Copy, B1, LDB1, N, NRHS);
  Epetra_SerialDenseMatrix LHS1(Copy, X1, LDX1, N, NRHS);
  assert(solver.SetVectors(LHS1, RHS1)==0);
  assert(!solver.Solved());

  assert(solver.Solve()>-1);
	


  OK = Residual(N, NRHS, A1, LDA1, solver.X(), solver.LDX(), B1, LDB1, resid);

  if (verbose) {
    if (!OK) std::cout << "************* Residual do not meet tolerance *************" << std::endl;
    std::cout << "Residuals using inverse to solve" << std::endl;
    for (i=0; i<NRHS; i++)
      std::cout << "Residual[" << i <<"] = "<< resid[i] << std::endl;
    std::cout  << std::endl;
  }
  delete [] resid;


  return(0);
}

 void GenerateHilbert(double *A, int LDA, int N) {
   for (int j=0; j<N; j++)
     for (int i=0; i<N; i++)
       A[i+j*LDA] = 1.0/((double)(i+j+1));
   return;
 }

bool Residual( int N, int NRHS, double * A, int LDA,
	       double * X, int LDX, double * B, int LDB, double * resid) {

  Epetra_BLAS Blas;
  char Transa = 'N';
  Blas.GEMM(Transa, 'N', N, NRHS, N, -1.0, A, LDA,
	    X, LDX, 1.0, B, LDB);
  bool OK = true;
  for (int i=0; i<NRHS; i++) {
    resid[i] = Blas.NRM2(N, B+i*LDB);
    if (resid[i]>1.0E-7) OK = false;
  }

  return(OK);
}
