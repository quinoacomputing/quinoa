/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"
#include "CTrilinos_enums.h"
#include "CTrilinos_flex_enums.h"

#include "CEpetra_Map.h"
#include "CEpetra_Time.h"
#include "CEpetra_CrsMatrix.h"
#include "CEpetra_JadMatrix.h"
#include "CEpetra_Vector.h"
#include "CEpetra_Flops.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_Operator.h"
#include "CEpetra_DistObject.h"
#include "CEpetra_Comm.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_RowMatrix.h"
#include "CEpetra_CompObject.h"
#include "CEpetra_SrcDistObject.h"
#ifdef HAVE_MPI
#include "CEpetra_MpiComm.h"
#include "mpi.h"
#else
#include "CEpetra_SerialComm.h"
#endif
#include "../../../../epetra/test/epetra_test_err.h"
#include "Epetra_Version.h"
#include <vector>
#include <algorithm>
#include <string>

// prototypes

int checkValues( double x, double y, string message = "", bool verbose = false) { 
  if (fabs((x-y)/x) > 0.01) {
    return(1); 
    if (verbose) cout << "********** " << message << " check failed.********** " << endl;
  }
  else {
    if (verbose) cout << message << " check OK." << endl;    
    return(0);
  }
}

int checkMultiVectors( CT_Epetra_MultiVector_ID_t & X, CT_Epetra_MultiVector_ID_t & Y, string message = "", bool verbose = false) {
  int numVectors = Epetra_MultiVector_NumVectors(X);
  int length = Epetra_MultiVector_MyLength(Y);
  int badvalue = 0;
  int globalbadvalue = 0;
  for (int j=0; j<numVectors; j++) {
    CT_Epetra_Vector_ID_t vecx = Epetra_MultiVector_getVector(X, j);
    CT_Epetra_Vector_ID_t vecy = Epetra_MultiVector_getVector(Y, j);
    for (int i=0; i< length; i++)
      if (checkValues(Epetra_Vector_getElement(vecx, i), Epetra_Vector_getElement(vecy, i))==1) badvalue = 1;
    Epetra_Vector_Destroy(&vecy);
    Epetra_Vector_Destroy(&vecx);
  }
  CT_Epetra_Comm_ID_t Comm = Epetra_DistObject_Comm(Epetra_DistObject_Degeneralize(
      Epetra_MultiVector_Generalize(X)));
  Epetra_Comm_MaxAll_Int(Comm, &badvalue, &globalbadvalue, 1);
  Epetra_Comm_Destroy(&Comm);

  if (verbose) {
    if (globalbadvalue==0) cout << message << " check OK." << endl;
    else cout << "********* " << message << " check failed.********** " << endl;
  }
  return(globalbadvalue);
}

int check(CT_Epetra_RowMatrix_ID_t & A, CT_Epetra_RowMatrix_ID_t & B, bool verbose);

int powerMethodTests(CT_Epetra_RowMatrix_ID_t & A, CT_Epetra_RowMatrix_ID_t & JadA,
                     CT_Epetra_Map_ID_Flex_t & Map, CT_Epetra_Vector_ID_Flex_t & q,
                     CT_Epetra_Vector_ID_Flex_t & z, CT_Epetra_Vector_ID_Flex_t & resid, bool verbose);

int power_method(boolean TransA, CT_Epetra_RowMatrix_ID_t& A, CT_Epetra_Vector_ID_Flex_t& q,
                 CT_Epetra_Vector_ID_Flex_t& z0, CT_Epetra_Vector_ID_Flex_t& resid, double* lambda,
                 int niters, double tolerance, bool verbose);

int main(int argc, char *argv[])
{
  int ierr = 0, i, forierr = 0;
#ifdef HAVE_MPI

  CT_Epetra_MpiComm_ID_Flex_t Comm;

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int rank; // My process ID

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Comm.Epetra_MpiComm = Epetra_MpiComm_Create( MPI_COMM_WORLD );

#else

  int rank = 0;
  CT_Epetra_SerialComm_ID_Flex_t Comm;
  Comm.Epetra_SerialComm = Epetra_SerialComm_Create();

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int verbose_int = verbose ? 1 : 0;
  Epetra_Comm_Broadcast_Int(Comm.Epetra_Comm, &verbose_int, 1, 0);
  verbose = verbose_int==1 ? true : false;


  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Epetra_Comm_Barrier(Comm.Epetra_Comm);

//  Epetra_Comm_SetTracebackMode(Comm.Epetra_Comm, 0); // This should shut down any error traceback reporting
  int MyPID = Epetra_Comm_MyPID(Comm.Epetra_Comm);
  int NumProc = Epetra_Comm_NumProc(Comm.Epetra_Comm);

  if(verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
		    << " is alive."<<endl;

  // Redefine verbose to only print on PE 0
  if(verbose && rank!=0) 
		verbose = false;

  int NumMyEquations = 10000;
  int NumGlobalEquations = (NumMyEquations * NumProc) + EPETRA_MIN(NumProc,3);
  if(MyPID < 3) 
    NumMyEquations++;

  // Construct a Map that puts approximately the same Number of equations on each processor

  CT_Epetra_Map_ID_Flex_t Map;
  Map.Epetra_Map = Epetra_Map_Create_Linear(NumGlobalEquations, NumMyEquations, 0, Comm.Epetra_Comm);
  
  // Get update list and number of local equations from newly created Map
  vector<int> MyGlobalElements(Epetra_BlockMap_NumMyElements(Map.Epetra_BlockMap));
  Epetra_BlockMap_MyGlobalElements_Fill(Map.Epetra_BlockMap, &MyGlobalElements[0]);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor

  vector<int> NumNz(NumMyEquations);

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for(i = 0; i < NumMyEquations; i++)
    if((MyGlobalElements[i] == 0) || (MyGlobalElements[i] == NumGlobalEquations - 1))
      NumNz[i] = 1;
    else
      NumNz[i] = 2;

  // Create a Epetra_Matrix

  CT_Epetra_CrsMatrix_ID_Flex_t A;
  A.Epetra_CrsMatrix = Epetra_CrsMatrix_Create_VarPerRow(
      CT_Epetra_DataAccess_E_Copy, Map.Epetra_Map, &NumNz[0], FALSE);
  EPETRA_TEST_ERR(Epetra_CrsMatrix_IndicesAreGlobal(A.Epetra_CrsMatrix),ierr);
  EPETRA_TEST_ERR(Epetra_CrsMatrix_IndicesAreLocal(A.Epetra_CrsMatrix),ierr);
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  vector<double> Values(2);
  Values[0] = -1.0; 
	Values[1] = -1.0;
	vector<int> Indices(2);
  double two = 2.0;
  int NumEntries;

  forierr = 0;
  for(i = 0; i < NumMyEquations; i++) {
    if(MyGlobalElements[i] == 0) {
			Indices[0] = 1;
			NumEntries = 1;
		}
    else if (MyGlobalElements[i] == NumGlobalEquations-1) {
			Indices[0] = NumGlobalEquations-2;
			NumEntries = 1;
		}
    else {
			Indices[0] = MyGlobalElements[i]-1;
			Indices[1] = MyGlobalElements[i]+1;
			NumEntries = 2;
		}
		forierr += !(Epetra_CrsMatrix_InsertGlobalValues(A.Epetra_CrsMatrix,
                             MyGlobalElements[i], NumEntries, &Values[0], &Indices[0])==0);
		forierr += !(Epetra_CrsMatrix_InsertGlobalValues(A.Epetra_CrsMatrix,
                             MyGlobalElements[i], 1, &two, &MyGlobalElements[i])>0); // Put in the diagonal entry
  }
  EPETRA_TEST_ERR(forierr,ierr);

  // Finish up
  Epetra_CrsMatrix_FillComplete(A.Epetra_CrsMatrix, TRUE);
  Epetra_CrsMatrix_OptimizeStorage(A.Epetra_CrsMatrix);

  CT_Epetra_JadMatrix_ID_Flex_t JadA, JadA1, JadA2;
  JadA.Epetra_JadMatrix = Epetra_JadMatrix_Create(A.Epetra_RowMatrix);
  JadA1.Epetra_JadMatrix = Epetra_JadMatrix_Create(A.Epetra_RowMatrix);
  JadA2.Epetra_JadMatrix = Epetra_JadMatrix_Create(A.Epetra_RowMatrix);

  // Create vectors for Power method

  CT_Epetra_Vector_ID_Flex_t q, z, resid;
  q.Epetra_Vector = Epetra_Vector_Create(Map.Epetra_BlockMap, TRUE);
  z.Epetra_Vector = Epetra_Vector_Create(Map.Epetra_BlockMap, TRUE);
  Epetra_MultiVector_Random(z.Epetra_MultiVector);
  resid.Epetra_Vector = Epetra_Vector_Create(Map.Epetra_BlockMap, TRUE);

  CT_Epetra_Flops_ID_t flopcounter = Epetra_Flops_Create();
  Epetra_CompObject_SetFlopCounter(A.Epetra_CompObject, flopcounter);

  Epetra_CompObject_SetFlopCounter_Matching(q.Epetra_CompObject, A.Epetra_CompObject);
  Epetra_CompObject_SetFlopCounter_Matching(z.Epetra_CompObject, A.Epetra_CompObject);
  Epetra_CompObject_SetFlopCounter_Matching(resid.Epetra_CompObject, A.Epetra_CompObject);

  Epetra_CompObject_SetFlopCounter_Matching(JadA.Epetra_CompObject, A.Epetra_CompObject);
  Epetra_CompObject_SetFlopCounter_Matching(JadA1.Epetra_CompObject, A.Epetra_CompObject);
  Epetra_CompObject_SetFlopCounter_Matching(JadA2.Epetra_CompObject, A.Epetra_CompObject);
  

  if (verbose) cout << "=======================================" << endl
		    << "Testing Jad using CrsMatrix as input..." << endl
		    << "=======================================" << endl;

  Epetra_CompObject_ResetFlops(A.Epetra_CompObject);
  powerMethodTests(A.Epetra_RowMatrix, JadA.Epetra_RowMatrix, Map, q, z, resid, verbose);

  // Increase diagonal dominance

  if (verbose) cout << "\n\nIncreasing the magnitude of first diagonal term and solving again\n\n"
		    << endl;

  
  if (Epetra_CrsMatrix_MyGlobalRow(A.Epetra_CrsMatrix, 0)) {
    int numvals = Epetra_CrsMatrix_NumGlobalEntries(A.Epetra_CrsMatrix, 0);
    vector<double> Rowvals(numvals);
    vector<int> Rowinds(numvals);
    Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices(A.Epetra_CrsMatrix, 0,
        numvals, &numvals, &Rowvals[0], &Rowinds[0]); // Get A[0,0]

    for (i=0; i<numvals; i++) if (Rowinds[i] == 0) Rowvals[i] *= 10.0;
    
    Epetra_CrsMatrix_ReplaceGlobalValues(A.Epetra_CrsMatrix, 0, numvals, &Rowvals[0], &Rowinds[0]);
  }
  Epetra_JadMatrix_UpdateValues(JadA.Epetra_JadMatrix, A.Epetra_RowMatrix, FALSE);
  Epetra_CompObject_ResetFlops(A.Epetra_CompObject);
  powerMethodTests(A.Epetra_RowMatrix, JadA.Epetra_RowMatrix, Map, q, z, resid, verbose);

  if (verbose) cout << "================================================================" << endl
		          << "Testing Jad using Jad matrix as input matrix for construction..." << endl
		          << "================================================================" << endl;
  Epetra_CompObject_ResetFlops(JadA1.Epetra_CompObject);
  powerMethodTests(JadA1.Epetra_RowMatrix, JadA2.Epetra_RowMatrix, Map, q, z, resid, verbose);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

return ierr ;
}

int powerMethodTests(CT_Epetra_RowMatrix_ID_t & A, CT_Epetra_RowMatrix_ID_t & JadA,
                     CT_Epetra_Map_ID_Flex_t & Map, CT_Epetra_Vector_ID_Flex_t & q,
                     CT_Epetra_Vector_ID_Flex_t & z, CT_Epetra_Vector_ID_Flex_t & resid, bool verbose) {

  // variable needed for iteration
  double lambda = 0.0;
  // int niters = 10000;
  int niters = 300;
  double tolerance = 1.0e-2;
  int ierr = 0;

  /////////////////////////////////////////////////////////////////////////////////////////////////
	
  // Iterate

  CT_Epetra_Comm_ID_t Comm = Epetra_BlockMap_Comm(Map.Epetra_BlockMap);
  CT_Epetra_Time_ID_t timer = Epetra_Time_Create(Comm);
  Epetra_Comm_Destroy(&Comm);
	
  double startTime = Epetra_Time_ElapsedTime(timer);
  EPETRA_TEST_ERR(power_method(FALSE, A, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  double elapsed_time = Epetra_Time_ElapsedTime(timer) - startTime;
  double total_flops = Epetra_CompObject_Flops(q.Epetra_CompObject);
  double MFLOPs = total_flops/elapsed_time/1000000.0;
  double lambdaref = lambda;
  double flopsref = total_flops;

  if (verbose) 
	  cout << "\n\nTotal MFLOPs for reference first solve = " << MFLOPs << endl
		  <<     "Total FLOPS                            = " <<total_flops <<endl<<endl;

  lambda = 0.0;
  startTime = Epetra_Time_ElapsedTime(timer);
  EPETRA_TEST_ERR(power_method(FALSE, JadA, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  elapsed_time = Epetra_Time_ElapsedTime(timer) - startTime;
  total_flops = Epetra_CompObject_Flops(q.Epetra_CompObject);
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) 
	  cout << "\n\nTotal MFLOPs for candidate first solve = " << MFLOPs << endl
		  <<     "Total FLOPS                            = " <<total_flops <<endl<<endl;

  EPETRA_TEST_ERR(checkValues(lambda,lambdaref," No-transpose Power Method result", verbose),ierr);
  EPETRA_TEST_ERR(checkValues(total_flops,flopsref," No-transpose Power Method flop count", verbose),ierr);

  /////////////////////////////////////////////////////////////////////////////////////////////////
	
  // Solve transpose problem

  if (verbose) cout << "\n\nUsing transpose of matrix and solving again (should give same result).\n\n"
		    << endl;
  // Iterate
  lambda = 0.0;
  startTime = Epetra_Time_ElapsedTime(timer);
  EPETRA_TEST_ERR(power_method(TRUE, A, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  elapsed_time = Epetra_Time_ElapsedTime(timer) - startTime;
  total_flops = Epetra_CompObject_Flops(q.Epetra_CompObject);
  MFLOPs = total_flops/elapsed_time/1000000.0;
  lambdaref = lambda;
  flopsref = total_flops;

  if (verbose) 
	 cout << "\n\nTotal MFLOPs for reference transpose solve = " << MFLOPs << endl
		 <<     "Total FLOPS                                = " <<total_flops <<endl<<endl;

  lambda = 0.0;
  startTime = Epetra_Time_ElapsedTime(timer);
  EPETRA_TEST_ERR(power_method(TRUE, JadA, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  elapsed_time = Epetra_Time_ElapsedTime(timer) - startTime;
  total_flops = Epetra_CompObject_Flops(q.Epetra_CompObject);
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) 
	  cout << "\n\nTotal MFLOPs for candidate transpose solve = " << MFLOPs << endl
		  <<     "Total FLOPS                                = " <<total_flops <<endl<<endl;

  EPETRA_TEST_ERR(checkValues(lambda,lambdaref,"Transpose Power Method result", verbose),ierr);
  EPETRA_TEST_ERR(checkValues(total_flops,flopsref,"Transpose Power Method flop count", verbose),ierr);

  EPETRA_TEST_ERR(check(A, JadA, verbose),ierr);

  return(0);
}
int power_method(boolean TransA, CT_Epetra_RowMatrix_ID_t& A, CT_Epetra_Vector_ID_Flex_t& q,
                 CT_Epetra_Vector_ID_Flex_t& z0, CT_Epetra_Vector_ID_Flex_t& resid, double* lambda,
                 int niters, double tolerance, bool verbose) 
{  
	
  // Fill z with random Numbers
  CT_Epetra_Vector_ID_Flex_t z;
  z.Epetra_Vector = Epetra_Vector_Duplicate(z0.Epetra_Vector);
	
  // variable needed for iteration
  double normz, residual;

  int ierr = 1;
	
  for(int iter = 0; iter < niters; iter++) {
		Epetra_MultiVector_Norm2(z.Epetra_MultiVector, &normz); // Compute 2-norm of z
		Epetra_MultiVector_Scale(q.Epetra_MultiVector, 1.0/normz, z.Epetra_MultiVector);
		Epetra_RowMatrix_Multiply(A, TransA, q.Epetra_MultiVector, z.Epetra_MultiVector); // Compute z = A*q // SEGFAULT HAPPENS HERE
		Epetra_MultiVector_Dot(q.Epetra_MultiVector, z.Epetra_MultiVector, lambda); // Approximate maximum eigenvaluE
		if(iter%100==0 || iter+1==niters) {
			Epetra_MultiVector_Update_WithAB(resid.Epetra_MultiVector, 1.0, z.Epetra_MultiVector, -(*lambda), q.Epetra_MultiVector, 0.0); // Compute A*q - lambda*q
			Epetra_MultiVector_Norm2(resid.Epetra_MultiVector, &residual);
			if(verbose) cout << "Iter = " << iter << "  Lambda = " << *lambda 
											 << "  Residual of A*q - lambda*q = " << residual << endl;
		}
		if(residual < tolerance) {
			ierr = 0;
			break;
		}
	}

  Epetra_Vector_Destroy(&z.Epetra_Vector);

  return(ierr);
}

int check(CT_Epetra_RowMatrix_ID_t& A, CT_Epetra_RowMatrix_ID_t & B, bool verbose)  {  

  int ierr = 0;

  CT_Epetra_Operator_ID_t tmp_oA = Epetra_Operator_Degeneralize(Epetra_RowMatrix_Generalize(A));
  CT_Epetra_Operator_ID_t tmp_oB = Epetra_Operator_Degeneralize(Epetra_RowMatrix_Generalize(B));

  CT_Epetra_Comm_ID_t CommA = Epetra_Operator_Comm(tmp_oA);
  CT_Epetra_Comm_ID_t CommB = Epetra_Operator_Comm(tmp_oB);
  EPETRA_TEST_ERR(!Epetra_Comm_NumProc(CommA)==Epetra_Comm_NumProc(CommB),ierr);
  EPETRA_TEST_ERR(!Epetra_Comm_MyPID(CommA)==Epetra_Comm_MyPID(CommB),ierr);
  Epetra_Comm_Destroy(&CommB);
  Epetra_Comm_Destroy(&CommA);

  EPETRA_TEST_ERR(!Epetra_RowMatrix_Filled(A)==Epetra_RowMatrix_Filled(B),ierr);
  EPETRA_TEST_ERR(!Epetra_Operator_HasNormInf(tmp_oA)==Epetra_Operator_HasNormInf(tmp_oB),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_LowerTriangular(A)==Epetra_RowMatrix_LowerTriangular(B),ierr);

  CT_Epetra_BlockMap_ID_t mapA = Epetra_SrcDistObject_Map(Epetra_SrcDistObject_Degeneralize(
      Epetra_RowMatrix_Generalize(A)));
  CT_Epetra_BlockMap_ID_t mapB = Epetra_SrcDistObject_Map(Epetra_SrcDistObject_Degeneralize(
      Epetra_RowMatrix_Generalize(B)));
  EPETRA_TEST_ERR(!Epetra_BlockMap_SameAs(mapA, mapB),ierr);
  Epetra_BlockMap_Destroy(&mapB);
  Epetra_BlockMap_Destroy(&mapA);

  EPETRA_TEST_ERR(!Epetra_RowMatrix_MaxNumEntries(A)==Epetra_RowMatrix_MaxNumEntries(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumGlobalCols(A)==Epetra_RowMatrix_NumGlobalCols(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumGlobalDiagonals(A)==Epetra_RowMatrix_NumGlobalDiagonals(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumGlobalNonzeros(A)==Epetra_RowMatrix_NumGlobalNonzeros(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumGlobalRows(A)==Epetra_RowMatrix_NumGlobalRows(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumMyCols(A)==Epetra_RowMatrix_NumMyCols(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumMyDiagonals(A)==Epetra_RowMatrix_NumMyDiagonals(B),ierr);
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumMyNonzeros(A)==Epetra_RowMatrix_NumMyNonzeros(B),ierr);
  for (int i=0; i<Epetra_RowMatrix_NumMyRows(A); i++) {
    int nA, nB;
    Epetra_RowMatrix_NumMyRowEntries(A,i,&nA);
    Epetra_RowMatrix_NumMyRowEntries(B,i,&nB);
    EPETRA_TEST_ERR(!nA==nB,ierr);
  }
  EPETRA_TEST_ERR(!Epetra_RowMatrix_NumMyRows(A)==Epetra_RowMatrix_NumMyRows(B),ierr);

  CT_Epetra_BlockMap_ID_t bodmA = Epetra_BlockMap_Degeneralize(Epetra_Map_Generalize(
      Epetra_Operator_OperatorDomainMap(tmp_oA)));
  CT_Epetra_BlockMap_ID_t bodmB = Epetra_BlockMap_Degeneralize(Epetra_Map_Generalize(
      Epetra_Operator_OperatorDomainMap(tmp_oB)));
  EPETRA_TEST_ERR(!Epetra_BlockMap_SameAs(bodmA, bodmB),ierr);

  CT_Epetra_BlockMap_ID_t bormA = Epetra_BlockMap_Degeneralize(Epetra_Map_Generalize(
      Epetra_Operator_OperatorRangeMap(tmp_oA)));
  CT_Epetra_BlockMap_ID_t bormB = Epetra_BlockMap_Degeneralize(Epetra_Map_Generalize(
      Epetra_Operator_OperatorRangeMap(tmp_oB)));
  EPETRA_TEST_ERR(!Epetra_BlockMap_SameAs(bormA, bormB),ierr);

  CT_Epetra_BlockMap_ID_t brcmA = Epetra_BlockMap_Degeneralize(Epetra_Map_Generalize(
      Epetra_RowMatrix_RowMatrixColMap(A)));
  CT_Epetra_BlockMap_ID_t brcmB = Epetra_BlockMap_Degeneralize(Epetra_Map_Generalize(
      Epetra_RowMatrix_RowMatrixColMap(B)));
  EPETRA_TEST_ERR(!Epetra_BlockMap_SameAs(brcmA, brcmB),ierr);

  CT_Epetra_BlockMap_ID_t brrmA = Epetra_BlockMap_Degeneralize(Epetra_Map_Generalize(
      Epetra_RowMatrix_RowMatrixRowMap(A)));
  CT_Epetra_BlockMap_ID_t brrmB = Epetra_BlockMap_Degeneralize(Epetra_Map_Generalize(
      Epetra_RowMatrix_RowMatrixRowMap(B)));
  EPETRA_TEST_ERR(!Epetra_BlockMap_SameAs(brrmA, brrmB),ierr);

  EPETRA_TEST_ERR(!Epetra_RowMatrix_UpperTriangular(A)==Epetra_RowMatrix_UpperTriangular(B),ierr);
  EPETRA_TEST_ERR(!Epetra_Operator_UseTranspose(tmp_oA)==Epetra_Operator_UseTranspose(tmp_oB),ierr);

  int NumVectors = 5;
  { // No transpose case
    CT_Epetra_MultiVector_ID_t X = Epetra_MultiVector_Create(bodmA, NumVectors, TRUE);
    CT_Epetra_MultiVector_ID_t YA1 = Epetra_MultiVector_Create(bormA, NumVectors, TRUE);
    CT_Epetra_MultiVector_ID_t YA2 = Epetra_MultiVector_Duplicate(YA1);
    CT_Epetra_MultiVector_ID_t YB1 = Epetra_MultiVector_Duplicate(YA1);
    CT_Epetra_MultiVector_ID_t YB2 = Epetra_MultiVector_Duplicate(YA1);
    Epetra_MultiVector_Random(X);
    
    boolean transA = FALSE;
    Epetra_Operator_SetUseTranspose(tmp_oA, transA);
    Epetra_Operator_SetUseTranspose(tmp_oB, transA);
    Epetra_Operator_Apply(tmp_oA,X,YA1);
    Epetra_RowMatrix_Multiply(A, transA, X, YA2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YA2,"A Multiply and A Apply", verbose),ierr);
    Epetra_Operator_Apply(tmp_oB,X,YB1);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB1,"A Multiply and B Multiply", verbose),ierr);
    Epetra_RowMatrix_Multiply(B, transA, X, YB2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB2,"A Multiply and B Apply", verbose), ierr);

    Epetra_MultiVector_Destroy(&X);
    Epetra_MultiVector_Destroy(&YB2);
    Epetra_MultiVector_Destroy(&YB1);
    Epetra_MultiVector_Destroy(&YA2);
    Epetra_MultiVector_Destroy(&YA1);
    
  }
  {// transpose case
    CT_Epetra_MultiVector_ID_t X = Epetra_MultiVector_Create(bormA, NumVectors, TRUE);
    CT_Epetra_MultiVector_ID_t YA1 = Epetra_MultiVector_Create(bodmA, NumVectors, TRUE);
    CT_Epetra_MultiVector_ID_t YA2 = Epetra_MultiVector_Duplicate(YA1);
    CT_Epetra_MultiVector_ID_t YB1 = Epetra_MultiVector_Duplicate(YA1);
    CT_Epetra_MultiVector_ID_t YB2 = Epetra_MultiVector_Duplicate(YA1);
    Epetra_MultiVector_Random(X);
    
    boolean transA = TRUE;
    Epetra_Operator_SetUseTranspose(tmp_oA, transA);
    Epetra_Operator_SetUseTranspose(tmp_oB, transA);
    Epetra_Operator_Apply(tmp_oA,X,YA1);
    Epetra_RowMatrix_Multiply(A, transA, X, YA2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YA2, "A Multiply and A Apply (transpose)", verbose),ierr);
    Epetra_Operator_Apply(tmp_oB,X,YB1);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB1, "A Multiply and B Multiply (transpose)", verbose),ierr);
    Epetra_RowMatrix_Multiply(B, transA, X, YB2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB2, "A Multiply and B Apply (transpose)", verbose),ierr);
    
    Epetra_MultiVector_Destroy(&X);
    Epetra_MultiVector_Destroy(&YB2);
    Epetra_MultiVector_Destroy(&YB1);
    Epetra_MultiVector_Destroy(&YA2);
    Epetra_MultiVector_Destroy(&YA1);
    
  }

  CT_Epetra_Vector_ID_Flex_t diagA;
  diagA.Epetra_Vector = Epetra_Vector_Create(brrmA, TRUE);
  EPETRA_TEST_ERR(Epetra_RowMatrix_ExtractDiagonalCopy(A, diagA.Epetra_Vector),ierr);
  CT_Epetra_Vector_ID_Flex_t diagB;
  diagB.Epetra_Vector = Epetra_Vector_Create(brrmB, TRUE);
  EPETRA_TEST_ERR(Epetra_RowMatrix_ExtractDiagonalCopy(B, diagB.Epetra_Vector),ierr);
  EPETRA_TEST_ERR(checkMultiVectors(diagA.Epetra_MultiVector,diagB.Epetra_MultiVector,
      "ExtractDiagonalCopy", verbose),ierr);

  CT_Epetra_Vector_ID_Flex_t rowA;
  rowA.Epetra_Vector = Epetra_Vector_Create(brrmA, TRUE);
  EPETRA_TEST_ERR(Epetra_RowMatrix_InvRowSums(A, rowA.Epetra_Vector),ierr);
  CT_Epetra_Vector_ID_Flex_t rowB;
  rowB.Epetra_Vector = Epetra_Vector_Create(brrmB, TRUE);
  EPETRA_TEST_ERR(Epetra_RowMatrix_InvRowSums(B, rowB.Epetra_Vector),ierr)
  EPETRA_TEST_ERR(checkMultiVectors(rowA.Epetra_MultiVector,rowB.Epetra_MultiVector,
      "InvRowSums", verbose),ierr);

  CT_Epetra_Vector_ID_Flex_t colA;
  colA.Epetra_Vector = Epetra_Vector_Create(brcmA, TRUE);
  EPETRA_TEST_ERR(Epetra_RowMatrix_InvColSums(A, colA.Epetra_Vector),ierr);
  CT_Epetra_Vector_ID_Flex_t colB;
  colB.Epetra_Vector = Epetra_Vector_Create(brcmB, TRUE);
  EPETRA_TEST_ERR(Epetra_RowMatrix_InvColSums(B, colB.Epetra_Vector),ierr);
  EPETRA_TEST_ERR(checkMultiVectors(colA.Epetra_MultiVector,colB.Epetra_MultiVector,
      "InvColSums", verbose),ierr);

  EPETRA_TEST_ERR(checkValues(Epetra_RowMatrix_NormInf(A), Epetra_RowMatrix_NormInf(B), "NormInf before scaling", verbose), ierr);
  EPETRA_TEST_ERR(checkValues(Epetra_RowMatrix_NormOne(A), Epetra_RowMatrix_NormOne(B), "NormOne before scaling", verbose),ierr);

  EPETRA_TEST_ERR(Epetra_RowMatrix_RightScale(A, colA.Epetra_Vector),ierr);  
  EPETRA_TEST_ERR(Epetra_RowMatrix_RightScale(B, colB.Epetra_Vector),ierr);


  EPETRA_TEST_ERR(Epetra_RowMatrix_LeftScale(A, rowA.Epetra_Vector),ierr);
  EPETRA_TEST_ERR(Epetra_RowMatrix_LeftScale(B, rowB.Epetra_Vector),ierr);

  Epetra_Vector_Destroy(&colB.Epetra_Vector);
  Epetra_Vector_Destroy(&colA.Epetra_Vector);
  Epetra_Vector_Destroy(&rowB.Epetra_Vector);
  Epetra_Vector_Destroy(&rowA.Epetra_Vector);
  Epetra_Vector_Destroy(&diagB.Epetra_Vector);
  Epetra_Vector_Destroy(&diagA.Epetra_Vector);

  Epetra_BlockMap_Destroy(&brrmA);
  Epetra_BlockMap_Destroy(&brrmB);
  Epetra_BlockMap_Destroy(&brcmA);
  Epetra_BlockMap_Destroy(&brcmB);
  Epetra_BlockMap_Destroy(&bormA);
  Epetra_BlockMap_Destroy(&bormB);
  Epetra_BlockMap_Destroy(&bodmA);
  Epetra_BlockMap_Destroy(&bodmB);

  EPETRA_TEST_ERR(checkValues(Epetra_RowMatrix_NormInf(A), Epetra_RowMatrix_NormInf(B), "NormInf after scaling", verbose), ierr);
  EPETRA_TEST_ERR(checkValues(Epetra_RowMatrix_NormOne(A), Epetra_RowMatrix_NormOne(B), "NormOne after scaling", verbose),ierr);

  vector<double> valuesA(Epetra_RowMatrix_MaxNumEntries(A));
  vector<int> indicesA(Epetra_RowMatrix_MaxNumEntries(A));  
  vector<double> valuesB(Epetra_RowMatrix_MaxNumEntries(B));
  vector<int> indicesB(Epetra_RowMatrix_MaxNumEntries(B));
  return(0);
  for (int i=0; i<Epetra_RowMatrix_NumMyRows(A); i++) {
    int nA, nB;
    EPETRA_TEST_ERR(Epetra_RowMatrix_ExtractMyRowCopy(A, i, Epetra_RowMatrix_MaxNumEntries(A), &nA, &valuesA[0], &indicesA[0]),ierr); 
    EPETRA_TEST_ERR(Epetra_RowMatrix_ExtractMyRowCopy(B, i, Epetra_RowMatrix_MaxNumEntries(B), &nB, &valuesB[0], &indicesB[0]),ierr);
    EPETRA_TEST_ERR(!nA==nB,ierr);
    for (int j=0; j<nA; j++) {
      double curVal = valuesA[j];
      int curIndex = indicesA[j];
      bool notfound = true;
      int jj = 0;
      while (notfound && jj< nB) {
	if (!checkValues(curVal, valuesB[jj])) notfound = false;
	jj++;
      }
      EPETRA_TEST_ERR(notfound, ierr);
      vector<int>::iterator p = find(indicesB.begin(),indicesB.end(),curIndex);  // find curIndex in indicesB
      EPETRA_TEST_ERR(p==indicesB.end(), ierr);
    }

  }
  if (verbose) cout << "RowMatrix Methods check OK" << endl;

  return (ierr);
}
