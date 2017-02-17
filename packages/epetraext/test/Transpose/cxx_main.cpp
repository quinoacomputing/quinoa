//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

// Epetra_BlockMap Test routine

#include "EpetraExt_Version.h"

#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Time.h"
#include "EpetraExt_Transpose_RowMatrix.h"
#include "Trilinos_Util.h"

int checkResults(Epetra_RowMatrix * A, Epetra_CrsMatrix * transA, 
                 Epetra_Vector * xexact, bool verbose);

int main(int argc, char *argv[])
{
  int i;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // Uncomment to debug in parallel int tmp; if (comm.MyPID()==0) cin >> tmp; comm.Barrier();

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (!verbose) comm.SetTracebackMode(0); // This should shut down any error traceback reporting

  if (verbose) std::cout << comm << std::endl << std::flush;

  if (verbose) verbose = (comm.MyPID()==0);

  if (verbose)
    std::cout << EpetraExt::EpetraExt_Version() << std::endl << std::endl;

  int nx = 128;
  int ny = comm.NumProc()*nx; // Scale y grid with number of processors

  // Create funky stencil to make sure the matrix is non-symmetric (transpose non-trivial):

  // (i-1,j-1) (i-1,j  )
  // (i  ,j-1) (i  ,j  ) (i  ,j+1)
  // (i+1,j-1) (i+1,j  )

  int npoints = 7;

  int xoff[] = {-1,  0,  1, -1,  0,  1,  0};
  int yoff[] = {-1, -1, -1,  0,  0,  0,  1};

  Epetra_Map * map;
  Epetra_CrsMatrix * A;
  Epetra_Vector * x, * b, * xexact;
	
  Trilinos_Util_GenerateCrsProblem(nx, ny, npoints, xoff, yoff, comm, map, A, x, b, xexact);

  if (nx<8)
  {
    std::cout << *A << std::endl;
    std::cout << "X exact = " << std::endl << *xexact << std::endl;
    std::cout << "B       = " << std::endl << *b << std::endl;
  }

  // Construct transposer 
  Epetra_Time timer(comm);

  double start = timer.ElapsedTime();

  //bool IgnoreNonLocalCols = false;
  EpetraExt::RowMatrix_Transpose transposer;

  if (verbose) std::cout << "\nTime to construct transposer  = " << timer.ElapsedTime() - start << std::endl;
  
  Epetra_CrsMatrix & transA = dynamic_cast<Epetra_CrsMatrix&>(transposer(*A));

  start = timer.ElapsedTime();
  if (verbose) std::cout << "\nTime to create transpose matrix  = " << timer.ElapsedTime() - start << std::endl;
 	
  // Now test output of transposer by performing matvecs
  int ierr = 0;
  ierr += checkResults(A, &transA, xexact, verbose);


  // Now change values in original matrix and test update facility of transposer
  // Add 2 to the diagonal of each row
  double Value = 2.0;
  for (i=0; i< A->NumMyRows(); i++)
  A->SumIntoMyValues(i, 1, &Value, &i);

  start = timer.ElapsedTime();
  transposer.fwd();

  if (verbose) std::cout << "\nTime to update transpose matrix  = " << timer.ElapsedTime() - start << std::endl;
 	
  ierr += checkResults(A, &transA, xexact, verbose);

  delete A;
  delete b;
  delete x;
  delete xexact;
  delete map;

  if (verbose) std::cout << std::endl << "Checking transposer for VbrMatrix objects" << std::endl<< std::endl;

  int nsizes = 4;
  int sizes[] = {4, 6, 5, 3};

  Epetra_VbrMatrix * Avbr;
  Epetra_BlockMap * bmap;

  Trilinos_Util_GenerateVbrProblem(nx, ny, npoints, xoff, yoff, nsizes, sizes,
                                   comm, bmap, Avbr, x, b, xexact);

  if (nx<8)
  {
    std::cout << *Avbr << std::endl;
    std::cout << "X exact = " << std::endl << *xexact << std::endl;
    std::cout << "B       = " << std::endl << *b << std::endl;
  }

  start = timer.ElapsedTime();
  EpetraExt::RowMatrix_Transpose transposer1;

  Epetra_CrsMatrix & transA1 = dynamic_cast<Epetra_CrsMatrix&>(transposer1(*Avbr));
  if (verbose) std::cout << "\nTime to create transpose matrix  = " << timer.ElapsedTime() - start << std::endl;
 	
  // Now test output of transposer by performing matvecs
;
  ierr += checkResults(Avbr, &transA1, xexact, verbose);

  // Now change values in original matrix and test update facility of transposer
  // Scale matrix on the left by rowsums

  Epetra_Vector invRowSums(Avbr->RowMap());

  Avbr->InvRowSums(invRowSums);
  Avbr->LeftScale(invRowSums);

  start = timer.ElapsedTime();
  transposer1.fwd();
  if (verbose) std::cout << "\nTime to update transpose matrix  = " << timer.ElapsedTime() - start << std::endl;
 	
  ierr += checkResults(Avbr, &transA1, xexact, verbose);

  delete Avbr;
  delete b;
  delete x;
  delete xexact;
  delete bmap;

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ierr;
}

int checkResults(Epetra_RowMatrix * A, Epetra_CrsMatrix * transA, 
                 Epetra_Vector * xexact, bool verbose) {

  int n = A->NumGlobalRows();

  if (n<100) std::cout << "A transpose = " << std::endl << *transA << std::endl;

  Epetra_Vector x1(View,A->OperatorDomainMap(), &((*xexact)[0]));
  Epetra_Vector b1(A->OperatorRangeMap());

  A->SetUseTranspose(true);

  Epetra_Time timer(A->Comm());
  double start = timer.ElapsedTime();
  A->Apply(x1, b1);
  if (verbose) std::cout << "\nTime to compute b1: matvec with original matrix using transpose flag  = " << timer.ElapsedTime() - start << std::endl;

  if (n<100) std::cout << "b1 = " << std::endl << b1 << std::endl;
  Epetra_Vector x2(View,transA->OperatorRangeMap(), &((*xexact)[0]));
  Epetra_Vector b2(transA->OperatorDomainMap());
  start = timer.ElapsedTime();
  transA->Multiply(false, x2, b2);
  if (verbose) std::cout << "\nTime to compute b2: matvec with transpose matrix                      = " << timer.ElapsedTime() - start << std::endl;

  if (n<100) std::cout << "b1 = " << std::endl << b1 << std::endl;

  double residual;
  Epetra_Vector resid(A->OperatorDomainMap());

  resid.Update(1.0, b1, -1.0, b2, 0.0);
  resid.Norm2(&residual);
  if (verbose) std::cout << "Norm of b1 - b2 = " << residual << std::endl;

  int ierr = 0;

  if (residual > 1.0e-10) ierr++;

  if (ierr!=0 && verbose) std::cerr << "Status: Test failed" << std::endl;
  else if (verbose) std::cerr << "Status: Test passed" << std::endl;

  return(ierr);
}
