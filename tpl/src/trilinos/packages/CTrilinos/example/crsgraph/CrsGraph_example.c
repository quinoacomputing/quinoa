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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifdef HAVE_MALLOC_H
#ifndef __APPLE__
#include <malloc.h>
#else
#include <sys/malloc.h>
#endif
#endif

#include "CTrilinos_enums.h"
#include "CTrilinos_flex_enums.h"
#include "CEpetra_Comm.h"
#ifdef HAVE_MPI
#include "CEpetra_MpiComm.h"
#include <mpi.h>
#else
#include "CEpetra_SerialComm.h"
#endif
#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_CrsGraph.h"

/*! @file CrsGraph_example.c
 * This is an example of how to use the CTrilinos interface to Epetra.
 * This example follows from the Epetra example my_example.
 */

/* This function is to be used when first identifying an error. */
#define CTRILINOS_TEST_ERR(a,b) { { \
  if ((int)(a) != 0) {\
    fprintf(stderr, "Non zero error code %d, file: %s, line: %d\n", \
            (int)(a), __FILE__, __LINE__); \
    if ((int)(a) < 0) (b)+=1; \
  }\
  }\
}

/* Prototype */
int check(CT_Epetra_CrsGraph_ID_t A, int NumMyRows1, int NumGlobalRows1, int NumMyNonzeros1,
	  int NumGlobalNonzeros1, int* MyGlobalElements, boolean verbose);

int main(int argc, char *argv[])
{
  int ierr, forierr, tempierr, One, i;
  int MyPID, NumProc, NumMyEquations, NumGlobalEquations, NumEntries, NumMyNonzeros;
  int NumMyElements1, NumMyEquations1, NumGlobalEquations1, NumEntries1;
  boolean debug, verbose, verbose1;
  int Indices[2], Indices1[2];
  int *MyGlobalElements, *NumNz, *MyGlobalElements1, *NumNz1;

  CT_Epetra_Map_ID_Flex_t Map, Map1;
  CT_Epetra_CrsGraph_ID_t A, AA, A1, B;

#ifdef HAVE_MPI

  /* Initialize MPI */
  CT_Epetra_MpiComm_ID_Flex_t Comm;
  MPI_Init(&argc,&argv);
  Comm.Epetra_MpiComm = Epetra_MpiComm_Create( MPI_COMM_WORLD );

#else

  CT_Epetra_SerialComm_ID_Flex_t Comm;
  Comm.Epetra_SerialComm = Epetra_SerialComm_Create();

#endif

  ierr = 0;
  forierr = 0;
  tempierr = 0;
  One = 1;

  debug = TRUE;
  verbose = FALSE;
  verbose1 = FALSE;

  /* Check if we should print results to standard out */
  if(argc > 1) {
    if(argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = TRUE;
    }
  }

  MyPID = Epetra_Comm_MyPID(Comm.Epetra_Comm);
  NumProc = Epetra_Comm_NumProc(Comm.Epetra_Comm);
  if(verbose) printf("Processor %d of %d is alive.\n", MyPID, NumProc);

  verbose1 = verbose;

  /* Redefine verbose to only print on PE 0 */
  if(verbose && MyPID != 0) verbose = FALSE;

  NumMyEquations = 5;
  NumGlobalEquations = NumMyEquations*NumProc+(NumProc < 3 ? NumProc : 3);
  if(MyPID < 3) 
    NumMyEquations++;

  /* Construct a Map that puts approximately the same Number of equations
   * on each processor */

  Map.Epetra_Map = Epetra_Map_Create_Linear(NumGlobalEquations, NumMyEquations, 0, Comm.Epetra_Comm);
  
  /* Get update list and number of local equations from newly created Map */
  MyGlobalElements = (int*)malloc(NumMyEquations*sizeof(int));
  if (MyGlobalElements == NULL) {
    fprintf(stderr,"Couldn't malloc for MyGlobalElements\n");
    printf( "\nEnd Result: TEST FAILED\n" );
#ifdef HAVE_MPI
    MPI_Finalize();
#endif /* HAVE_MPI */
    return 1;
  }
  Epetra_BlockMap_MyGlobalElements_Fill(Map.Epetra_BlockMap, MyGlobalElements);

  /* Create an integer vector NumNz that is used to build the Petra Matrix.
   * NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
   *on this processor */

  NumNz = (int*)malloc(NumMyEquations*sizeof(int));
  if (NumNz == NULL) {
    fprintf(stderr,"Couldn't malloc for NumNz\n");
    printf( "\nEnd Result: TEST FAILED\n" );
#ifdef HAVE_MPI
    MPI_Finalize();
#endif /* HAVE_MPI */
    return 1;
  }

  /* We are building a tridiagonal matrix where each row has (-1 2 -1)
   * So we need 2 off-diagonal terms (except for the first and last equation) */

  for(i=0; i<NumMyEquations; i++)
    if(MyGlobalElements[i]==0 || MyGlobalElements[i] == NumGlobalEquations-1)
      NumNz[i] = 1;
    else
      NumNz[i] = 2;

  /* Create a Epetra_CrsGraph */

  A = Epetra_CrsGraph_Create_VarPerRow(CT_Epetra_DataAccess_E_Copy,
         Map.Epetra_BlockMap, NumNz, FALSE);
  CTRILINOS_TEST_ERR(Epetra_CrsGraph_IndicesAreGlobal(A),ierr);
  CTRILINOS_TEST_ERR(Epetra_CrsGraph_IndicesAreLocal(A),ierr);
  
  /* Add  rows one-at-a-time
   * Need some vectors to help
   * Off diagonal Values will always be -1 */

  forierr = 0;
  for(i = 0; i < NumMyEquations; i++) {
    if(MyGlobalElements[i] == 0) {
      Indices[0] = 1;
      NumEntries = 1;
    } else if(MyGlobalElements[i] == NumGlobalEquations-1) {
      Indices[0] = NumGlobalEquations-2;
      NumEntries = 1;
    } else {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      NumEntries = 2;
    }
    forierr += !(Epetra_CrsGraph_InsertGlobalIndices(A, MyGlobalElements[i],
                 NumEntries, Indices)==0);
    forierr += !(Epetra_CrsGraph_InsertGlobalIndices(A, MyGlobalElements[i],
                 1, MyGlobalElements+i)>0); /* Put in the diagonal entry (should cause realloc) */
  }
  CTRILINOS_TEST_ERR(forierr,ierr);
  
  /* Finish up */
  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_IndicesAreGlobal(A)),ierr);
  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_FillComplete(A)==0),ierr);
  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_IndicesAreLocal(A)),ierr);
  CTRILINOS_TEST_ERR(Epetra_CrsGraph_StorageOptimized(A),ierr);

  Epetra_CrsGraph_OptimizeStorage(A);

  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_StorageOptimized(A)),ierr);
  CTRILINOS_TEST_ERR(Epetra_CrsGraph_UpperTriangular(A),ierr);
  CTRILINOS_TEST_ERR(Epetra_CrsGraph_LowerTriangular(A),ierr);

  if(verbose) printf("\n*****Testing variable entry constructor\n\n");

  NumMyNonzeros = 3 * NumMyEquations;
  if(Epetra_CrsGraph_LRID(A, 0) >= 0) 
    NumMyNonzeros--; /* If I own first global row, then there is one less nonzero */
  if(Epetra_CrsGraph_LRID(A, NumGlobalEquations-1) >= 0) 
    NumMyNonzeros--; /* If I own last global row, then there is one less nonzero */

  CTRILINOS_TEST_ERR(check(A, NumMyEquations, NumGlobalEquations,
      NumMyNonzeros, 3*NumGlobalEquations-2, MyGlobalElements, verbose),ierr);
  forierr = 0;
  for(i = 0; i < NumMyEquations; i++) 
    forierr += !(Epetra_CrsGraph_NumGlobalIndices(A, MyGlobalElements[i])==NumNz[i]+1);
  CTRILINOS_TEST_ERR(forierr,ierr);
  for(i = 0; i < NumMyEquations; i++) 
    forierr += !(Epetra_CrsGraph_NumMyIndices(A, i)==NumNz[i]+1);
  CTRILINOS_TEST_ERR(forierr,ierr);

  if(verbose) printf("NumIndices function check OK\n");

  Epetra_CrsGraph_Destroy(&A);

  if(debug) Epetra_Comm_Barrier(Comm.Epetra_Comm);

  if(verbose) printf("\n*****Testing constant entry constructor\n\n");

  AA = Epetra_CrsGraph_Create(CT_Epetra_DataAccess_E_Copy, Map.Epetra_BlockMap, 5, FALSE);
  
  if(debug) Epetra_Comm_Barrier(Comm.Epetra_Comm);

  for(i = 0; i < NumMyEquations; i++) 
    Epetra_CrsGraph_InsertGlobalIndices(AA, MyGlobalElements[i], 1,
         MyGlobalElements+i);

  /* Note:  All processors will call the following Insert routines, but only the processor
   *        that owns it will actually do anything */

  if(Epetra_CrsGraph_MyGlobalRow(AA, 0)) {
    CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_InsertGlobalIndices(AA, 0, 0, &One)==0),ierr);
  } else {
   CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_InsertGlobalIndices(AA, 0, 1, &One)==-2),ierr);
  }
  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_FillComplete(AA)==0),ierr);
  CTRILINOS_TEST_ERR(Epetra_CrsGraph_StorageOptimized(AA),ierr);
  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_UpperTriangular(AA)),ierr);
  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_LowerTriangular(AA)),ierr);

  if(debug) Epetra_Comm_Barrier(Comm.Epetra_Comm);
  CTRILINOS_TEST_ERR(check(AA, NumMyEquations, NumGlobalEquations,
      NumMyEquations, NumGlobalEquations, MyGlobalElements, verbose),ierr);

  if(debug) Epetra_Comm_Barrier(Comm.Epetra_Comm);

  forierr = 0;
  for(i = 0; i < NumMyEquations; i++) 
    forierr += !(Epetra_CrsGraph_NumGlobalIndices(AA, MyGlobalElements[i])==1);
  CTRILINOS_TEST_ERR(forierr,ierr);

  if(verbose) printf("NumIndices function check OK\n");

  if(debug) Epetra_Comm_Barrier(Comm.Epetra_Comm);

  if(verbose) printf("\n*****Testing copy constructor\n\n");

  B = Epetra_CrsGraph_Duplicate(AA);

  Epetra_CrsGraph_Destroy(&AA);

  CTRILINOS_TEST_ERR(check(B, NumMyEquations, NumGlobalEquations, NumMyEquations,
      NumGlobalEquations, MyGlobalElements, verbose),ierr);

  forierr = 0;
  for(i = 0; i < NumMyEquations; i++) 
    forierr += !(Epetra_CrsGraph_NumGlobalIndices(B, MyGlobalElements[i])==1);
  CTRILINOS_TEST_ERR(forierr,ierr);

  if(verbose) printf("NumIndices function check OK\n");

  if(debug) Epetra_Comm_Barrier(Comm.Epetra_Comm);

  if(verbose) printf("\n*****Testing post construction modifications\n\n");

  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_InsertGlobalIndices(B, 0, 1, &One)==-2),ierr);

  /* Release all objects */
  free(MyGlobalElements);
  free(NumNz);

  Epetra_BlockMap_Destroy(&Map.Epetra_BlockMap);
  Epetra_CrsGraph_Destroy(&B);


  if (verbose1) {
    /* Test ostream << operator (if verbose1)
     * Construct a Map that puts 2 equations on each PE */
    
    NumMyElements1 = 4;
    NumMyEquations1 = NumMyElements1;
    NumGlobalEquations1 = NumMyEquations1*NumProc;

    Map1.Epetra_Map = Epetra_Map_Create_Linear(NumGlobalEquations1, NumMyEquations1, 1, Comm.Epetra_Comm);
    
    /* Get update list and number of local equations from newly created Map */
    MyGlobalElements1 = (int*)malloc(NumMyEquations1*sizeof(int));
    if (MyGlobalElements1 == NULL) {
      fprintf(stderr,"Couldn't malloc for MyGlobalElements1\n");
      printf( "\nEnd Result: TEST FAILED\n" );
#ifdef HAVE_MPI
      MPI_Finalize();
#endif /* HAVE_MPI */
      return 1;
    }
    Epetra_BlockMap_MyGlobalElements_Fill(Map1.Epetra_BlockMap, MyGlobalElements1);
    
    /* Create an integer vector NumNz that is used to build the Petra Matrix.
     * NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
     * on this processor */
    
    NumNz1 = (int*)malloc(NumMyEquations1*sizeof(int));
    if (NumNz1 == NULL) {
      fprintf(stderr,"Couldn't malloc for NumNz1\n");
      printf( "\nEnd Result: TEST FAILED\n" );
#ifdef HAVE_MPI
      MPI_Finalize();
#endif /* HAVE_MPI */
      return 1;
    }
    
    /* We are building a tridiagonal matrix where each row has (-1 2 -1)
     * So we need 2 off-diagonal terms (except for the first and last equation) */
    
    for(i = 0; i < NumMyEquations1; i++)
      if(MyGlobalElements1[i]==1 || MyGlobalElements1[i] == NumGlobalEquations1)
	NumNz1[i] = 1;
      else
	NumNz1[i] = 2;
    
    /* Create a Epetra_Graph using 1-based arithmetic */
    
    A1 = Epetra_CrsGraph_Create_VarPerRow(CT_Epetra_DataAccess_E_Copy, Map1.Epetra_BlockMap,
              NumNz1, FALSE);
    
    /* Add  rows one-at-a-time
     * Need some vectors to help
     * Off diagonal Values will always be -1 */
    
    forierr = 0;
    for(i = 0; i < NumMyEquations1; i++) {
      if(MyGlobalElements1[i]==1) {
	Indices1[0] = 2;
	NumEntries1 = 1;
      } else if(MyGlobalElements1[i] == NumGlobalEquations1) {
	Indices1[0] = NumGlobalEquations1-1;
	NumEntries1 = 1;
      } else {
	Indices1[0] = MyGlobalElements1[i]-1;
	Indices1[1] = MyGlobalElements1[i]+1;
	NumEntries1 = 2;
      }
      forierr += !(Epetra_CrsGraph_InsertGlobalIndices(A1, MyGlobalElements1[i],
          NumEntries1, Indices1)==0);
      forierr += !(Epetra_CrsGraph_InsertGlobalIndices(A1, MyGlobalElements1[i],
          1, MyGlobalElements1+i)>0); /* Put in the diagonal entry */
    }
    CTRILINOS_TEST_ERR(forierr,ierr);
		
    /* Finish up */
    CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_FillComplete(A1)==0),ierr);
    
    /* if(verbose)
     *   printf("Print out tridiagonal matrix, each part on each processor. Index base is one.\n\n"); */
    /* cout << A1 << endl; */
    
    /* Release all objects */
    free(MyGlobalElements1);
    free(NumNz1);

    Epetra_CrsGraph_Destroy(&A1);
    Epetra_BlockMap_Destroy(&Map1.Epetra_BlockMap);
  }

  forierr = 0;
  if (ierr < 0) ierr *= -1;
  Epetra_Comm_MaxAll_Int(Comm.Epetra_Comm, &ierr, &forierr, 1);

  if (forierr == 0)
    printf( "\nEnd Result: TEST PASSED\n" );
  else
    printf( "\nEnd Result: TEST FAILED\n" );
  
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
/* end main */

  return((forierr == 0) ? 0 : 1);
}

/* ============================================================================== */
int check(CT_Epetra_CrsGraph_ID_t A, int NumMyRows1, int NumGlobalRows1, int NumMyNonzeros1,
	  int NumGlobalNonzeros1, int* MyGlobalElements, boolean verbose)
{  
  int ierr, forierr, i, j, Row, maxmygid, minmygid;
  int NumGlobalIndices, NumMyIndices, MaxNumIndices;
  int NumMyRows, NumMyNonzeros, NumGlobalRows, NumGlobalNonzeros;
  int *MyViewIndices, *MyCopyIndices, *GlobalCopyIndices;

  CT_Epetra_BlockMap_ID_t rowMap;

  (void)MyGlobalElements;

  ierr = 0;
  forierr = 0;

  MaxNumIndices = Epetra_CrsGraph_MaxNumIndices(A);
  MyCopyIndices = (int*)malloc(MaxNumIndices*sizeof(int));
  assert(MyCopyIndices != NULL);
  GlobalCopyIndices = (int*)malloc(MaxNumIndices*sizeof(int));
  assert(GlobalCopyIndices != NULL);

  /* Test query functions */

  NumMyRows = Epetra_CrsGraph_NumMyRows(A);
  if(verbose) printf("Number of local Rows = %d\n", NumMyRows);

  CTRILINOS_TEST_ERR(!(NumMyRows==NumMyRows1),ierr);

  NumMyNonzeros = Epetra_CrsGraph_NumMyNonzeros(A);
  if(verbose) printf("Number of local Nonzero entries = %d\n", NumMyNonzeros);

  CTRILINOS_TEST_ERR(!(NumMyNonzeros==NumMyNonzeros1),ierr);

  NumGlobalRows = Epetra_CrsGraph_NumGlobalRows(A);
  if(verbose) printf("Number of global Rows = %d\n", NumGlobalRows);

  CTRILINOS_TEST_ERR(!(NumGlobalRows==NumGlobalRows1),ierr);

  NumGlobalNonzeros = Epetra_CrsGraph_NumGlobalNonzeros(A);
  if(verbose) printf("Number of global Nonzero entries = %d\n", NumGlobalNonzeros);

  CTRILINOS_TEST_ERR(!(NumGlobalNonzeros==NumGlobalNonzeros1),ierr);

  /* GlobalRowView should be illegal (since we have local indices) */

  rowMap = Epetra_CrsGraph_RowMap(A);
  maxmygid = Epetra_BlockMap_MaxMyGID(rowMap);
  minmygid = Epetra_BlockMap_MinMyGID(rowMap);

  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_ExtractGlobalRowView(A, maxmygid, &NumGlobalIndices, &GlobalCopyIndices)==-2),ierr);

  /* Other binary tests */

  CTRILINOS_TEST_ERR(Epetra_CrsGraph_NoDiagonal(A),ierr);
  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_Filled(A)),ierr);
  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_MyGRID(A, maxmygid)),ierr);
  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_MyGRID(A, minmygid)),ierr);
  CTRILINOS_TEST_ERR(Epetra_CrsGraph_MyGRID(A, 1+maxmygid),ierr);
  CTRILINOS_TEST_ERR(Epetra_CrsGraph_MyGRID(A, -1+minmygid),ierr);
  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_MyLRID(A, 0)),ierr);
  CTRILINOS_TEST_ERR(!(Epetra_CrsGraph_MyLRID(A, NumMyRows-1)),ierr);
  CTRILINOS_TEST_ERR(Epetra_CrsGraph_MyLRID(A, -1),ierr);
  CTRILINOS_TEST_ERR(Epetra_CrsGraph_MyLRID(A, NumMyRows),ierr);
    
  forierr = 0;
  for(i = 0; i < NumMyRows; i++) {
    Row = Epetra_CrsGraph_GRID(A, i);
    Epetra_CrsGraph_ExtractGlobalRowCopy(A, Row, MaxNumIndices, &NumGlobalIndices, GlobalCopyIndices);
    Epetra_CrsGraph_ExtractMyRowView(A, i, &NumMyIndices, &MyViewIndices);
    forierr += !(NumGlobalIndices==NumMyIndices);
    for(j = 1; j < NumMyIndices; j++) CTRILINOS_TEST_ERR(!(MyViewIndices[j-1]<MyViewIndices[j]),ierr);
    for(j = 0; j < NumGlobalIndices; j++) {
      forierr += !(GlobalCopyIndices[j]==Epetra_CrsGraph_GCID(A, MyViewIndices[j]));
      forierr += !(Epetra_CrsGraph_LCID(A, GlobalCopyIndices[j])==MyViewIndices[j]);
    }
  }
  CTRILINOS_TEST_ERR(forierr,ierr);
  forierr = 0;
  for(i = 0; i < NumMyRows; i++) {
    Row = Epetra_CrsGraph_GRID(A, i);
    Epetra_CrsGraph_ExtractGlobalRowCopy(A, Row, MaxNumIndices, &NumGlobalIndices, GlobalCopyIndices);
    Epetra_CrsGraph_ExtractMyRowCopy(A, i, MaxNumIndices, &NumMyIndices, MyCopyIndices);
    forierr += !(NumGlobalIndices==NumMyIndices);
    for(j = 1; j < NumMyIndices; j++) 
      CTRILINOS_TEST_ERR(!(MyCopyIndices[j-1]<MyCopyIndices[j]),ierr);
    for(j = 0; j < NumGlobalIndices; j++) {
      forierr += !(GlobalCopyIndices[j]==Epetra_CrsGraph_GCID(A, MyCopyIndices[j]));
      forierr += !(Epetra_CrsGraph_LCID(A, GlobalCopyIndices[j])==MyCopyIndices[j]);
    }

  }
  CTRILINOS_TEST_ERR(forierr,ierr);

  free(MyCopyIndices);
  free(GlobalCopyIndices);

  Epetra_BlockMap_Destroy(&rowMap);

  if(verbose) printf("Rows sorted check OK\n");

  return(ierr);
}

