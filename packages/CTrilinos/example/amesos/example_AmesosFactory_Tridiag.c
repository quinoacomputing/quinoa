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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "CTrilinos_enums.h"
#include "CTrilinos_flex_enums.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "CEpetra_MpiComm.h"
#else
#include "CEpetra_SerialComm.h"
#endif
#include "CEpetra_Comm.h"
#include "CAmesos.h"
#include "CAmesos_BaseSolver.h"
#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_CrsMatrix.h"
#include "CEpetra_Vector.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_Object.h"
#include "CEpetra_LinearProblem.h"
#include "CTeuchos_ParameterList.h"

/*! @file example_AmesosFactory_Tridiag,c
 * This is an example of how to use the CTrilinos interface to Amesos.
 * This example follows from the Amesos example by the same name.
 */

/*
 * ====================
 * M A I N  D R I V E R
 * ====================
 *
 * This example will:
 * 1.- Create an tridiagonal matrix;
 * 2.- Call SymbolicFactorization();
 * 3.- Change the numerical values of the matrix;
 * 4.- Call NumericFactorization();
 * 5.- Set the entries of the RHS;
 * 6.- Call Solve().
 *
 * This example is intended to show the required data
 * for each phase. Phase (2) requires the matrix structure only. 
 * Phase (4) requires the matrix structure (supposed unchanged 
 * from phase (2)) and the matrix data. Phase (6) requires the 
 * RHS and solution vector.
 *
 * This example can be run with any number of processors.
 *
 * Author: Marzio Sala, SNL 9214
 * Last modified: Apr-05.
 */

int main(int argc, char *argv[]) 
{
#ifdef HAVE_CTRILINOS_AMESOS

  int NumGlobalElements, NumMyElements, NumEntries, Indices[3], i, check;
  double Values[3], residual;
  double sfact_time, nfact_time, solve_time, mtx_conv_time, mtx_redist_time, vec_redist_time;
  int *MyGlobalElements = NULL;

  CT_Epetra_Map_ID_Flex_t Map;
  CT_Epetra_CrsMatrix_ID_Flex_t A;
  CT_Epetra_LinearProblem_ID_t Problem;
  CT_Amesos_ID_t Factory;
  CT_Amesos_BaseSolver_ID_t Solver;
  CT_Epetra_Vector_ID_Flex_t b, x, Ax;
  CT_Teuchos_ParameterList_ID_t TimingsList;
  
#ifdef HAVE_MPI
  CT_Epetra_MpiComm_ID_Flex_t Comm;
  MPI_Init(&argc, &argv);
  Comm.Epetra_MpiComm = Epetra_MpiComm_Create(MPI_COMM_WORLD);
#else
  CT_Epetra_SerialComm_ID_Flex_t Comm;
  Comm.Epetra_SerialComm = Epetra_SerialComm_Create();
#endif

  NumGlobalElements = 100; /* global dimension of the problem. */

  /*
   * =======================================================
   * B E G I N N I N G   O F   M A T R I X   C R E A T I O N
   * =======================================================
   */

  /*
   * Construct a Map that puts approximatively the same number of 
   * equations on each processor. `0' is the index base (that is,
   * numbering starts from 0.
   */
  Map.Epetra_Map = Epetra_Map_Create(NumGlobalElements, 0, Comm.Epetra_Comm);

  /* Create an empty EpetraCrsMatrix */
  A.Epetra_CrsMatrix = Epetra_CrsMatrix_Create(CT_Epetra_DataAccess_E_Copy, Map.Epetra_Map, 0, FALSE);

  /* Create the structure of the matrix (tridiagonal) */
  NumMyElements = Epetra_BlockMap_NumMyElements(Map.Epetra_BlockMap);

  /* Add  rows one-at-a-time
   * Need some vectors to help
   */

  /* Right now, we put zeros only in the matrix. */
  Values[0] = 0.0;
  Values[1] = 0.0;
  Values[2] = 0.0;
  /* global ID's of local ID's */
  MyGlobalElements = Epetra_BlockMap_MyGlobalElements(Map.Epetra_BlockMap);

  /* At this point we simply set the nonzero structure of A.
   * Actual values will be inserted later (now all zeros)
   */
  for (i = 0; i < NumMyElements; i++) 
  {
    if (MyGlobalElements[i] == 0) 
    {
      Indices[0] = 0; 
      Indices[1] = 1;
      NumEntries = 2;
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1) 
    {
      Indices[0] = NumGlobalElements-1;
      Indices[1] = NumGlobalElements-2;
      NumEntries = 2;
    }
    else 
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i];
      Indices[2] = MyGlobalElements[i]+1;
      NumEntries = 3;
    }

    check = Epetra_CrsMatrix_InsertGlobalValues(A.Epetra_CrsMatrix, MyGlobalElements[i],
        NumEntries, Values, Indices);
    assert(check == 0);
  }

  /* Finish up. */
  Epetra_CrsMatrix_FillComplete(A.Epetra_CrsMatrix, TRUE);

  /*
   * ===========================================
   * E N D   O F   M A T R I X   C R E A T I O N
   * ===========================================
   */

  /*
   * Now the matrix STRUCTURE is set. We cannot add
   * new nonzero elements, but we can still change the
   * numerical values of all inserted elements (as we will
   * do later).
   */

  /*
   * =====================================================
   * B E G I N N I N G   O F  T H E   AM E S O S   P A R T
   * =====================================================
   */

  /*
   * For comments on the commands in this section, please
   * see file example_AmesosFactory.cpp.
   */

  Problem = Epetra_LinearProblem_Create();  

  Epetra_LinearProblem_SetOperator_Matrix(Problem, A.Epetra_RowMatrix);

  /* Initializes Amesos solver. Here we solve with Amesos_Klu. */

  Factory = Amesos_Create();
  Solver = Amesos_CreateSolver(Factory, "Amesos_Klu", Problem);

  /* Factory.Create() returns 0 if the requested solver
   * is not available
   */

/* 
  if (Solver == 0) {
    std::cerr << "Selected solver is not available" << std::endl;
    // return ok not to break the test harness
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }
*/

  /* At this point we can perform the numeric factorization.
   * Note that the matrix contains 0's only.
   */

  Amesos_BaseSolver_SymbolicFactorization(Solver);

  /* Now, we repopulate the matrix with entries corresponding
   * to a 1D Laplacian. LHS and RHS are still untouched.
   */

  for (i = 0; i < NumMyElements; i++) 
  {
    if (MyGlobalElements[i] == 0) 
    {
      Indices[0] = 0;   
      Indices[1] = 1;
      Values[0]  = 2.0; 
      Values[1]  = -1.0;
      NumEntries = 2;
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1) 
    {
      Indices[0] = NumGlobalElements - 1;
      Indices[1] = NumGlobalElements - 2;
      Values[0]  = 2.0; 
      Values[1]  = -1.0;
      NumEntries = 2;
    }
    else 
    {
      Indices[0] = MyGlobalElements[i] - 1;
      Indices[1] = MyGlobalElements[i];
      Indices[2] = MyGlobalElements[i] + 1;
      Values[0] = -1.0; 
      Values[1] = 2.0; 
      Values[2] = -1.0;
      NumEntries = 3;
    }

    check = Epetra_CrsMatrix_ReplaceGlobalValues(A.Epetra_CrsMatrix, MyGlobalElements[i],
        NumEntries, Values, Indices);
    assert(check == 0);
  }

  /* ... and we can compute the numeric factorization. */
  Amesos_BaseSolver_NumericFactorization(Solver);

  /* Finally, we set up the LHS and the RHS vector (Random(). */
  b.Epetra_Vector = Epetra_Vector_Create(Map.Epetra_BlockMap, TRUE);
  Epetra_MultiVector_Random(b.Epetra_MultiVector);

  x.Epetra_Vector = Epetra_Vector_Create(Map.Epetra_BlockMap, TRUE);
  Epetra_MultiVector_PutScalar(x.Epetra_MultiVector, 0.0);

  Epetra_LinearProblem_SetLHS(Problem, x.Epetra_MultiVector);
  Epetra_LinearProblem_SetRHS(Problem, b.Epetra_MultiVector);
  
  Amesos_BaseSolver_Solve(Solver);

  /* Print out the timing information and get it from the solver */
  Amesos_BaseSolver_PrintTiming(Solver);
  
  TimingsList = Teuchos_ParameterList_Create();
  Amesos_BaseSolver_GetTiming(Solver, TimingsList);
  
  /* You can find out how much time was spent in ...
   * sfact_time, nfact_time, solve_time,
   * mtx_conv_time, mtx_redist_time, vec_redist_time */

  /* 1) The symbolic factorization
   *    (parameter doesn't always exist)
   */
  sfact_time = Teuchos_ParameterList_get_def_double( TimingsList, "Total symbolic factorization time", 0.0 );

  /* 2) The numeric factorization
   *    (always exists if NumericFactorization() is called)
   */
  nfact_time = Teuchos_ParameterList_get_double( TimingsList, "Total numeric factorization time" );

  /* 3) Solving the linear system
   *    (always exists if Solve() is called)
   */
  solve_time = Teuchos_ParameterList_get_double( TimingsList, "Total solve time" );

  /* 4) Converting the matrix to the accepted format for the solver
   *    (always exists if SymbolicFactorization() is called)
   */
  mtx_conv_time = Teuchos_ParameterList_get_double( TimingsList, "Total solve time" );

  /* 5) Redistributing the matrix for each solve to the accepted format for the solver */
  mtx_redist_time = Teuchos_ParameterList_get_def_double( TimingsList, "Total matrix redistribution time", 0.0 );

  /* 6) Redistributing the vector for each solve to the accepted format for the solver */
  vec_redist_time = Teuchos_ParameterList_get_def_double( TimingsList, "Total vector redistribution time", 0.0 );

  /*
   * ===========================================
   * E N D   O F   T H E   A M E S O S   P A R T
   * ===========================================
   */

  /*
   * ==================
   * compute ||Ax - b||
   * ==================
   */

  Ax.Epetra_Vector = Epetra_Vector_Create(Map.Epetra_BlockMap, TRUE);
  Epetra_CrsMatrix_Multiply_Vector(A.Epetra_CrsMatrix, FALSE, x.Epetra_Vector, Ax.Epetra_Vector);
  Epetra_MultiVector_Update_WithA(Ax.Epetra_MultiVector, 1.0, b.Epetra_MultiVector, -1.0);
  Epetra_MultiVector_Norm2(Ax.Epetra_MultiVector, &residual);

  if (!Epetra_Comm_MyPID(Comm.Epetra_Comm))
    printf("After AMESOS solution, ||b-Ax||_2 = %g\n", residual);

  /* delete Solver. Do this before MPI_Finalize()
   * as MPI calls can occur in the destructor.
   */
  Amesos_BaseSolver_Destroy(&Solver);
    
  if (residual > 1e-5)
    return(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

#endif /* HAVE_CTRILINOS_AMESOS */

  return(EXIT_SUCCESS);

} /* end of main() */

