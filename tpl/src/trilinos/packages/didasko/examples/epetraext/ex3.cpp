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


#include "Didasko_config.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_EPETRAEXT)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_RowMatrixOut.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // creates few simple objects
  int NumMyElements = 5 * Comm.NumProc();
  Epetra_Map Map(-1, NumMyElements, 0, Comm);
  Epetra_Vector X(Map);
  X.Random();
  Epetra_CrsMatrix A(Copy, Map, 0);

  // create a simple diagonal matrix
  for (int i = 0 ; i < NumMyElements ; ++i)
  {
    int j = Map.GID(i);
    double value = 1.0;
    EPETRA_CHK_ERR(A.InsertGlobalValues(j, 1, &value, &j));
  }
  A.FillComplete();

  EpetraExt::BlockMapToMatrixMarketFile("Map.mm", Map, "test map", "This is a test map");
  EpetraExt::VectorToMatrixMarketFile("X.mm", X, "test vector", "This is a test vector");
  EpetraExt::RowMatrixToMatrixMarketFile("A.mm", A, "test matrix", "This is a test matrix");

  // to read the output in MATLAB:
  // 1) download the mmread.m file from the Matrix Market web site
  // 2) use A = mmread('A.mm')

  // one can also read from Matrix Market

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(0);
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
      "--enable-epetraext\n");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

#endif

