
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

// Using IFPACK factorizations as Aztec's preconditioners

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_IFPACK) && defined(HAVE_DIDASKO_AZTECOO)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "AztecOO.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"

// function for fancy output

string toString(const int& x) {
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}

string toString(const double& x) {
  char s[100];
  sprintf(s, "%g", x);
  return string(s);
}

// main driver

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  //const bool verbose = (Comm.MyPID()==0);

  // B E G I N   O F   M A T R I X   C O N S T R U C T I O N

  // matrix downloaded from MatrixMarket
  char FileName[] = "../HBMatrices/fidap005.rua";

  Epetra_Map * readMap; // Pointers because of Trilinos_Util_ReadHb2Epetra
  Epetra_CrsMatrix * readA;
  Epetra_Vector * readx;
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;

  // Call routine to read in HB problem
  Trilinos_Util_ReadHb2Epetra(FileName, Comm, readMap, readA, readx,
      readb, readxexact);

  int NumGlobalElements = readMap->NumGlobalElements();

  // Create uniform distributed map
  Epetra_Map map(NumGlobalElements, 0, Comm);

  // Create Exporter to distribute read-in matrix and vectors

  Epetra_Export exporter(*readMap, map);
  Epetra_CrsMatrix A(Copy, map, 0);
  Epetra_Vector x(map);
  Epetra_Vector b(map);
  Epetra_Vector xexact(map);

  Epetra_Time FillTimer(Comm);
  A.Export(*readA, exporter, Add);
  x.Export(*readx, exporter, Add);
  b.Export(*readb, exporter, Add);
  xexact.Export(*readxexact, exporter, Add);

  A.FillComplete();

  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;

  // E N D   O F   M A T R I X   C O N S T R U C T I O N

  // ============================= //
  // Construct RILU preconditioner //
  // ---------------=------------- //

  //  modify those parameters
  int    LevelFill = 0;
  int    Overlap = 2;
  double Athresh = 0.0;
  double Rthresh = 1.0;

  Ifpack_IlukGraph * Graph = 0;
  Ifpack_CrsRiluk * RILU = 0;

  Graph = new Ifpack_IlukGraph(A.Graph(), LevelFill, Overlap);
  Graph->ConstructFilledGraph();

  RILU = new Ifpack_CrsRiluk(*Graph);
  int initerr = RILU->InitValues(A);
  if (initerr!=0) cout << Comm << "*ERR* InitValues = " << initerr;

  RILU->Factor();

  // Define label for printing out during the solve phase
  string label = "Ifpack_CrsRiluk Preconditioner: LevelFill = " + toString(LevelFill) +
    " Overlap = " + toString(Overlap) +
    " Athresh = " + toString(Athresh) +
    " Rthresh = " + toString(Rthresh);
  RILU->SetLabel(label.c_str());

  // Here we create an AztecOO object
  AztecOO solver;
  solver.SetUserMatrix(&A);
  solver.SetLHS(&x);
  solver.SetRHS(&b);

  // Here we set the IFPACK preconditioner and specify few parameters

  solver.SetPrecOperator(RILU);

  int Niters = 1200;
  solver.SetAztecOption(AZ_kspace, Niters);
  solver.Iterate(Niters, 5.0e-6);

  delete RILU;
  delete Graph;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return 0 ;
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
      "--enable-ifpack\n"
      "--enable-aztecoo");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
#endif
