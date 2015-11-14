
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

// define a vector of global dimension 5, split among 2 processes p0 and p1
// int the following way:
//  0 -- 1 -- 2 -- 3 -- 4
//  <-p0->    <----p1--->
// This code must be run with two processes

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  if (Comm.NumProc() != 2) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return(0);
  }

  int NumGlobalElements=5;   // global dimension of the problem
  int MyElements = 0;        // local dimension of the problem
  int *MyGlobalElements = 0; // local-to-global map
  int MyPID = Comm.MyPID();  // ID of this process

  if( Comm.NumProc() != 2 ) {
    cerr << "This code must be run with 2 processes\n";
    exit( EXIT_FAILURE );
  }

  // hardwired number of local elements and local-to-global map
  switch( MyPID ) {
    case 0:
      MyElements = 2;
      MyGlobalElements = new int[MyElements];
      MyGlobalElements[0] = 0;
      MyGlobalElements[1] = 1;
      break;
    case 1:
      MyElements = 3;
      MyGlobalElements = new int[MyElements];
      MyGlobalElements[0] = 2;
      MyGlobalElements[1] = 3;
      MyGlobalElements[2] = 4;
      break;
  }

  Epetra_Map Map(NumGlobalElements,MyElements,MyGlobalElements,0,Comm);

  cout << Map;

  delete[] MyGlobalElements;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);

} /* main */

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure Didasko with:\n"
      "--enable-epetra");

  return 0;
}

#endif
