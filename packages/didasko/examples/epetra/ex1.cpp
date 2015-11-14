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

// Basic definition of communicator.
// This code should be run with at least two processes
// However, note that the SAME code works even if Epetra
// has been configured without MPI!

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  // define an Epetra communicator
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // get the proc ID of this process
  int MyPID = Comm.MyPID();

  // get the total number of processes
  int NumProc = Comm.NumProc();

  // output some information to std output
  cout << Comm << endl;

  // ======================== //
  // now some basic MPI calls //
  // ------------------------ //

  int    ivalue;
  double dvalue, dvalue2;
  double* dvalues;  dvalues  = new double[NumProc];
  double* dvalues2; dvalues2 = new double[NumProc];
  int root = 0;

  // equivalent to MPI_Barrier

  Comm.Barrier();

  if (MyPID == root) dvalue = 12.0;

  // On input, the root processor contains the list of values
  // (in this case, a single value). On exit, all processes will
  // have he same list of values. Note that all values must be allocated
  // vefore the broadcast

  // equivalent to  MPI_Broadcast

  Comm.Broadcast(&dvalue, 1, root);

  // as before, but with integer values. As C++ can bind to the appropriate
  // interface based on argument typing, the type of data is not required.

  Comm.Broadcast(&ivalue, 1, root);

  // equivalent MPI_Allgather

  Comm.GatherAll(dvalues, dvalues2, 1);

  // equivalent to MPI_Allreduce with MPI_SUM

  dvalue = 1.0*MyPID;

  Comm.SumAll( &dvalue, dvalues, 1);

  // equivalent to MPI_Allreduce with MPI_SUM

  Comm.MaxAll( &dvalue, dvalues, 1);

  // equiavant to MPI_Scan with MPI_SUM

  dvalue = 1.0 * MyPID;

  Comm.ScanSum(&dvalue, &dvalue2, 1);

  cout << "On proc " << MyPID << " dvalue2  = " << dvalue2 << endl;

  delete[] dvalues;
  delete[] dvalues2;

  // ======================= //
  // Finalize MPI and return //
  // ----------------------- //

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

} /* main */

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
      "--enable-epetra");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
#endif

