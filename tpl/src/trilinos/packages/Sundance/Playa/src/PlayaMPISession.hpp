/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef PLAYA_MPISESSION_H
#define PLAYA_MPISESSION_H

/*! \file Playa_MPISession.hpp
  \brief A MPI utilities class, providing methods for initializing,
	finalizing, and querying the global MPI session
*/
#include "PlayaDefs.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

namespace Playa
{
/**
 * \brief This class provides methods for initializing, finalizing, 
 * and querying the global MPI session. 
 */
class MPISession
{
public:
  //! Initializer, calls MPI_Init() if necessary
  static void init(int* argc, void*** argv);

  //! Initializer, calls MPI_Init() if necessary
  static void init(int* argc, char*** argv);

  //! Returns the process rank relative to MPI_COMM_WORLD
  static int getRank() {return rank_;}

  //! Returns the number of processors in MPI_COMM_WORLD 
  static int getNProc() {return nProc_;}

  //! Finalizer, calls MPI_Finalize() if necessary
  static void finalize();

  /** Set to true if a message should be written by each processor
   * at startup. */
  static bool& showStartupMessage() {static bool rtn=false; return rtn;}

private:
  static int rank_;
  static int nProc_;

  /** private ctor, to be called only by init() */
  MPISession(int* argc, char*** argv);

public:
  /** dtor calls MPI finalize */
  ~MPISession();
};
}

#endif
