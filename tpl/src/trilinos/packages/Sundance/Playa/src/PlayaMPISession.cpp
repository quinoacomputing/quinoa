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

#include "PlayaMPISession.hpp"
#include "PlayaMPIDataType.hpp"
#include "PlayaMPIOp.hpp"
#include "Teuchos_Assert.hpp"

namespace Playa
{
using namespace Teuchos;

int MPISession::rank_ = 0 ;
int MPISession::nProc_ = 1 ;

MPISession::MPISession(int* argc, char*** argv)
{
  #ifdef HAVE_MPI
  /* initialize MPI */
	int mpiHasBeenStarted = 0;
	MPI_Initialized(& mpiHasBeenStarted);
	int mpierr = 0 ;
	if (!mpiHasBeenStarted)
		{
			mpierr = ::MPI_Init (argc, (char ***) argv);
      TEUCHOS_TEST_FOR_EXCEPTION(mpierr != 0, std::runtime_error,
                         "Error code=" << mpierr 
                         << " detected in MPI_Init()");
		}
	
	/* find rank */
	mpierr = ::MPI_Comm_rank (MPI_COMM_WORLD, &rank_);
	TEUCHOS_TEST_FOR_EXCEPTION(mpierr != 0, std::runtime_error,
                     "Error code=" << mpierr 
                     << " detected in MPI_Comm_rank()");

	/* find number of procs */
	mpierr = ::MPI_Comm_size (MPI_COMM_WORLD, &nProc_);

	TEUCHOS_TEST_FOR_EXCEPTION(mpierr != 0, std::runtime_error,
                     "Error code=" << mpierr 
                     << " detected in MPI_Comm_size()");

  /* get machine name */
  int nameLen;
	char procName[MPI_MAX_PROCESSOR_NAME];
  mpierr = ::MPI_Get_processor_name(procName,&nameLen);

  TEUCHOS_TEST_FOR_EXCEPTION(mpierr != 0, std::runtime_error,
                     "Error code=" << mpierr 
                     << " detected in MPI_Get_processor_name()");

  if (showStartupMessage())
    {
      std::cerr << "Playa::MPISession::init() started processor " 
           << procName << std::endl;
    }
  else
    {
#else
  std::cerr << "Playa::MPISession::init() started serial run" << std::endl;
#endif
#ifdef HAVE_MPI
    }
#endif
}


void MPISession::init(int* argc, char*** argv)
{
  static RCP<MPISession> session = rcp(new MPISession(argc, argv));
}

void MPISession::init(int* argc, void*** argv)
{
  MPISession::init(argc, (char***)argv);
}

void MPISession::finalize()
{
  // nothing to do here
}


MPISession::~MPISession()
{
#ifdef HAVE_MPI
  int mpiHasBeenFinalized = 0;
	MPI_Finalized(& mpiHasBeenFinalized);

  if (!mpiHasBeenFinalized)
  {
#define BROKEN_MPI_REGISTRIES
#ifndef BROKEN_MPI_REGISTRIES
    MPIDataType::clearTypeRegistry();
    MPIOp::clearOpRegistry();
#endif
    int mpierr = ::MPI_Finalize();
    
    TEUCHOS_TEST_FOR_EXCEPTION(mpierr != 0, std::runtime_error,
      "Error code=" << mpierr << " detected in MPI_Finalize()");
  }
#endif
}

}
