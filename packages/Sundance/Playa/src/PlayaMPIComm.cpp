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

#include "PlayaMPIComm.hpp"
#include "PlayaMPIDataType.hpp"
#include "PlayaMPIOp.hpp"
#include "PlayaErrorPolling.hpp"




namespace Playa
{
using namespace Teuchos;


MPIComm::MPIComm()
	:
#ifdef HAVE_MPI
	comm_(MPI_COMM_WORLD),
#endif
	nProc_(0), myRank_(0)
{
	init();
}

#ifdef HAVE_MPI
MPIComm::MPIComm(MPI_Comm comm)
	: comm_(comm), nProc_(0), myRank_(0)
{
	init();
}
#endif

int MPIComm::mpiIsRunning() const
{
  int mpiStarted = 0;
#ifdef HAVE_MPI
  MPI_Initialized(&mpiStarted);
#endif
  return mpiStarted;
}

void MPIComm::init()
{
#ifdef HAVE_MPI

  if (mpiIsRunning())
  {
    errCheck(MPI_Comm_rank(comm_, &myRank_), "Comm_rank");
    errCheck(MPI_Comm_size(comm_, &nProc_), "Comm_size");
  }
  else
  {
    nProc_ = 1;
    myRank_ = 0;
  }
	
#else
	nProc_ = 1;
	myRank_ = 0;
#endif
}

#ifdef USE_MPI_GROUPS /* we're ignoring groups for now */

MPIComm::MPIComm(const MPIComm& parent, const MPIGroup& group)
	:
#ifdef HAVE_MPI
	comm_(MPI_COMM_WORLD), 
#endif
	nProc_(0), myRank_(0)
{
#ifdef HAVE_MPI
	if (group.getNProc()==0)
  {
    rank_ = -1;
    nProc_ = 0;
  }
	else if (parent.containsMe())
  {
    MPI_Comm parentComm = parent.comm_;
    MPI_Group newGroup = group.group_;
			
    errCheck(MPI_Comm_create(parentComm, newGroup, &comm_), 
      "Comm_create");
			
    if (group.containsProc(parent.getRank()))
    {
      errCheck(MPI_Comm_rank(comm_, &rank_), "Comm_rank");
					
      errCheck(MPI_Comm_size(comm_, &nProc_), "Comm_size");
    }
    else
    {
      rank_ = -1;
      nProc_ = -1;
      return;
    }
  }
	else
  {
    rank_ = -1;
    nProc_ = -1;
  }
#endif
}

#endif /* USE_MPI_GROUPS */

MPIComm& MPIComm::world()
{
	static MPIComm w = MPIComm();
	return w;
}


MPIComm& MPIComm::self()
{
#ifdef HAVE_MPI
	static MPIComm w = MPIComm(MPI_COMM_SELF);
#else
	static MPIComm w = MPIComm();
#endif
	return w;
}


void MPIComm::synchronize() const 
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
    if (mpiIsRunning())
    {
      /* test whether errors have been detected on another proc before
       * doing the collective operation. */
      TEUCHOS_POLL_FOR_FAILURES(*this);
      /* if we're to this point, all processors are OK */
        
      errCheck(::MPI_Barrier(comm_), "Barrier");
    }
	}
	//mutex_.unlock();
#endif
}

void MPIComm::allToAll(void* sendBuf, int sendCount, 
  const MPIDataType& sendType,
  void* recvBuf, int recvCount, const MPIDataType& recvType) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Datatype mpiSendType = sendType.handle();
		MPI_Datatype mpiRecvType = recvType.handle();


    if (mpiIsRunning())
    {
      /* test whether errors have been detected on another proc before
       * doing the collective operation. */
      TEUCHOS_POLL_FOR_FAILURES(*this);
      /* if we're to this point, all processors are OK */
        
      errCheck(::MPI_Alltoall(sendBuf, sendCount, mpiSendType,
          recvBuf, recvCount, mpiRecvType,
          comm_), "Alltoall");
    }
	}
	//mutex_.unlock();
#else
  (void)sendBuf;
  (void)sendCount;
  (void)sendType;
  (void)recvBuf;
  (void)recvCount;
  (void)recvType;
#endif
}

void MPIComm::allToAllv(void* sendBuf, int* sendCount, 
  int* sendDisplacements, const MPIDataType& sendType,
  void* recvBuf, int* recvCount, 
  int* recvDisplacements, const MPIDataType& recvType) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Datatype mpiSendType = sendType.handle();
		MPI_Datatype mpiRecvType = recvType.handle();

    if (mpiIsRunning())
    {
      /* test whether errors have been detected on another proc before
       * doing the collective operation. */
      TEUCHOS_POLL_FOR_FAILURES(*this);
      /* if we're to this point, all processors are OK */		
        
      errCheck(::MPI_Alltoallv(sendBuf, sendCount, sendDisplacements, mpiSendType,
          recvBuf, recvCount, recvDisplacements, mpiRecvType,
          comm_), "Alltoallv");
    }
	}
	//mutex_.unlock();
#else
  (void)sendBuf;
  (void)sendCount;
  (void)sendDisplacements;
  (void)sendType;
  (void)recvBuf;
  (void)recvCount;
  (void)recvDisplacements;
  (void)recvType;
#endif
}

void MPIComm::gather(void* sendBuf, int sendCount, const MPIDataType& sendType,
  void* recvBuf, int recvCount, const MPIDataType& recvType,
  int root) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Datatype mpiSendType = sendType.handle();
		MPI_Datatype mpiRecvType = recvType.handle();


    if (mpiIsRunning())
    {
      /* test whether errors have been detected on another proc before
       * doing the collective operation. */
      TEUCHOS_POLL_FOR_FAILURES(*this);
      /* if we're to this point, all processors are OK */
        
      errCheck(::MPI_Gather(sendBuf, sendCount, mpiSendType,
          recvBuf, recvCount, mpiRecvType,
          root, comm_), "Gather");
    }
  }
	//mutex_.unlock();
#endif
}

void MPIComm::gatherv(void* sendBuf, int sendCount, const MPIDataType& sendType,
  void* recvBuf, int* recvCount, int* displacements, const MPIDataType& recvType,
  int root) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Datatype mpiSendType = sendType.handle();
		MPI_Datatype mpiRecvType = recvType.handle();
		
    if (mpiIsRunning())
    {
      /* test whether errors have been detected on another proc before
       * doing the collective operation. */
      TEUCHOS_POLL_FOR_FAILURES(*this);
      /* if we're to this point, all processors are OK */
        
      errCheck(::MPI_Gatherv(sendBuf, sendCount, mpiSendType,
          recvBuf, recvCount, displacements, mpiRecvType,
          root, comm_), "Gatherv");
    }
	}
	//mutex_.unlock();
#endif
}

void MPIComm::allGather(void* sendBuf, int sendCount, const MPIDataType& sendType,
  void* recvBuf, int recvCount, 
  const MPIDataType& recvType) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Datatype mpiSendType = sendType.handle();
		MPI_Datatype mpiRecvType = recvType.handle();
		
    if (mpiIsRunning())
    {
      /* test whether errors have been detected on another proc before
       * doing the collective operation. */
      TEUCHOS_POLL_FOR_FAILURES(*this);
      /* if we're to this point, all processors are OK */
        
      errCheck(::MPI_Allgather(sendBuf, sendCount, mpiSendType,
          recvBuf, recvCount, 
          mpiRecvType, comm_), 
        "AllGather");
    }
	}
	//mutex_.unlock();
#endif
}


void MPIComm::allGatherv(void* sendBuf, int sendCount, const MPIDataType& sendType,
  void* recvBuf, int* recvCount, 
  int* recvDisplacements,
  const MPIDataType& recvType) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Datatype mpiSendType = sendType.handle();
		MPI_Datatype mpiRecvType = recvType.handle();
    
    if (mpiIsRunning())
    {
      /* test whether errors have been detected on another proc before
       * doing the collective operation. */
      TEUCHOS_POLL_FOR_FAILURES(*this);
      /* if we're to this point, all processors are OK */
        
      errCheck(::MPI_Allgatherv(sendBuf, sendCount, mpiSendType,
          recvBuf, recvCount, recvDisplacements,
          mpiRecvType, 
          comm_), 
        "AllGatherv");
    }
	}
	//mutex_.unlock();
#endif
}


void MPIComm::bcast(void* msg, int length, 
  const MPIDataType& type, int src) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
    if (mpiIsRunning())
    {
      /* test whether errors have been detected on another proc before
       * doing the collective operation. */
      TEUCHOS_POLL_FOR_FAILURES(*this);
      /* if we're to this point, all processors are OK */
        
      MPI_Datatype mpiType = type.handle();
      errCheck(::MPI_Bcast(msg, length, mpiType, src, 
          comm_), "Bcast");
    }
	}
	//mutex_.unlock();
#endif
}

void MPIComm::allReduce(void* input, void* result, int inputCount, 
    const MPIDataType& type,
    const MPIOp& op) const 
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Op mpiOp = op.handle();
		MPI_Datatype mpiType = type.handle();
		
    if (mpiIsRunning())
    {
      errCheck(::MPI_Allreduce(input, result, inputCount, mpiType,
          mpiOp, comm_), 
        "Allreduce");
    }
	}
	//mutex_.unlock();
#endif
}


void MPIComm::errCheck(int errCode, const std::string& methodName)
{
  TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
    "MPI function MPI_" << methodName 
    << " returned error code=" << errCode);
}



}
