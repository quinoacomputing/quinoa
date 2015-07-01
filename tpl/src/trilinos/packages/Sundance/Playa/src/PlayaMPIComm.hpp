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

#ifndef PLAYA_MPICOMM_H
#define PLAYA_MPICOMM_H

/*! \file PlayaMPIComm.hpp
  \brief Object representation of a MPI communicator
*/

#include "PlayaDefs.hpp"
#include "PlayaMPIOp.hpp"
#include "PlayaMPIDataType.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#endif


namespace Playa
{


/**
 * \brief Object representation of an MPI communicator.
 *
 * At present, groups are not implemented so the only communicators
 * are MPI_COMM_WORLD and MPI_COMM_SELF
 */
class MPIComm
{
public:

  //! Empty constructor builds an object for MPI_COMM_WORLD
  MPIComm();

#ifdef HAVE_MPI
  //! Construct a MPIComm for a given MPI communicator
  MPIComm(MPI_Comm comm);
#endif

  //! Get an object representing MPI_COMM_WORLD 
  static MPIComm& world();
  //! Get an object representing MPI_COMM_SELF
  static MPIComm& self();

  //! Return process rank
  int getRank() const {return myRank_;}

  //! Return number of processors in the communicator
  int getNProc() const {return nProc_;}

  //! Synchronize all the processors in the communicator
  void synchronize() const ;

  //! @name Collective communications 
  //@{

  //! All-to-all gather-scatter
  void allToAll(void* sendBuf, int sendCount, const MPIDataType& sendType,
    void* recvBuf, int recvCount, 
    const MPIDataType& recvType) const ;

  //! Variable-length gather-scatter
  void allToAllv(void* sendBuf, int* sendCount, int* sendDisplacements,
    const MPIDataType& sendType,
    void* recvBuf, int* recvCount,
    int* recvDisplacements,
    const MPIDataType& recvType) const ;

  //! Do a collective operation, scattering the results to all processors
  void allReduce(void* input, void* result, int inputCount, 
    const MPIDataType& type,
    const MPIOp& op) const ;


  //! Gather to root 
  void gather(void* sendBuf, int sendCount, const MPIDataType& sendType,
    void* recvBuf, int recvCount, const MPIDataType& recvType,
    int root) const ;

  //! Gather variable-sized arrays to root 
  void gatherv(void* sendBuf, int sendCount, const MPIDataType& sendType,
    void* recvBuf, int* recvCount, int* displacements, 
    const MPIDataType& recvType, int root) const ;

  //! Gather to all processors
  void allGather(void* sendBuf, int sendCount, 
    const MPIDataType& sendType,
    void* recvBuf, int recvCount, 
    const MPIDataType& recvType) const ;

  //! Variable-length gather to all processors
  void allGatherv(void* sendBuf, int sendCount, 
    const MPIDataType& sendType,
    void* recvBuf, int* recvCount, int* recvDisplacements,
    const MPIDataType& recvType) const ;

  //! Broadcast 
  void bcast(void* msg, int length, 
    const MPIDataType& type, int src) const ;

  //@}

#ifdef HAVE_MPI
  //! Get the MPI_Comm communicator handle 
  MPI_Comm getComm() const {return comm_;}
#endif

  

  // errCheck() checks the return value of an MPI call and throws
  // a ParallelException upon failure.
  static void errCheck(int errCode, const std::string& methodName);

private:
#ifdef HAVE_MPI
  MPI_Comm comm_;
#endif

  int nProc_;
  int myRank_;

  /** common initialization function, called by all ctors */
  void init();

  /** Indicate whether MPI is currently running */
  int mpiIsRunning() const ;
};
}
#endif

