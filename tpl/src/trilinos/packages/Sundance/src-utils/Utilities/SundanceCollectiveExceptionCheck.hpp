/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
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

#ifndef SUNDANCE_COLLECTIVEEXCEPTIONCHECK_H
#define SUNDANCE_COLLECTIVEEXCEPTIONCHECK_H

#include "SundanceDefs.hpp"
#include "PlayaMPIComm.hpp"


namespace Sundance
{
using Playa::MPIComm;


  /** Call this function upon catching an exception at a point before a collective
   * operation. This function will do an AllReduce in conjunction with calls
   * to either this function or its partner, checkForFailures(), on the
   * other processors. This procedure has the effect of communicating to the other
   * processors that an exception has been detected on this one. */
  void reportFailure(const MPIComm& comm);

  /** Call this function after exception-free completion of a
   * try/catch block preceding a collective operation. This function
   * will do an AllReduce in conjunction with calls to either this
   * function or its partner, reportFailure(), on the other
   * processors. If a failure has been reported by another processor, the
   * call to checkForFailures() will return true and an exception can be thrown. */
  bool checkForFailures(const MPIComm& comm);
}

#endif
