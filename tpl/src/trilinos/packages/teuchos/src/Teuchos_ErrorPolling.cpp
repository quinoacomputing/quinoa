// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

// Make deprecated warnings go away for compling this source file.   All of this
// code is deprecated so there is little point in flagging it as such!
#define TEUCHOS_DEPRECATED

#include "Teuchos_ErrorPolling.hpp"
#include "Teuchos_MPIComm.hpp"

namespace Teuchos
{
  void ErrorPolling::reportFailure(const MPIComm& comm)
  {
    if (isActive())
      {
        int myBad = 1;
        int anyBad = 0;
        comm.allReduce((void*) &myBad, (void*) &anyBad, 1, MPIComm::INT,
                       MPIComm::SUM);
      }
  }

  bool ErrorPolling::pollForFailures(const MPIComm& comm)
  {
    /* bypass if inactive */
    if (!isActive()) return true;

    int myBad = 0;
    int anyBad = 0;
    try
      {
        comm.allReduce((void*) &myBad, (void*) &anyBad, 1, MPIComm::INT,
                       MPIComm::SUM);
      }
    catch(const std::exception&)
      {
        return true;
      }
    return anyBad > 0;
  }
}



	





