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

#include "PlayaSerialVectorSpace.hpp"
#include "PlayaSerialVector.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaOut.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorSpaceImpl.hpp"
#endif




namespace Playa
{
using Teuchos::RCP;


SerialVectorSpace::SerialVectorSpace(int dim)
  : dim_(dim), comm_(MPIComm::self())
{}


// Overridden from VectorSpaceBase

RCP<VectorBase<double> >
SerialVectorSpace::createMember(const VectorSpace<double>& self) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(self.ptr().get() != this,
    InternalError, 
    "inconsistency between space and self-reference in SerialVectorSpace::createMember()");

  return rcp(new SerialVector(self));
}

bool SerialVectorSpace::isCompatible(const VectorSpaceBase<double>* other) const
{
  const SerialVectorSpace* svs = dynamic_cast<const SerialVectorSpace*>(other);
  if (svs == 0) return false;
  return this->dim() == svs->dim();
}

string SerialVectorSpace::description() const
{
  std::string rtn = "SerialVS[d=" + Teuchos::toString(this->dim()) + "]";
  return rtn;
}

}


