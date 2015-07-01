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

#include "PlayaVectorDecl.hpp"
#include "PlayaSerialVector.hpp"
#include "PlayaSerialVectorSpace.hpp"
#include "Teuchos_Assert.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using namespace Teuchos;
using namespace Playa;

SerialVector::SerialVector(const VectorSpace<double>& vs)
  : SingleChunkVector<double>(),
    vecSpace_(vs),
    data_(vs.dim()),
    dim_(vs.dim())
{
  const SerialVectorSpace* rvs 
    = dynamic_cast<const SerialVectorSpace*>(vs.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(rvs==0, std::runtime_error,
    "could not cast vector space to SerialVectorSpace in "
    "SerialVector ctor");
}


void SerialVector::setElement(int index, const double& value)
{
  data_[index] = value;
}

void SerialVector::addToElement(int index, const double& value)
{
  data_[index] += value;
}

void SerialVector::setElements(int numElems, const int* globalIndices,
  const double* values)
{
  for (int i=0; i<numElems; i++)
  {
    data_[globalIndices[i]] = values[i];
  }
}

void SerialVector::addToElements(int numElems, const int* globalIndices,
  const double* values)
{
  for (int i=0; i<numElems; i++)
  {
    data_[globalIndices[i]] += values[i];
  }
}

const SerialVector* SerialVector::getConcrete(const Vector<double>& x)
{
  const SerialVector* rtn = dynamic_cast<const SerialVector*>(x.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPT(rtn==0);
  return rtn;
}

SerialVector* SerialVector::getConcrete(Vector<double>& x)
{
  SerialVector* rtn = dynamic_cast<SerialVector*>(x.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPT(rtn==0);
  return rtn;
}

void SerialVector::finalizeAssembly()
{
  // no-op
}

void SerialVector::getElements(const int* globalIndices, int numElems,
  Array<double>& elems) const
{
  elems.resize(numElems);
  for (int i=0; i<numElems; i++)
  {
    elems[i] = (*this)[globalIndices[i]];
  }
}

std::string SerialVector::description() const 
{
  std::ostringstream oss;
  oss << "SerialVector[dim=" << dim_ << "]" ;
  return oss.str();
}

