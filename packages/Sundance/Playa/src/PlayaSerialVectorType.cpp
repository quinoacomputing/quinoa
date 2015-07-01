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


#include "PlayaSerialVectorType.hpp"
#include "PlayaSerialVectorSpace.hpp"
#include "PlayaDenseSerialMatrixFactory.hpp"
#include "PlayaSerialGhostImporter.hpp"
#include "PlayaOut.hpp"

#include "Teuchos_RefCountPtr.hpp"

namespace Playa
{

using namespace Teuchos;

SerialVectorType::SerialVectorType()
{;}


RCP<const VectorSpaceBase<double> > 
SerialVectorType::createSpace(int dimension,
  int nLocal,
  const int* localIndices,
  const MPIComm& comm) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(nLocal < 0, std::runtime_error, "negative vector size n=" << nLocal);
  TEUCHOS_TEST_FOR_EXCEPTION(dimension != nLocal, std::runtime_error, 
    "nLocal=" << nLocal << " and dimension=" << dimension
    << " should be equal for a replicated space");

	return rcp(new SerialVectorSpace(dimension));
}

RCP<GhostImporter<double> > 
SerialVectorType::createGhostImporter(const VectorSpace<double>& space,
                                      int nGhost,
                                      const int* ghostIndices) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(dynamic_cast<const SerialVectorSpace*>(space.ptr().get())==0, std::runtime_error, "expected " 
    << space << " to be a SerialVectorSpace");
  return rcp(new SerialGhostImporter());
}

RCP<MatrixFactory<double> >
SerialVectorType::createMatrixFactory(const VectorSpace<double>& domain,
                                      const VectorSpace<double>& range) const
{
  RCP<MatrixFactory<double> > rtn 
    = rcp(new DenseSerialMatrixFactory(domain, range));

  return rtn;
}


VectorSpace<double> SerialVectorType
::createEvenlyPartitionedSpace(const MPIComm& /* comm */,
  int nLocal) const
{
  RCP<const VectorSpaceBase<double> > rtn = rcp(new SerialVectorSpace(nLocal));
  return rtn;
}

}




