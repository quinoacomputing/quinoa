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

#ifndef PLAYA_VECTORTYPEBASE_HPP
#define PLAYA_VECTORTYPEBASE_HPP

#include "PlayaHandle.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp" 
#include "PlayaMatrixFactory.hpp" 
#include "PlayaGhostImporter.hpp"

namespace Playa
{
using namespace Teuchos;

/**
 *
 */
template <class Scalar>
class VectorTypeBase
{
public:
  /** Virtual dtor */
  virtual ~VectorTypeBase() {;}

  /** create a distributed vector space.
   * @param dimension the dimension of the space 
   * @param nLocal number of indices owned by the local processor
   * @param locallyOwnedIndices array of indices owned by this processor  
   */
  virtual RCP<const VectorSpaceBase<Scalar> >
  createSpace(int dimension, 
    int nLocal,
    const int* locallyOwnedIndices,
    const MPIComm& comm) const = 0 ;
   

  /** Default implementation creates a vector space having 
   * nLocal elements on each processor. Serial types should override this
   * to produce a replicated space. */
  virtual VectorSpace<Scalar> 
  createEvenlyPartitionedSpace(const MPIComm& comm,
    int nLocal) const ;

  /**  
   * Create an importer for accessing ghost elements.
   * @param space the distributed vector space on which ghost elements
   * are to be shared
   * @param nGhost number of ghost elements needed by this processor
   * @param ghostIndices read-only C array of off-processor indices needed
   * by this processor.
   * @return A RCP to a GhostImporter object.
   */
  virtual RCP<GhostImporter<Scalar> > 
  createGhostImporter(const VectorSpace<Scalar>& space,
    int nGhost,
    const int* ghostIndices) const = 0 ;

    
  /**
   * Create a matrix factory of type compatible with this vector type,
   * sized according to the given domain and range spaces.
   */
  virtual RCP<MatrixFactory<Scalar> >
  createMatrixFactory(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range) const = 0 ;
    
    
};



/* Default implementation */
template <class Scalar> inline 
VectorSpace<Scalar> VectorTypeBase<Scalar>
::createEvenlyPartitionedSpace(const MPIComm& comm,
  int nLocal) const
{
  int rank = comm.getRank();
  int nProc = comm.getNProc();
  int dimension = nLocal * nProc;
  Array<int> locallyOwnedIndices(nLocal);
  int lowestLocalRow = rank*nLocal;
  for (int i=0; i<nLocal; i++)
  {
    locallyOwnedIndices[i] = lowestLocalRow + i;
  }
  return this->createSpace(dimension, nLocal, &(locallyOwnedIndices[0]), comm);
}

  
}

#endif
