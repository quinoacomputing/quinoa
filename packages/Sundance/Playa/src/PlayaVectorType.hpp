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

#ifndef PLAYA_VECTORTYPE_HPP
#define PLAYA_VECTORTYPE_HPP

#include "PlayaHandle.hpp"
#include "PlayaVectorTypeBase.hpp"
#include "PlayaVectorSpaceDecl.hpp"

namespace Playa
{
using namespace Teuchos;

/**
 * \brief Vector type objects are used by the application code to create
 * vector spaces and operators of a given type.
 */
template <class Scalar>
class VectorType : public Playa::Handle<VectorTypeBase<Scalar> >
{
public:
  HANDLE_CTORS(VectorType<Scalar>, VectorTypeBase<Scalar>);
   

  /** \brief Create a vector space having nLocal elements on each 
      processor of the specified communicator. */
  VectorSpace<Scalar> createEvenlyPartitionedSpace(const MPIComm& comm,
    int nLocal) const ;

  /** \brief Create a distributed vector space with an arbitrary
   * user-specified distribution of elements.
   * @param dimension the dimension of the space
   * @param nLocal number of indices owned by the local processor
   * @param locallyOwnedIndices array of indices owned by this processor  
   * @param comm the MPI communicator over which the space is to be distributed.
   */
  VectorSpace<Scalar> createSpace(int dimension, 
    int nLocal,
    const int* locallyOwnedIndices,
    const MPIComm& comm) const ;



  /** 
   * \brief Create an importer for ghost (off-processor) elements.
   **/
  RCP<GhostImporter<Scalar> > 
  createGhostImporter(const VectorSpace<Scalar>& space,
    int nGhost,
    const int* ghostIndices) const ;

  /**
   * \brief Create a matrix factory of type compatible with this vector type,
   * sized according to the given domain and range spaces.
   */
  virtual RCP<MatrixFactory<Scalar> >
  createMatrixFactory(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range) const ;
                                                      
};




template <class Scalar> inline 
VectorSpace<Scalar> VectorType<Scalar>::createSpace(int dimension,
  int nLocal,
  const int* locallyOwnedIndices, const MPIComm& comm) const
{
  return this->ptr()->createSpace(dimension, nLocal, locallyOwnedIndices, comm);
}

template <class Scalar> inline 
VectorSpace<Scalar> VectorType<Scalar>
::createEvenlyPartitionedSpace(const MPIComm& comm,
  int nLocal) const
{
  return this->ptr()->createEvenlyPartitionedSpace(comm, nLocal);
}


template <class Scalar> inline 
RCP<GhostImporter<Scalar> > 
VectorType<Scalar>::createGhostImporter(const VectorSpace<Scalar>& space,
  int nGhost,
  const int* ghostIndices) const
{
  return this->ptr()->createGhostImporter(space, nGhost, ghostIndices);
}

template <class Scalar> inline
RCP<MatrixFactory<Scalar> >
VectorType<Scalar>::createMatrixFactory(const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range) const
{
  return this->ptr()->createMatrixFactory(domain, range);
}

  
}

#endif
