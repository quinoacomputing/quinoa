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

#ifndef PLAYA_VECTORSPACEDECL_HPP
#define PLAYA_VECTORSPACEDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaVectorSpaceBaseDecl.hpp"
#include "PlayaBlockIteratorDecl.hpp"

namespace Playa
{
using namespace Teuchos;

/**
 *  User-level VectorSpace class.
 */
template <class Scalar>
class VectorSpace : public Playa::Handle< const VectorSpaceBase<Scalar> >
{
public:
  HANDLE_CTORS(VectorSpace<Scalar>, const VectorSpaceBase<Scalar>);
    
  /** Create a new element of this vector space */
  Vector<Scalar>  createMember() const ;

  /** Return the dimension of the space */
  int dim() const {return this->ptr()->dim();}

  /** Return the lowest global index accessible on this processor */
  int baseGlobalNaturalIndex() const ;

  /** Return the number of elements owned by this processor */
  int numLocalElements() const ;

  /** Return the MPI communicator */
  const MPIComm& comm() const {return this->ptr()->comm();}

  /** Check compatibility with another space. */
  bool isCompatible(const VectorSpace<Scalar>& vecSpc) const; 


  /** test equality between two spaces */
  bool operator==(const VectorSpace<Scalar>& other) const ;


  /** test inequality of two spaces */
  bool operator!=(const VectorSpace<Scalar>& other) const ;


  /** test whether the space contains a given vector */
  bool contains(const Vector<Scalar>& vec) const ;


  /** return the number of subblocks at the highest level. */
  int numBlocks() const ;

  /** indicate whether I am a block vector space */
  bool isBlockSpace() const ;

  /** get the i-th subblock */
  const VectorSpace<Scalar>& getBlock(int i) const ;

  /** get a subblock as specified by a block iterator */
  const VectorSpace<Scalar>& getBlock(const BlockIterator<Scalar>& iter) const ;

  /** get a subblock as specified by a deque of indices */
  const VectorSpace<Scalar>& getBlock(const std::deque<int>& iter) const ;

  /** */
  BlockIterator<Scalar> beginBlock() const ;

  /** */
  BlockIterator<Scalar> endBlock() const ;

  /** */
  int mapToGNI(const BlockIterator<Scalar>& b, int indexWithinBlock) const ;

  /** */
  bool containsGNI(int gni) const ;

  /** */
  void getBlockAndOffsetFromGNI(int gni,
    BlockIterator<Scalar>& block, int& indexWithinBlock) const ;
  
protected:
  

};

#define PLAYA_CHECK_SPACES(space1, space2) \
  TEUCHOS_TEST_FOR_EXCEPTION(!space1.isCompatible(space2), std::runtime_error, \
    "incompatible spaces " << space1 << " and " << space2)


template <class Scalar>
STREAM_OUT(VectorSpace<Scalar>)



}


#endif
