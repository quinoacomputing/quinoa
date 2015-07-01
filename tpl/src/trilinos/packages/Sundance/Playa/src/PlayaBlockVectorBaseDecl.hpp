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

#ifndef PLAYA_BLOCKVECTORBASEDECL_HPP
#define PLAYA_BLOCKVECTORBASEDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorBaseDecl.hpp"
#include "Teuchos_Describable.hpp"

namespace Playa
{

template <class Scalar> class Vector;
using Teuchos::Describable;

/** 
 * Base class for blocked vectors 
 */
template <class Scalar>
class BlockVectorBase : public VectorBase<Scalar>,
                        public Describable
{
public:
  /** */
  BlockVectorBase() : currentBlock_(0) {}

  /** */
  virtual ~BlockVectorBase(){}

  /** */
  virtual void setBlock(int b, const Vector<Scalar>& block) = 0 ;

  /** */
  virtual const Vector<Scalar>& getBlock(int b) const = 0 ;

  /** */
  virtual Vector<Scalar> getNonConstBlock(int b) = 0 ;

  /** */
  virtual ConstDataChunk<Scalar> nextConstChunk() const ;

  /** */
  virtual NonConstDataChunk<Scalar> nextChunk() ;

  /** */
  virtual bool hasMoreChunks() const ;

  /** */
  virtual void rewind() const ;

  /** */
  virtual std::string description() const ;

  /** */
  virtual void update(const Scalar& alpha, const VectorBase<Scalar>* other,
    const Scalar& gamma);

  /** */
  virtual void update(
    const Scalar& alpha, const VectorBase<Scalar>* x,
    const Scalar& beta, const VectorBase<Scalar>* y,
    const Scalar& gamma) ;

  /** */
  virtual void update(
    const Scalar& alpha, const VectorBase<Scalar>* x,
    const Scalar& beta, const VectorBase<Scalar>* y,
    const Scalar& gamma, const VectorBase<Scalar>* z,
    const Scalar& delta) ;

  /** */
  virtual Scalar dot(const VectorBase<Scalar>* other) const ;

  /** */
  virtual Scalar norm2() const ;
  
private:
  mutable int currentBlock_;
};







/** */





      
      
  
  
}

#endif
