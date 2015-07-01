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

#ifndef PLAYA_BLOCKVECTORSPACEBASEDECL_HPP
#define PLAYA_BLOCKVECTORSPACEBASEDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorSpaceBaseDecl.hpp"

namespace Playa
{

/**
 * This is a base class for a blocked vector space. It assumes nothing about
 * the physical storage of the blocks. 
 *
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar>
class BlockVectorSpaceBase : public VectorSpaceBase<Scalar>
{
public:
  /** virtual dtor */
  virtual ~BlockVectorSpaceBase() {;}

  /** Compute dimension. The default implementation sums the dimensions
   * of all blocks. */
  virtual int dim() const ;

  /** Count the locally owned elements */
  virtual int numLocalElements() const ;

  /** Check compatibility with another space. The default implementation
   * loops over blocks checking compatiblilty at each index. */
  virtual bool isCompatible(const VectorSpaceBase<Scalar>* other) const ;

  /** Get a block specified by an integer index. This function
   * should hrow an exception if the index is out of range 
   */
  virtual const VectorSpace<Scalar>& getBlock(int b) const = 0 ; 

  /** Write a description by recursivle describing the blocks */
  virtual std::string description() const ;

};






}

#endif
