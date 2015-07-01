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

#ifndef PLAYA_BLOCKITERATORDECL_HPP
#define PLAYA_BLOCKITERATORDECL_HPP

#include "PlayaDefs.hpp"
#include "Teuchos_Array.hpp"
#include <deque> 

namespace Playa
{

using Teuchos::RCP;
using Teuchos::Array;

/* Forward declare VectorSpaceBase */
template <class Scalar> 
class VectorSpaceBase;

/* Forward declare VectorSpace */
template <class Scalar> 
class VectorSpace;

/* Forward declare Vector */
template <class Scalar> 
class Vector;


/** 
 * BlockIterator can locate a block within an arbitrarily nested block
 * vector space.
 */
template <class Scalar>
class BlockIterator
{
public:

  friend class VectorSpaceBase<Scalar>;
  friend class VectorSpace<Scalar>;
  friend class Vector<Scalar>;

  /** Empty ctor */
  BlockIterator() : space_(),
                    index_(), 
                    atEnd_(true){}

  /** Compare two iterators */
  bool operator==(const BlockIterator<Scalar>& other) const ;

  /** Compare two iterators */
  bool operator!=(const BlockIterator<Scalar>& other) const 
    {
      return !operator==(other);
    }

  /** Compare two iterators */
  bool operator<(const BlockIterator<Scalar>& other) const ;

  /** Advance the block iterator */
  BlockIterator<Scalar> operator++(int);

  /** */
  const VectorSpace<Scalar>& space() const ;
  
  /** Print the iterator */
  std::ostream& print(std::ostream& os) const ;

  /** */
  const std::deque<int>& blockIndex() const {return index_;}

  /** */
  bool atEnd() const {return atEnd_;}

  /** Set to true (by doing BlockIterator<Scalar>::debug()=true) to
   * trace the iterator */
  static bool& debug() 
    {static bool rtn = false; return rtn;}

protected:

  /** Build an index pointing to the leftmost entry in the 
   * given block space */
  void goToStart(const VectorSpace<Scalar>& space,
    std::deque<int>& pos) const ;

  /** Constructor is private: the construction is always done inside
   * the begin and end methods of vector space. */
  BlockIterator(const VectorSpace<Scalar>& space,
    bool atEnd);

  /** Advance the index at the specified level. Return false if
   * no further advance is possible */
  bool advance(int level);

private:
  /** Store the VC in an RCP so we can do forward declaration of VS */
  RCP<const VectorSpace<Scalar> > space_;
  std::deque<int> index_;
  bool atEnd_;
};

}

namespace Playa
{
/** \relates BlockIterator */
template <class Scalar> inline
std::ostream& operator<<(std::ostream& os, 
  const Playa::BlockIterator<Scalar>& i)
{
  return i.print(os);
}
}


#endif
