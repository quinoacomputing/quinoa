/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
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

#ifndef SUNDANCE_CELLITERATOR_H
#define SUNDANCE_CELLITERATOR_H


#include "SundanceDefs.hpp"
#include "SundanceSet.hpp"
#include "SundanceMap.hpp"
#include "SundanceCellType.hpp"
#include "SundanceCellReordererImplemBase.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMesh.hpp"

namespace Sundance
{
using namespace Teuchos;

/**
 * CellIterator is an iterator for walking through cell sets.
 * It satisfies the requirements for an input iterator in STL.
 *
 * This class design violates the usual rules of good OO style:
 * it has polymorphic behavior packed into a single class. The
 * justification for this decision is to avoid the expense
 * of the clone() operations that would be required by the
 * copy ctor for iterators were 
 * a polymorphic class heirarchy used. 
 *
 * Two cell set types exist: explicit, where the member cells LIDs
 * are enumerated in a physical Set<int> object, and implicit,
 * where no physical Set is made, rather, the sequence of cell LIDs
 * is obtained through some scheme of walking the mesh. 
 *
 * The only cell sets that can represented implicitly are
 * the set of all cells of a given dimension.
 * 
 * \see CellSet, CellFilter, CellIteratorPos
 */
class CellIterator : public std::iterator<std::input_iterator_tag, int>
{
public:

      
  /** 
   * CellIteratorPos is used to specify whether a new CellIterator
   * is positioned at the beginning or end of a set.
   */
  enum CellIteratorPos {Begin, End};

  /** Empty ctor */
  CellIterator();

  /** Copy ctor */
  CellIterator(const CellIterator& other);

  /** Construct an implicit iterator for walking all cells of a given
   * dimension on the given mesh. */
  CellIterator(const Mesh& mesh, int cellDim, CellIteratorPos pos);

  /** Construct an explicit iterator for walking an explicitly
   * enumerated set of cells. */
  CellIterator(const Set<int>* cells, CellIteratorPos pos);

  /** */
  CellIterator& operator=(const CellIterator& other);
      
  /** Dereferencing operator */
  const int& operator*() const 
    {
      if (isImplicit_) return currentLID_;
      else return *iter_;
    }
      
  /** Postfix increment: advances iterator and returns previous value  */
  CellIterator operator++(int) 
    {
      CellIterator old = *this;
      advance();
      return old;
    }
      

  /** Prefix increment: advances iterator, returning new value */
  CellIterator& operator++()
    {
      advance();
      return *this;
    }

  /** */
  bool operator==(const CellIterator& other) const 
    {
      if (isImplicit_)
      {
        return currentLID_ == other.currentLID_;
      }
      else
      {
        return iter_ == other.iter_;
      }
    }

  /** */
  bool operator!=(const CellIterator& other) const 
    {
      return !(*this == other);
    }

      
private:

  /** Advance the iterator */
  void advance()
    {
      if (isImplicit_) 
      {
        if (reorderer_ != 0) 
        {
          currentLID_ = reorderer_->advance(currentLID_);
        }
        else currentLID_++;
      }
      else iter_++;
    }
      
  /** Flag indicating whether this iterator is implicit */
  bool isImplicit_;
      
  /** The LID to which this iterator is currently pointing.
   * Used only for implicit iterators. */
  int currentLID_;

  /** Unmanaged pointer to the reorderer used for walking 
   * implicit cell sets. Used only for implicit iterators. */
  const CellReordererImplemBase* reorderer_; 

  /** iterator for enumerated cells.
   * Used only for explicit iterators. */
  Set<int>::const_iterator iter_;

  /** */
  static Set<int> dummy() {static Set<int> rtn; return rtn;}
};
}


#endif
