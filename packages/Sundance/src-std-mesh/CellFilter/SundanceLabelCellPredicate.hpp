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


#ifndef SUNDANCE_LABELCELLPREDICATE_H
#define SUNDANCE_LABELCELLPREDICATE_H

#include "SundanceDefs.hpp"
#include "SundanceCellPredicateBase.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
using namespace Teuchos;
  
/** 
 * LabelCellPredicate tests whether a cell's label is equal to
 * a reference label. 
 */
class LabelCellPredicate : public CellPredicateBase 
{
public:
  /** Construct with an integer label */
  LabelCellPredicate(int label) 
    : CellPredicateBase(), labelIndices_(makeSet(label)){;}
  /** Construct with an array of labels */
  LabelCellPredicate(const Array<int>& labels) 
    : CellPredicateBase(), labelIndices_(makeSet(labels)){;}
  /** Construct with a set of labels */
  LabelCellPredicate(const Set<int>& labels) 
    : CellPredicateBase(), labelIndices_(labels){;}

  /** virtual dtor */
  virtual ~LabelCellPredicate(){;}
      
  /** Test the predicate on a batch of cells */
  virtual void testBatch(const Array<int>& cellLID,
    Array<int>& results) const ;

  /** Write to XML */
  virtual XMLObject toXML() const ;

  /** comparison */
  virtual bool lessThan(const CellPredicateBase* other) const ;

  /** */
  virtual std::string description() const 
    {return "Label(" + labelIndices_.toString() + ")";}

  /* */
  GET_RCP(CellPredicateBase);

private:

  Set<int> labelIndices_;

};

}

#endif
