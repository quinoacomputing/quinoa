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

#ifndef SUNDANCE_TEMPSTACK_H
#define SUNDANCE_TEMPSTACK_H

#include "SundanceDefs.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceNoncopyable.hpp"
#include <stack>

namespace Sundance
{
using namespace Sundance;
/**
 * TempStack provides a stack of temporary variables for use during
 * evaluation. 
 *
 * During the course of evaluating an expression, it is often necessary
 * to create temporary variables. For example, in evaluating
 * \code
 * a += b*(c+d)
 * \endcode
 * it is required to create three temporaries (ignoring for
 * explanatory purposes any copies made in cases where one of
 * the vectors will be used elsewhere). We can see this
 * by breaking the operation
 * down into the following steps:
 * \code
 * 1. Create a temporary variable t1
 * 2. Evaluate expression b into t1
 * 3. Create a temporary variable t2
 * 4. Evaluate expression c into t2
 * 3. Create a temporary variable t3
 * 4. Evaluate expression d into t3
 * 5. Carry out t2 += t3
 * 6. Carry out t1 *= t2
 * 7. Carry out a += t1
 * \endcode
 * The number of temporaries required for a given expression
 * will depend on the graph of the expression. In general, we want to
 * create exactly as many temporaries as are needed, and reuse any
 * temporaries that are no longer needed. This is a well-known problem
 * in compiler design, and can be accomplished by maintaining a
 * stack of temporaries. When a new temporary is needed, it is popped
 * from the stack; if the stack is empty, a new temporary is allocated.
 * When a step of a calculation is done, any temporaries used are
 * put back on the stack for further use.
 */
class TempStack : public Noncopyable
{
public:
  /** Empty ctor */
  TempStack();

  /** Construct with an initial vector size */
  TempStack(int vecSize);

  /** Push vector data onto the stack */
  void pushVectorData(const RCP<Array<double> >& vecData) ;

  /** Pop vector data from the stack */
  RCP<Array<double> > popVectorData() ;

  /** Get a new vector (which will often reuse stack data) */
  RCP<EvalVector> popVector() 
    {return rcp(new EvalVector(this));}

  /** */
  void setVecSize(int vecSize) {vecSize_ = vecSize;}

  /** */
  void resetCounter() ;

  /** */
  int numVecsAccessed() const {return numVecsAccessed_;}

  /** */
  int numVecsAllocated() const {return numVecsAllocated_;}

  /** */
  int vecSize() const {return vecSize_;}

private:
          
  int vecSize_;

  std::stack<RCP<Array<double> > > stack_;

  int numVecsAllocated_;

  int numVecsAccessed_;
};
}

#endif
