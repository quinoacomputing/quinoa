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

#ifndef SUNDANCE_EVALUATOR_H
#define SUNDANCE_EVALUATOR_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceSparsitySubset.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultiIndex.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

namespace Sundance 
{
class CoordExpr;

class EvalContext;
  


class EvalManager;

/**
 * Base class for evaluator objects. Each EvaluatableExpr type will 
 * have an associated Evaluator subtype.
 */
class Evaluator : public ObjectWithClassVerbosity<Evaluator>
{
public:
  /** */
  Evaluator();

  /** */
  virtual ~Evaluator(){;}

  /** 
   * Client-level evaluation method. Computes new results on the
   * first call, makes copies on subsequent calls up to the last client, 
   * and finally returns the original result vector upon the 
   * last client's call. 
   */
  void eval(const EvalManager& mgr,
    Array<double>& constantResults,
    Array<RCP<EvalVector> >& vectorResults) const ;

  /** Reset the number of calls to zero. This should be called
   * at the beginning of every new evaluation cycle. */
  virtual void resetNumCalls() const {numCalls_=0;}

  /** */
  virtual void 
  internalEval(const EvalManager& mgr,
    Array<double>& constantResults,
    Array<RCP<EvalVector> >& vectorResults) const = 0 ;

  /** Add one to the number of clients. */
  void addClient() {numClients_++;}

  /** */
  void addConstantIndex(int index, int constantIndex);

  /** */
  void addVectorIndex(int index, int vectorIndex);

      

  /** */
  const Sundance::Map<int, int>& constantIndexMap() const 
    {return constantIndexMap_;}

  /** */
  const Sundance::Map<int, int>& vectorIndexMap() const 
    {return vectorIndexMap_;}
protected:

  /** Return the number of clients that will require results
   * from this evaluator */
  int numClients() const {return numClients_;}

  /** */
  bool isOne(int x) const {return x==1;}

  /** */
  bool isOne(const double& x) const {return isZero(x-1.0);}

  /** */
  bool isZero(const double& x) const {return fabs(x-0.0)<1.0e-15;}

  /** */
  const Array<int>& constantIndices() const {return constantIndices_;}

  /** */
  const Array<int>& vectorIndices() const {return vectorIndices_;}


private:
  int numClients_;

  mutable int numCalls_;

  mutable Array<RCP<EvalVector> > vectorResultCache_;

  mutable Array<double> constantResultCache_;

  Sundance::Map<int, int> constantIndexMap_;

  Sundance::Map<int, int> vectorIndexMap_;

  Array<int> vectorIndices_;

  Array<int> constantIndices_;
};


    

}

#endif
