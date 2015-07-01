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

#ifndef SUNDANCE_EVALCONTEXT_H
#define SUNDANCE_EVALCONTEXT_H


#include "SundanceDefs.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_Utils.hpp"
#include <algorithm>


namespace Sundance
{
using namespace Teuchos;
using namespace Sundance;
using Sundance::Set;

/** 
 * Different contexts might require the same expression to be
 * evaluated to different orders of functional differentiation; for
 * example, in setting up a linear system, second-order derivatives
 * are required, but in evaluating a functional only zeroth derivs
 * are required. 
 * An EvaluationContext is used as a key to associate an evaluator and
 * its corresponding set of
 * functional derivatives with a context.
 *
 * They key consists of three parts: first, an integer identifier
 * indicating the caller, e.g., an assembler or functional evaluator,
 * second, a set indicating which orders of 
 * differentiation are required by the top level caller, and third,
 a region-quadrature combination.  
*/
class EvalContext
{
public:
  /** Empty ctor */
  EvalContext() : setupVerbosity_(0), evalSetupVerbosity_(0), 
                  maxDiffOrder_(0),
                  data_() {;}

  /** Construct with a region-quadrature combination and
   * an identifier of the construcing context. */
  EvalContext(const RegionQuadCombo& rqc,
    const Set<int>& needsDiffOrder,
    int contextID)
    : setupVerbosity_(0), evalSetupVerbosity_(0),
      maxDiffOrder_(*std::max_element(needsDiffOrder.begin(), needsDiffOrder.end())),
      data_(rcp(new OrderedTriple<Set<int>, int, RegionQuadCombo>(needsDiffOrder, contextID, rqc)))
    {}

  /** Set the verbosity level to be used during preprocessing 
   * of expressions in this context */
  void setSetupVerbosity(int v) const {setupVerbosity_ = v;}

  /** Return the verbosity level to be used during preprocessing 
   * of expressions in this context */
  int setupVerbosity() const {return setupVerbosity_;}

  /** Set the verbosity level to be used during setup of evaluators
   * for expressions in this context */
  void setEvalSetupVerbosity(int v) const {evalSetupVerbosity_ = v;}

  /** Get the verbosity level to be used during setup of evaluators
   * for expressions in this context */
  int evalSetupVerbosity() const {return evalSetupVerbosity_;}

  /** Comparison operator for use in maps */
  bool operator<(const EvalContext& other) const 
    {return *data_ < *other.data_;}
          
  /** Write to a std::string */
  std::string toString() const
    {return "EvalContext[diffOrder=" 
        + Teuchos::toString(data_->a())
        + ", id=" 
        + Teuchos::toString(data_->b())
        + ", " + data_->c().toString() + "]";}
          
  /** Write a short description to a std::string */
  std::string brief() const
    {return "EvalContext[diffOrder=" 
        + Teuchos::toString(data_->a())
        + ", id=" 
        + Teuchos::toString(data_->b())
        + "]";}

  /** */
  int topLevelDiffOrder() const {return maxDiffOrder_;}

  /** Indicate whether or not a given order of differentiation 
   * is needed in this context */
  bool needsDerivOrder(int order) const {return data_->a().contains(order);}
  

  /** Return a unique context ID */
  static int nextID() {static int rtn=0; return rtn++;}
private:
  mutable int setupVerbosity_;
  mutable int evalSetupVerbosity_;
  int maxDiffOrder_;
  RCP<OrderedTriple<Set<int>, int, RegionQuadCombo> > data_;
};

}


namespace std
{
/** \relates Sundance::EvalContext */
inline ostream& operator<<(std::ostream& os, 
  const Sundance::EvalContext& c)
{
  os << c.toString();
  return os;
}
}

namespace Teuchos
{
/** \relates Sundance::EvalContext */
inline std::string toString(const Sundance::EvalContext& h)
{return h.toString();}

}


#endif
