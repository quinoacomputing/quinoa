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

#ifndef SUNDANCE_EVALMANAGER_H
#define SUNDANCE_EVALMANAGER_H

#include "SundanceDefs.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceAbstractEvalMediator.hpp"
#include "SundanceTempStack.hpp"
#include "SundanceNoncopyable.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Sundance
{
using namespace Sundance;
class CoordExpr;
class MultiIndex;
class DiscreteFuncElement;


/**
 * EvalManager provides methods for interfacing to the framework
 * through an AbstractEvalMediator and managing temporary variables
 * through a TempStack.
 *
 * If no mediator is set, std::string evaluations will be done 
 */
class EvalManager : public Noncopyable
{
public:
  /** Empty ctor */
  EvalManager();

  /** */
  void evalCoordExpr(const CoordExpr* expr,
    RCP<EvalVector>&  result) const ;

  /** */
  void evalCellDiameterExpr(const CellDiameterExpr* expr,
    RCP<EvalVector>&  result) const ;

  /** */
  void evalCurveNormExpr(const CurveNormExpr* expr,
    RCP<EvalVector>&  result) const ;

  /** */
  void evalCellVectorExpr(const CellVectorExpr* expr,
    RCP<EvalVector>&  result) const ;

  /** */
  void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
    const Array<MultiIndex>& mi,
    Array<RCP<EvalVector> >& result) const ;

  void showResults(std::ostream& os,
		   const RCP<SparsitySuperset>& sparsity,
		   const Array<RCP<EvalVector> >& vecResults,
		   const Array<double>& constantResults) const ;

  /** */
  void setMediator(const RCP<AbstractEvalMediator>& med) 
    {mediator_ = med;}

  /** */
  void setVerb(int verb) ;


  /** */
  int verb() const {return verb_;}

  /** */
  void setVecSize(int vecSize) {stack().setVecSize(vecSize);}
          

  /** Return a pointer to the mediator. We'll need the
   * mediator for computing framework-specific functions.
   */
  const AbstractEvalMediator* mediator() const {return mediator_.get();}

  /** */
  void setRegion(const EvalContext& region)
    {region_ = region;}

  /** */
  const EvalContext& getRegion() const {return region_;}

  /** */
  static TempStack& stack();

  /** */
  int getMaxDiffOrder() const ;


  /** */
  RCP<EvalVector> popVector() const ;

  /** */
  TEUCHOS_TIMER(coordEvalTimer, "coord function evaluation");

  /** */
  TEUCHOS_TIMER(discFuncEvalTimer, "discrete function evaluation");

private:
  int verb_;

  EvalContext region_;

  RCP<AbstractEvalMediator> mediator_;

};

}

#endif
