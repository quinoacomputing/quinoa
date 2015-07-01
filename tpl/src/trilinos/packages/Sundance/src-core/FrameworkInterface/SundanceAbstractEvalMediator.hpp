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

#ifndef SUNDANCE_EVALMEDIATOR_H
#define SUNDANCE_EVALMEDIATOR_H


#include "SundanceDefs.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceObjectWithVerbosity.hpp"

namespace Sundance
{

using namespace Sundance;

class CoordExpr;
class CellDiameterExpr;
class CurveNormExpr;
class CellVectorExpr;
class MultiIndex; 
class DiscreteFuncElement;
class SparsitySuperset;

/**
 * Base class for evaluation mediator objects. 
 * Evaluation mediators are responsible
 * for evaluating those expressions whose
 * calculation must be delegated to the framework.
 */
class AbstractEvalMediator 
{
public:
  /** */
  AbstractEvalMediator(int verb=0);

  /** */
  virtual ~AbstractEvalMediator(){;}

  /** */
  void setVerb(int verb, int dfVerb) const 
    {verb_ = verb; dfVerb_=dfVerb;}

  /** */
  int verb() const {return verb_;}

  /** */
  int dfVerb() const {return dfVerb_;}

          

  /** Evaluate the given coordinate expression, putting
   * its numerical values in the given EvalVector. */
  virtual void evalCoordExpr(const CoordExpr* expr,
    RCP<EvalVector>& vec) const = 0 ;

  /** Evaluate the given discrete function, putting
   * its numerical values in the given EvalVector. */
  virtual void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
    const Array<MultiIndex>& mi,
    Array<RCP<EvalVector> >& vec) const = 0 ;

  /** Evaluate the given cell diameter expression, putting
   * its numerical values in the given EvalVector. */
  virtual void evalCellDiameterExpr(const CellDiameterExpr* expr,
    RCP<EvalVector>& vec) const = 0 ;

  /** Evaluates one component of the normal vector to a given parameterized curve
   * i.e. x,y or z component of that vector in 3D <br>
   * , this method is only in the CurveEvalMediator class implemented */
  virtual void evalCurveNormExpr(const CurveNormExpr* expr,
    RCP<EvalVector>& vec) const
    {
	  TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error,
		" EvalMediator::evalCurveNormExpr , not possible with the current EvalMediator (only with CurveEvalMediator)");
    }

  /** Evaluate the given cell vector expression, putting
   * its numerical values in the given EvalVector. */
  virtual void evalCellVectorExpr(const CellVectorExpr* expr,
    RCP<EvalVector>& vec) const = 0 ;


  /** Print evaluation results */
  virtual void showResults(std::ostream& os,
			   const RCP<SparsitySuperset>& sparsity,
			   const Array<RCP<EvalVector> >& vecResults,
			   const Array<double>& constantResults) const  ;
private:
  mutable int verb_;
  mutable int dfVerb_;
};
}


#endif
