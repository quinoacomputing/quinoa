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

#ifndef SUNDANCE_FUNCTIONALEVALUATOR_H
#define SUNDANCE_FUNCTIONALEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaVectorType.hpp"
#include "SundanceFunctionalEvaluatorBase.hpp"

namespace Sundance
{

using namespace Playa;
using namespace Teuchos;

class Mesh;
class Assembler;

/** 
 *
 */
class FunctionalEvaluator : public FunctionalEvaluatorBase
{
public:
  /** */
  FunctionalEvaluator();

  /** */
  FunctionalEvaluator(const Mesh& mesh, 
    const Expr& integral);

  /** */
  FunctionalEvaluator(const Mesh& mesh, 
    const Expr& integral,
    const Expr& bcs,
    const Expr& var,
    const Expr& varEvalPts,
    const VectorType<double>& vectorType);

  /** */
  FunctionalEvaluator(const Mesh& mesh, 
    const Expr& integral,
    const Expr& bcs,
    const Expr& vars,
    const Expr& varEvalPts,
    const Expr& fields,
    const Expr& fieldValues,
    const VectorType<double>& vectorType);


  /** */
  double evaluate() const ;

  /** */
  Expr evalGradient(double& value) const ;

  /** */
  double fdGradientCheck(double h) const ;
          

private:

  /** */
  Vector<double> evalGradientVector(double& value) const ;
      
  /** */
  RCP<Assembler> assembler_;
      
  /** */
  mutable Expr varValues_;

  /** */
  VectorType<double> vecType_;
      
  /** */
  mutable Array<Vector<double> > gradient_;
      
};

/** */
double evaluateIntegral(const Mesh& mesh, const Expr& expr);


}


#endif
