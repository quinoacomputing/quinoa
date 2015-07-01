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

#include "SundanceFunctional.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceAssembler.hpp"


using namespace Sundance;
using namespace Teuchos;
using namespace Playa;


Functional::Functional(const Mesh& mesh, const Expr& integral, 
  const Playa::VectorType<double>& vecType)
  : mesh_(mesh),
    integral_(integral),
    bc_(),
    vecType_(vecType)
{
  
}

Functional::Functional(
  const Mesh& mesh, 
  const Expr& integral,  
  const Expr& essentialBC,
  const Playa::VectorType<double>& vecType)
  : mesh_(mesh),
    integral_(integral),
    bc_(essentialBC),
    vecType_(vecType)
{
  
}


LinearProblem Functional::linearVariationalProb(
  const Expr& var,
  const Expr& varEvalPts,
  const Expr& unk,
  const Expr& fixed,
  const Expr& fixedEvalPts) const 
{

  Array<Expr> zero(unk.size());
  for (int i=0; i<unk.size(); i++) 
    {
      Expr z = new ZeroExpr();
      zero[i] = z;
    }

  Expr unkEvalPts = new ListExpr(zero);

  Expr unkParams;
  Expr fixedParams;
  Expr unkParamValues;
  Expr fixedParamValues;

  RCP<EquationSet> eqn 
    = rcp(new EquationSet(integral_, bc_, 
                          tuple(var), tuple(varEvalPts),
                          tuple(unk), tuple(unkEvalPts), 
                          unkParams, unkParamValues,
                          tuple(fixed), tuple(fixedEvalPts)));

  RCP<Assembler> assembler 
    = rcp(new Assembler(mesh_, eqn, tuple(vecType_), tuple(vecType_), false));

  return LinearProblem(assembler);
}

NonlinearProblem Functional
::nonlinearVariationalProb(const Expr& var,
                           const Expr& varEvalPts,
                           const Expr& unk,
                           const Expr& unkEvalPts,
                           const Expr& fixed,
                           const Expr& fixedEvalPts) const
{

  Expr unkParams;
  Expr fixedParams;
  Expr unkParamValues;
  Expr fixedParamValues;

  RCP<EquationSet> eqn 
    = rcp(new EquationSet(integral_, bc_, 
                          tuple(var), tuple(varEvalPts),
                          tuple(unk), tuple(unkEvalPts), 
                          fixedParams, fixedParamValues,
                          tuple(fixed), tuple(fixedEvalPts)));

  RCP<Assembler> assembler 
    = rcp(new Assembler(mesh_, eqn, tuple(vecType_), tuple(vecType_), false));

  return NonlinearProblem(assembler, unkEvalPts);
}

FunctionalEvaluator Functional::evaluator(const Expr& var,
                                          const Expr& varEvalPts,
                                          const Expr& fixed,
                                          const Expr& fixedEvalPts) const 
{

  Expr unkParams;
  Expr fixedParams;
  Expr unkParamValues;
  Expr fixedParamValues;

  return FunctionalEvaluator(mesh_, integral_, bc_,
                             var, 
                             varEvalPts, 
                             fixed, fixedEvalPts,
                             vecType_);
}


FunctionalEvaluator Functional::evaluator(const Expr& var,
                                          const Expr& varEvalPts) const 
{
  return FunctionalEvaluator(mesh_, integral_, bc_,
                             var, 
                             varEvalPts,
                             vecType_);
}


namespace Sundance
{

double L2Norm(const Mesh& mesh, const CellFilter& domain,
  const Expr& f, const QuadratureFamily& quad,
  const WatchFlag& watch)
{
  Expr I2 = Integral(domain, f*f, quad, watch);

  return sqrt(evaluateIntegral(mesh, I2));
}


double H1Seminorm(
  const Mesh& mesh,
  const CellFilter& filter,
  const Expr& f,
  const QuadratureFamily& quad,
  const WatchFlag& watch)
{
  Expr grad = gradient(mesh.spatialDim());
  return L2Norm(mesh, filter, grad*f, quad, watch);
}
  
double H1Norm(
  const Mesh& mesh,
  const CellFilter& filter,
  const Expr& f,
  const QuadratureFamily& quad,
  const WatchFlag& watch)
{
  Expr grad = gradient(mesh.spatialDim());
  Expr g = grad*f;
  Expr I2 = Integral(filter, f*f + g*g, quad, watch);

  return sqrt(evaluateIntegral(mesh, I2));  
}
}
