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

#include "SundanceFunctionalEvaluator.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace Playa;


using std::endl;
using std::setw;


namespace Sundance
{
double evaluateIntegral(const Mesh& mesh, const Expr& expr)
{
  FunctionalEvaluator eval(mesh, expr);
  return eval.evaluate();
}
}

FunctionalEvaluator::FunctionalEvaluator()
  : assembler_(),
    varValues_(),
    vecType_(),
    gradient_()
{}

FunctionalEvaluator::FunctionalEvaluator(const Mesh& mesh,
  const Expr& integral)
  : assembler_(),
    varValues_(),
    vecType_(),
    gradient_(1)
{
  Array<Expr> fields;
  Expr bcs;
  Expr params;

  
  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(integral, bcs, params, params, fields, fields));
  
  
  assembler_ = rcp(new Assembler(mesh, eqnSet));
}


FunctionalEvaluator::FunctionalEvaluator(const Mesh& mesh,
  const Expr& integral,
  const Expr& bcs,
  const Expr& var,
  const Expr& varValues,
  const VectorType<double>& vectorType)
  : assembler_(),
    varValues_(varValues),
    vecType_(vectorType),
    gradient_(1)
{
  Array<Expr> v = tuple(var.flatten());
  Array<Expr> v0 = tuple(varValues.flatten());
  Array<Expr> fixed;
  Expr params;
  
  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(integral, bcs, v, v0, params, params, fixed, fixed));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vectorType), tuple(vectorType), false));
}


FunctionalEvaluator::FunctionalEvaluator(const Mesh& mesh,
  const Expr& integral,
  const Expr& bcs,
  const Expr& vars,
  const Expr& varEvalPts,
  const Expr& fields,
  const Expr& fieldValues,
  const VectorType<double>& vectorType)
  : assembler_(),
    varValues_(varEvalPts),
    vecType_(vectorType),
    gradient_(1)
{
  Array<Expr> f = tuple(fields.flatten());
  Array<Expr> f0 = tuple(fieldValues.flatten());
  Array<Expr> v = tuple(vars.flatten());
  Array<Expr> v0 = tuple(varEvalPts.flatten());

  Expr params;
  
  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(integral, bcs, v, v0, params, params, f, f0));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vectorType), tuple(vectorType), false));
}

double FunctionalEvaluator::evaluate() const
{
  double value;
  assembler_->evaluate(value);
  return value;
}


Vector<double> FunctionalEvaluator::evalGradientVector(double& value) const 
{
  assembler_->evaluate(value, gradient_);

  return gradient_[0];
}

Expr FunctionalEvaluator::evalGradient(double& value) const 
{
  Vector<double> g = evalGradientVector(value);

  Array<Expr> rtn(assembler_->rowSpace().size());
  for (int i=0; i<rtn.size(); i++)
  {
    std::string name = "gradient";
    if (rtn.size() > 1) name += "[" + Teuchos::toString(i) + "]";
    rtn[i] = new DiscreteFunction(*(assembler_->rowSpace()[i]),
      g.getBlock(i), name);
  }
  if ((int)rtn.size()==1)
  {
    return rtn[0];
  }
  else
  {
    return new ListExpr(rtn);
  }
}


double FunctionalEvaluator::fdGradientCheck(double h) const
{
  bool showAll = false;
  Tabs tabs;
  double f0, fPlus, fMinus;
  Expr gradF0 = evalGradient(f0);

  FancyOStream& os = Out::os();


  DiscreteFunction* df = DiscreteFunction::discFunc(varValues_);
  DiscreteFunction* dg = DiscreteFunction::discFunc(gradF0);
  Vector<double> x = df->getVector();
  Vector<double> x0 = x.copy();
  Vector<double> gf = dg->getVector();

  RCP<GhostView<double> > xView = df->ghostView();
  RCP<GhostView<double> > gradF0View = dg->ghostView();


  TEUCHOS_TEST_FOR_EXCEPTION(xView.get() == 0, std::runtime_error, 
    "bad pointer in FunctionalEvaluator::fdGradientCheck");
  TEUCHOS_TEST_FOR_EXCEPTION(gradF0View.get() == 0, std::runtime_error, 
    "bad pointer in FunctionalEvaluator::fdGradientCheck");

  int nTot = x.space().dim();
  int n = x.space().numLocalElements();
  int lowestIndex = x.space().baseGlobalNaturalIndex();

  os << tabs << "doing fd check:  h=" << h << std::endl;
  Array<double> df_dx(n);

  int localIndex = 0;
  for (int globalIndex=0; globalIndex<nTot; globalIndex++)
  {
    double tmp=0.0;
    bool isLocal = globalIndex >= lowestIndex 
      && globalIndex < (lowestIndex+n);
    if (isLocal)
    {
      tmp = xView->getElement(globalIndex);
      loadable(x)->setElement(globalIndex, tmp + h);
    }

    df->setVector(x);
    fPlus = evaluate();
    if (isLocal)
    {
      loadable(x)->setElement(globalIndex, tmp - h);
    }

    df->setVector(x);
    fMinus = evaluate();
      
    if (isLocal)
    {
      df_dx[localIndex] = (fPlus - fMinus)/2.0/h;
      os << "g=" << setw(5) << globalIndex << ", l=" << setw(5) << localIndex << " f0="
         << setw(12) << f0 
         << " fPlus=" << setw(12) << fPlus << " fMinus=" << setw(12) << fMinus << " df_dx="
         << setw(12) << df_dx[localIndex] << std::endl;
      if (showAll)
      {
        os << "i " << globalIndex << " x_i=" << tmp 
           << " f(x)=" << f0 
           << " f(x+h)=" << fPlus 
           << " f(x-h)=" << fMinus << std::endl;
      }
      loadable(x)->setElement(globalIndex, tmp);
      localIndex++;
    }
    df->setVector(x);
  }
  
  double localMaxErr = 0.0;

  showAll = true;
  VectorSpace<double> space = x.space();
  for (int i=0; i<space.numLocalElements(); i++)
  {
    double num =  fabs(df_dx[i]-gf[i]);
    double den = fabs(df_dx[i]) + fabs(gf[i]) + 1.0e-14;
    double r = 0.0;
    if (fabs(den) > 1.0e-16) r = num/den;
    else r = 1.0;
    if (showAll)
    {
      os << "i " << i;
      os << " FD=" << df_dx[i] 
         << " grad=" << gf[i]
         << " r=" << r << std::endl;
    }
    if (localMaxErr < r) localMaxErr = r;
  }
  os << "local max error = " << localMaxErr << std::endl;
  
  double maxErr = localMaxErr;
  df->mesh().comm().allReduce((void*) &localMaxErr, (void*) &maxErr, 1, 
    MPIDataType::doubleType(), MPIOp::maxOp());
  os << tabs << "fd check: max error = " << maxErr << std::endl;

  return maxErr;
}

