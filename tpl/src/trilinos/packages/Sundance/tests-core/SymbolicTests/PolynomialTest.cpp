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


#include "SundanceExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceParameter.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceFunctionalPolynomial.hpp"

#include "SundanceSymbolicTransformation.hpp"
#include "SundanceProductTransformation.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEquationSet.hpp"

using Sundance::List;
using namespace Sundance;
using namespace Teuchos;


RCP<ScalarExpr> expr2scalar(const Expr& e)
{
  return rcp_dynamic_cast<ScalarExpr>(e[0].ptr());
}

RCP<FunctionalPolynomial> expr2poly(const Expr& e)
{
  return FunctionalPolynomial::toPoly(expr2scalar(e));
}

int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      Expr::showAllParens() = true; 
      ProductTransformation::optimizeFunctionDiffOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

      Expr u = new UnknownFunctionStub("u");
      Expr v = new UnknownFunctionStub("v");
      Expr w = new UnknownFunctionStub("w");
      
      RCP<FunctionalPolynomial> p 
        = rcp(new FunctionalPolynomial(expr2scalar(u)));
      
      RCP<FunctionalPolynomial> q 
        = rcp(new FunctionalPolynomial(expr2scalar(w)));
      
      RCP<FunctionalPolynomial> r
        = rcp(new FunctionalPolynomial(expr2scalar(w)));

      p = p->multiplyPoly(expr2poly(v).get());
      p = p->multiplyPoly(expr2poly(v).get());
      p = p->multiplyPoly(expr2poly(u).get());

      q = q->multiplyPoly(expr2poly(u).get());
      q = q->multiplyPoly(expr2poly(v).get());

      p = p->addPoly(q.get(), 1);

      p = p->addPoly(r.get(), 1);

          
      Out::os() << "p = " << p->toXML() << std::endl;
      p->toText(Out::os(), false) << std::endl;

      TimeMonitor::summarize();      
    }
	catch(std::exception& e)
		{
      Out::os() << e.what() << std::endl;
		}

  
}
