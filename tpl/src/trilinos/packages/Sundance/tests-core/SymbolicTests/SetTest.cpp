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
#include "SundanceStdMathOps.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceUnknownFuncElement.hpp"
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
#include "SundanceEvalVector.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceStringEvalMediator.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

using Sundance::List;


static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}

int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      TimeMonitor t(totalTimer());

      int maxDiffOrder = 0;

      verbosity<SymbolicTransformation>() = 0;
      verbosity<Evaluator>() = 0;
      verbosity<EvalVector>() = 0;
      verbosity<EvaluatableExpr>() = 5;
      Expr::showAllParens() = true;

      EvalVector::shadowOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

			Expr u = new UnknownFunctionStub("u");
			Expr w = new UnknownFunctionStub("w");
			Expr v = new TestFunctionStub("v");
			Expr alpha = new UnknownFunctionStub("alpha");

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      Expr u0 = new DiscreteFunctionStub("u0");
      Expr w0 = new DiscreteFunctionStub("w0");
      Expr alpha0 = new DiscreteFunctionStub("alpha0");
      Expr dum;

      Array<Expr> tests;

      tests.append( dx*u );


      for (int i=0; i<tests.length(); i++)
        {
          for (int m=0; m<=maxDiffOrder; m++)
          {
            Out::os() << "=========  diff order "  << m << "========"<<endl;
            RegionQuadCombo rqc(rcp(new CellFilterStub()), 
              rcp(new QuadratureFamilyStub(1)));
            EvalContext context(rqc, m, EvalContext::nextID());
            DerivSet ds = SymbPreprocessor::setupFwdProblem(
              tests[i],
              v, u, u0,
              dum, dum, dum, dum,
              alpha, alpha0,
              context);
            Out::os() << "derivs=" << ds << std::endl;
          }
        }

      TimeMonitor::summarize();
    }
	catch(std::exception& e)
		{
			Out::println(e.what());
		}

  return 0;
  
}
