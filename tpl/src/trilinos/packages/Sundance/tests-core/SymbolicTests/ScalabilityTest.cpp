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
#include "SundanceEvaluationTester.hpp"

using namespace Sundance;
using namespace SundanceTesting;
using namespace Teuchos;

using std::cout;
using std::exception;

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
      Tabs tabs;
      TimeMonitor timer(totalTimer());

      //       verbosity<SymbolicTransformation>() = 0;
      //       verbosity<EvaluationTester>() = 5;
      //       verbosity<Evaluator>() = 5;
      //       verbosity<EvalVector>() = 5;
      //       verbosity<EvaluatableExpr>() = 5;
      //       verbosity<AbstractEvalMediator>() = 5;
      Expr::showAllParens() = true;

      EvalVector::shadowOps() = false;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);

      ADCoord X(0);
      ADCoord Y(1);
      ADCoord Z(2);

      ADDerivative Dx(0);
      ADDerivative Dy(1);
      ADDerivative Dz(2);

      Time stopwatch("test");

      for (int n=2; n<30; n++)
        {
          double t0 = stopwatch.wallTime();
          stopwatch.start();
          int numFields = n;

          Array<ADField> U(numFields);
          Array<Expr> u(numFields);

          Expr big = 0.0;
          ADReal value;
          for (int i=0; i<numFields; i++)
            {
              U[i] = ADField(ADBasis(i+1), ::sqrt(i+1.0));
              u[i] = new TestUnknownFunction(U[i], 
                                             "u" + Teuchos::toString(i) + "_");
              big = big + u[i]*u[i]*x;
              value = value + U[i]*U[i]*X;
              if (i > 0) 
                {
                  big = big + (dx*u[i-1])*(dy*u[i]) 
                    + (dy*u[i-1])*(dz*u[i])
                    + (dx*u[i-1])*(dz*u[i]);
                  value = value + (Dx*U[i-1])*(Dy*U[i])
                    + (Dy*U[i-1])*(Dz*U[i])
                    + (Dx*U[i-1])*(Dz*U[i]);
                }
            }

          EvaluationTester tester(big);
          Array<double> df;
          Array<Array<double> > df2;
          tester.evaluate(df, df2);
          stopwatch.stop();
          double t1 = stopwatch.wallTime();
          int nNodes = tester.numNodes();
          int nnz = tester.numNonzeros();
          std::cerr << n << "   " << nNodes << "   " 
               << nnz << "    " << t1-t0 << std::endl;
          TimeMonitor::summarize();
        }
    }
	catch(std::exception& e)
		{
			Out::println(e.what());
		}


  
}
