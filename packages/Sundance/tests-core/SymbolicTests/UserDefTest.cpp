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
#include "SundanceUserDefOp.hpp"
#include "SundancePointwiseUserDefFunctor.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceProductTransformation.hpp"
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
using namespace Sundance;
using Sundance::List;
using namespace Sundance;
using namespace Teuchos;


static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}


#define LOUD()                                          \
  {                                                     \
    verbosity<EvaluationTester>() = 5;        \
    verbosity<Evaluator>() = 5;               \
    verbosity<SparsitySuperset>() = 0;         \
    verbosity<EvalVector>() = 0;               \
    verbosity<EvaluatableExpr>() = 5;         \
    verbosity<AbstractEvalMediator>() = 5;    \
  }

#define QUIET()                                         \
  {                                                     \
    verbosity<EvaluationTester>() = 0;         \
    verbosity<Evaluator>() = 0;                \
    verbosity<SparsitySuperset>() = 0;         \
    verbosity<EvalVector>() = 0;               \
    verbosity<EvaluatableExpr>() = 0;          \
    verbosity<AbstractEvalMediator>() = 0;     \
  }


class MyFunc : public PointwiseUserDefFunctor1
{
public:
  MyFunc() : PointwiseUserDefFunctor1("MyFunc", 3, 1){;}
  virtual ~MyFunc(){;}
  void eval1(const double* vars, double* f, double* df) const ;
  void eval0(const double* vars, double* f) const ;
};


void MyFunc::eval1(const double* vars, double* f, double* df) const
{
  f[0] = vars[0]*sin(vars[1]) + 2.0*vars[2]*vars[0];
  df[0] = sin(vars[1]) + 2.0*vars[2];
  df[1] = vars[0]*cos(vars[1]);
  df[2] = 2.0*vars[0];
}

void MyFunc::eval0(const double* vars, double* f) const
{
  f[0] = vars[0]*sin(vars[1]) + 2.0*vars[2]*vars[0];
}



class MyFunc2 : public PointwiseUserDefFunctor1
{
public:
  MyFunc2() : PointwiseUserDefFunctor1("MyFunc2", 3, 2){;}
  virtual ~MyFunc2(){;}
  void eval1(const double* vars, double* f, double* df) const ;
  void eval0(const double* vars, double* f) const ;
};


void MyFunc2::eval1(const double* vars, double* f, double* df) const
{
  f[0] = vars[0]*sin(vars[1]) + 2.0*vars[2]*vars[0];
  f[1] = vars[0]*vars[1]*vars[2];
  df[0] = sin(vars[1]) + 2.0*vars[2];
  df[1] = vars[0]*cos(vars[1]);
  df[2] = 2.0*vars[0];
  df[3] = vars[1]*vars[2];
  df[4] = vars[0]*vars[2];
  df[5] = vars[0]*vars[1];
}

void MyFunc2::eval0(const double* vars, double* f) const
{
  f[0] = vars[0]*sin(vars[1]) + 2.0*vars[2]*vars[0];
  f[1] = vars[0]*vars[1]*vars[2];
}




class MyFunc3 : public PointwiseUserDefFunctor2
{
public:
  MyFunc3() : PointwiseUserDefFunctor2("MyFunc3", 3, 2){;}
  virtual ~MyFunc3(){;}
  void eval2(const double* vars, double* f, double* df, double* d2f) const ;
  void eval1(const double* vars, double* f, double* df) const ;
  void eval0(const double* vars, double* f) const ;
};


void MyFunc3::eval2(const double* vars, double* f, double* df,
                    double* d2f) const
{
  f[0] = vars[0]*sin(vars[1]) + 2.0*vars[2]*vars[0];
  f[1] = vars[0]*vars[1]*vars[2];

  df[0] = sin(vars[1]) + 2.0*vars[2];
  df[1] = vars[0]*cos(vars[1]);
  df[2] = 2.0*vars[0];

  df[3] = vars[1]*vars[2];
  df[4] = vars[0]*vars[2];
  df[5] = vars[0]*vars[1];

  d2f[0] = 0.0;                          // (0,0)
  d2f[1] = cos(vars[1]);                 // (1,0)  
  d2f[2] = -vars[0]*sin(vars[1]);        // (1,1)
  d2f[3] = 2.0;                          // (2,0)
  d2f[4] = 0.0;                          // (2,1)
  d2f[5] = 0.0;                          // (2,2)

  d2f[7] = 0.0;                          // (0,0)
  d2f[8] = vars[2];                      // (1,0)  
  d2f[9] = 0.0;                          // (1,1)
  d2f[10] = 0.0;                         // (2,0)
  d2f[11] = vars[0];                     // (2,1)
  d2f[12] = 0.0;                         // (2,2)
}

void MyFunc3::eval1(const double* vars, double* f, double* df) const
{
  f[0] = vars[0]*sin(vars[1]) + 2.0*vars[2]*vars[0];
  f[1] = vars[0]*vars[1]*vars[2];
  df[0] = sin(vars[1]) + 2.0*vars[2];
  df[1] = vars[0]*cos(vars[1]);
  df[2] = 2.0*vars[0];
  df[3] = vars[1]*vars[2];
  df[4] = vars[0]*vars[2];
  df[5] = vars[0]*vars[1];
}

void MyFunc3::eval0(const double* vars, double* f) const
{
  f[0] = vars[0]*sin(vars[1]) + 2.0*vars[2]*vars[0];
  f[1] = vars[0]*vars[1]*vars[2];
}


int main(int argc, char** argv)
{
  int stat = 0;
  try
		{
      GlobalMPISession session(&argc, &argv);
      Tabs tabs;
      TimeMonitor timer(totalTimer());

      Expr::showAllParens() = true;

      EvalVector::shadowOps() = true;

      ProductTransformation::optimizeFunctionDiffOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);

      ADField U(ADBasis(1), sqrt(2.0));
      ADField W(ADBasis(2), sqrt(3.0));

      ADCoord X(0);
      ADCoord Y(1);
      ADCoord Z(2);

      ADDerivative Dx(0);
      ADDerivative Dy(1);
      ADDerivative Dz(2);

      ADReal C_old = sin(X)*sin(Y);

			Expr u = new TestUnknownFunction(U, "u");
			Expr w = new TestUnknownFunction(W, "w");

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      Expr p = 2.0*(dx*u) + w;
      Expr q = 2.0*u*w + x*w;
      Expr f = new UserDefOp(List(p, q, x), rcp(new MyFunc()));
      Expr F2 = new UserDefOp(List(p, q, x), rcp(new MyFunc2()));


      ADReal P = 2.0*(Dx*U) + W;
      ADReal Q = 2.0*U*W + X*W;
      ADReal adF = P*sin(Q) + 2.0*P*X;
      ADReal adF2 = (P*sin(Q) + 2.0*P*X)*X + (P*Q*X)*Y;
      


      Array<double> df1(2);
      Array<double> df2(2);
      Array<double> df5(2);
      Array<double> df6(2);

      cout << "============== Function #2: R3 --> R1 =======================" << std::endl;
      cout << "evaluating symb expr, eval1()" << std::endl;
      EvaluationTester fTester2(p*sin(q) + 2.0*p*x, 1);
      double f2 = fTester2.evaluate(df2);
      cout << "value (symb)    = " << f2 << std::endl;
      cout << "deriv (symb)    = " << df2 << std::endl;

      //verbosity<Evaluator>() = 5;
      cout << "evaluating user def expr, eval1()" << std::endl;

      EvaluationTester fTester1(f, 1);
      double f1 = fTester1.evaluate(df1);

      cout << "value (functor) = " << f1 << std::endl;
      cout << "deriv (functor) = " << df1 << std::endl;

      cout << "evaluating symb expr, eval0()" << std::endl;
      EvaluationTester fTester3(p*sin(q) + 2.0*p*x, 0);
      double f3 = fTester3.evaluate();
      cout << "value (symb)    = " << f3 << std::endl;

      cout << "value (AD)      = " << adF.value() << std::endl;


      cout << "evaluating user def expr, eval0()" << std::endl;
      EvaluationTester fTester4(f, 0);
      double f4 = fTester4.evaluate();
      cout << "value (functor)    = " << f4 << std::endl;

      cout << "============== Function #2: R3 --> R2 =======================" << std::endl;
      cout << "evaluating symb expr, eval1()" << std::endl;
      EvaluationTester fTester5(x*(p*sin(q) + 2.0*p*x)+y*(p*q*x), 1);
      double f5 = fTester5.evaluate(df5);
      cout << "value (symb)    = " << f5 << std::endl;
      cout << "deriv (symb)    = " << df5 << std::endl;

      //verbosity<Evaluator>() = 5;
      cout << "evaluating user def expr, eval1()" << std::endl;

      EvaluationTester fTester6(F2*List(x,y), 1);
      double f6 = fTester6.evaluate(df6);

      cout << "value (functor) = " << f6 << std::endl;
      cout << "deriv (functor) = " << df6 << std::endl;

      cout << "evaluating symb expr, eval0()" << std::endl;
      EvaluationTester fTester7((p*sin(q) + 2.0*p*x)*x + y*(p*q*x), 0);
      double f7 = fTester7.evaluate();
      cout << "value (symb)    = " << f7 << std::endl;


      cout << "evaluating user def expr, eval0()" << std::endl;
      EvaluationTester fTester8(F2*List(x,y), 0);
      double f8 = fTester8.evaluate();
      cout << "value (functor)    = " << f8 << std::endl;




      

      double tol = 1.0e-10;

      bool isOK = ::fabs(f2 - f1) < tol;
      isOK = ::fabs(f3-f4) < tol && isOK;
      isOK = ::fabs(f5-f6) < tol && isOK;
      isOK = ::fabs(f7-f8) < tol && isOK;
      for (int i=0; i<df1.size(); i++)
        {
          isOK = isOK && ::fabs(df2[i]-df1[i]) < tol;
          isOK = isOK && ::fabs(df5[i]-df6[i]) < tol;
        }


      if (isOK)
        {
          std::cerr << "all tests PASSED!" << std::endl;
        }
      else
        {
          stat = -1;
          std::cerr << "overall test FAILED!" << std::endl;
        }
      TimeMonitor::summarize();
    }
	catch(std::exception& e)
		{
      stat = -1;
      std::cerr << "overall test FAILED!" << std::endl;
      std::cerr << "detected exception: " << e.what() << std::endl;
		}

  return stat;
}
