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
#include "SundanceUserDefOp.hpp"
#include "SundancePointwiseUserDefFunctor.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceFunctionalDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceCellVectorExpr.hpp"
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

using namespace SundanceTesting;
using Sundance::List;
using namespace Teuchos;
using std::exception;

/* a^b */
Expr myPow( Expr a , Expr b )
{
  return exp( log( a ) * b );
}

static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}



class MyScalarFunc : public PointwiseUserDefFunctor1
{
public:
  MyScalarFunc() : PointwiseUserDefFunctor1("MyScalarFunc", 3, 1){;}
  virtual ~MyScalarFunc(){;}
  void eval1(const double* vars, double* f, double* df) const ;
  void eval0(const double* vars, double* f) const ;
};




class MyVectorFunc1 : public PointwiseUserDefFunctor1
{
public:
  MyVectorFunc1() : PointwiseUserDefFunctor1("MyVectorFunc1", 3, 2){;}
  virtual ~MyVectorFunc1(){;}
  void eval1(const double* vars, double* f, double* df) const ;
  void eval0(const double* vars, double* f) const ;
};


class MyVectorFunc2 : public PointwiseUserDefFunctor2
{
public:
  MyVectorFunc2() : PointwiseUserDefFunctor2("MyVectorFunc2", 3, 2){;}
  virtual ~MyVectorFunc2(){;}
  void eval2(const double* vars, double* f, double* df, double* d2f) const ;
};

class MyScalarFunc2 : public PointwiseUserDefFunctor2
{
public:
  MyScalarFunc2() : PointwiseUserDefFunctor2("F", 3, 1){;}
  virtual ~MyScalarFunc2(){;}
  void eval2(const double* vars, double* f, double* df, double* d2f) const ;
};


inline Expr FD(const Expr& f, const Expr& u) 
{return FunctionalDerivative(f, u);}



#define LOUD()                                          \
  {                                                     \
    verbosity<EvaluationTester>() = 6;        \
    verbosity<Evaluator>() = 6;               \
    verbosity<SparsitySuperset>() = 6;         \
    verbosity<EvalVector>() = 6;               \
    verbosity<EvaluatableExpr>() = 6;         \
    verbosity<AbstractEvalMediator>() = 6;    \
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




#define TESTER_N(expr, adExpr, order)                                   \
  {                                                                     \
    Tabs tabs1;                                                         \
    Out::os() << tabs1 << std::endl << tabs1                                      \
              << "------------- Testing " << #expr << " -----------"        \
              << std::endl << tabs1 << std::endl;                                      \
    bool thisTestIsOK = true;                                           \
    try                                                                 \
    {                                                                 \
      EvaluationTester tester((expr), order);                               \
      double f = tester.fdEvaluate(fdStep, tol1, tol2, thisTestIsOK);   \
      if (!thisTestIsOK)                                                \
      {                                                               \
        isOK = false;                                                 \
      }                                                               \
      double adf = (adExpr).value();                                    \
      Out::os() << tabs1 << "expr value = " << f << " check=" << adf         \
                << " |f-check|=" << fabs(f-adf) << std::endl;                     \
      double fError = fabs(f-adf);                                      \
      if (fError > tol1)                                                \
      {                                                               \
        thisTestIsOK=false;                                           \
        Out::os() << "value computation FAILED" << std::endl;                   \
        isOK = false;                                                 \
      }                                                               \
    }                                                                 \
    catch(std::exception& ex)                                           \
    {                                                                 \
      thisTestIsOK = false;                                           \
      isOK=false;                                                     \
      Out::os() << "exception: " << ex.what() << std::endl ;                    \
    }                                                                 \
    if (!thisTestIsOK)                                                  \
    {                                                                 \
      failures.append(#expr);                                         \
      Out::os() << "test " << (expr).toString() << " FAILED" << std::endl << std::endl;\
    }                                                                 \
    else                                                                \
    {                                                                 \
      Out::os() << "test " << (expr).toString() << " PASSED" << std::endl << std::endl; \
    }\
  }

#define TESTER(expr, adExpr) TESTER_N(expr, adExpr, 2)
#define TESTER1(expr, adExpr) TESTER_N(expr, adExpr, 1)



#define XTESTER(X,Y)                            \
  LOUD();                                       \
  TESTER(X,Y);                                   \
  QUIET();                                      \
  goto finish;

#define XTESTER1(X,Y)                            \
  LOUD();                                       \
  TESTER1(X,Y);                                   \
  QUIET();                                      \
  goto finish;





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
    ADField V(ADBasis(2), sqrt(2.5));
    ADField W(ADBasis(2), sqrt(3.0));
    ADField B(ADBasis(2), 0.0);
    ADField C(ADBasis(2), 0.0);

    ADCoord X(0);
    ADCoord Y(1);
    ADCoord Z(2);

    ADDerivative Dx(0);
    ADDerivative Dy(1);
    ADDerivative Dz(2);

    ADReal C_old = sin(X)*sin(Y);

    Expr u = new TestUnknownFunction(U, "u");
    Expr v = new TestUnknownFunction(V, "v");
    Expr w = new TestUnknownFunction(W, "w");
    Expr b = new TestUnknownFunction(B, "b");
    Expr c = new TestUnknownFunction(C, "c");

    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);

    Expr g = x*x + y*y;
    Expr f = x*x;
    Expr h = x+y;
    Expr c_old = sin(x)*sin(y);
    double dt = 0.01;

     

    Expr grad = List(dx, dy);

    double tol1 = 1.0e-5;
    double tol2 = 1.0e-5;
    double fdStep = 1.0e-5;
    bool isOK = true;
    Array<string> failures;



#ifdef BLARF
    XTESTER(b, B);
#endif



    TESTER(u, U);



    TESTER(-u, -U);

      
    /* ----------- tests of symbolic simplifications -------------*/
    TESTER( u - u, U - U );

    TESTER( u + u, 2.0*U );



    /* ----------- distinct cases of sum expressions ------------------- */

    /* tests const-const and vec-vec sums */
    TESTER( u + w, U + W );

    /* tests const-const and vec-vec subtractions */
    TESTER( u - w, U - W );

    /* tests vec-vec and const-0 sums */
    TESTER( u + x, U + X );

    /* tests vec-vec and const-0 subtractions */
    TESTER( u - x, U - X );

    /* tests vec-vec and const-0 subtractions */
    TESTER( u - (w + v), U - (W + V) );

    /* tests vec-vec and 0-const sums */
    TESTER( x + u, X + U );

    /* tests vec-vec and 0-const subtractions */
    TESTER( x - u, X - U );

    /* tests vec-vec and const-vec sums */
    TESTER( u + x*u, U + X*U );

    /* tests vec-vec and const-vec subtractions */
    TESTER( u - x*u, U - X*U );

    /* tests vec-vec and vec-const sums */
    TESTER( x*u + u, X*U + U );

    /* tests vec-vec and vec-const subtractions */
    TESTER( x*u - u, X*U - U );

    /* ----------- cases of product expressions ------------------- */

    /* */

    TESTER( u*u, U*U );

    /* */
    TESTER( u*u*u, U*U*U );

    /* */
    TESTER( u*w, U*W );

    /* */
    TESTER( w*(u*u-2.0), W*(U*U-2.0) );

    /* */
    TESTER( u*u*w, U*U*W );

    /* */
    TESTER( w*u*w, W*U*W );

    /* */
    TESTER( w*u*u, W*U*U );

    /* */
    TESTER( (u+w)*u, (U+W)*U );

    /* */
    TESTER( u*(u+w), U*(U+W) );

    /* */
    TESTER( (u+w*u)*u, (U+W*U)*U );

    /* */
    TESTER( u*(u+w*u), U*(U+W*U) );

    /* */
    TESTER( (u+w)*(u+w), (U+W)*(U+W) );

    /* */
    TESTER( (u+x)*(u+w), (U+X)*(U+W) );

    /* */
    TESTER( (u+w*u)*(u+w), (U+W*U)*(U+W) );

    /* */
    TESTER( (u+w)*(u+w*u), (U+W)*(U+W*U) );


    /* */
    TESTER( (x*u)*(y*u), (X*U)*(Y*U) );

    /* */
    TESTER( (x*u)*(y*w), (X*U)*(Y*W) );


    /* */
    TESTER( (2.0*u)*(y*u), (2.0*U)*(Y*U) );

    /* */
    TESTER( (x*u)*(2.0*u), (X*U)*(2.0*U) );

    /* */
    TESTER( (2.0*u)*(4.0*u), (2.0*U)*(4.0*U) );

    /* */
    TESTER( (2.0*u*u)*(y*u), (2.0*U*U)*(Y*U) );

    /* */
    TESTER( (x*u)*(2.0*u*u), (X*U)*(2.0*U*U) );

    /* */
    TESTER( 2.0*(y*u), 2.0*(Y*U) );

    /* */
    TESTER( (x*u)*2.0, (X*U)*2.0 );

    /* */
    TESTER( u*(y*u), U*(Y*U) );

    /* */
    TESTER( (x*u)*u, (X*U)*U );


    /* -------------- tests of diff ops ----------------------- */


    //#endif     

    /* */
    TESTER1((dx*u), (Dx*U));

    /* */
    TESTER((dx*x), (Dx*X));

    TESTER((dx*(x*x)), (Dx*(X*X)));

    TESTER((dx*(y*y)), (Dx*(Y*Y)));

    TESTER((dx*(x*x + y*y)), (Dx*(X*X + Y*Y)));


// disabled by KL
//    TESTER((b*x*dx*(b) + b*b*(dx*x)), (B*(Dx*(X*B))));
//    TESTER((b*b*(dx*x)), (B*B)*(Dx*X));
    TESTER((dx*x), (Dx*X));
    TESTER((dx*(u+x)), (Dx*U + Dx*X));
//    TESTER((b*b), (B*B));
      
    TESTER((c*dx*(x*b)), (C*(Dx*(X*B))));

    TESTER((u*dx*(x*b)), (U*(Dx*(X*B))));
    TESTER((u*dx*(x*b)), (U*(Dx*(X*B))));

// disabled by KL
//    TESTER((dx*(sin(x)*b)), (Dx*(sin(X)*B)));



    /* */
    TESTER((dx*y), (Dx*Y));

    TESTER((dx*u+dx*w), (Dx*U+Dx*W));


      
    TESTER((dx*(u+w)), (Dx*(U+W)));


    TESTER((dx*(u-w)), (Dx*(U-W)));

    TESTER((dx*(x+w)), (Dx*(X+W)));

    TESTER((dx*(x-w)), (Dx*(X-W)));

    TESTER((dx*(x*w)), (Dx*(X*W)));


    TESTER((dx*(x*w+w)), (Dx*(X*W+W)));

    TESTER((dx*(x*w-w)), (Dx*(X*W-W)));

    TESTER((dx*(w + x*w)), (Dx*(W + X*W)));


    TESTER((dx*(w - x*w)), (Dx*(W - X*W)));

    TESTER((dx*(y*w)), (Dx*(Y*W)));

    TESTER((dx*(u*w)), (Dx*(U*W)));

    TESTER((dx*(3.0*u*w)), (Dx*(3.0*U*W)));

    TESTER((u*(dx*u)), (U*(Dx*U)));

    TESTER((dx*(u*u)), (Dx*(U*U)));


    TESTER((dx*(h)), (Dx*(X+Y)));
    TESTER((dx*(x*h)), (Dx*(X*(X+Y))));

    TESTER((dx*(h) + dy*(h)), (Dx*(X+Y) + Dy*(X+Y)));

    TESTER((dx*(x*h) + dy*(y*h)), (Dx*(X*(X+Y)) + Dy*(Y*(X+Y))));

    TESTER((dx*(y*h) + dy*(x*h)), (Dx*(Y*(X+Y)) + Dy*(X*(X+Y))));

    TESTER((u*(dx*(f) + dy*(f))), (U*(Dx*(X*X) + Dy*(X*X))));


    TESTER((u*(dx*(x*x+y*y) + dy*(x*x+y*y))), (U*(Dx*(X*X + Y*Y) + Dy*(X*X + Y*Y))));

    TESTER((u*(dx*(g) + dy*(g))), (U*(Dx*(X*X + Y*Y) + Dy*(X*X + Y*Y))));
       


    TESTER((u*g + w*g), (U*(X*X + Y*Y) + W*(X*X + Y*Y)));



      
    TESTER((dx*(2.0*u+4.0*w)), (Dx*(2.0*U+4.0*W)));

    TESTER((dx*(2.0*u-4.0*w)), (Dx*(2.0*U-4.0*W)));

    TESTER((dx*(2.0*x+4.0*w)), (Dx*(2.0*X+4.0*W)));

    TESTER((dx*(2.0*x-4.0*w)), (Dx*(2.0*X-4.0*W)));



    TESTER((dx*(x*w+u)), (Dx*(X*W+U)));

    TESTER((dx*(w + x*w)), (Dx*(W + X*W)));

    TESTER((dx*(x*w+u*y)), (Dx*(X*W+U*Y))); 

    TESTER(((dx*u)+(dx*w) + w), ((Dx*U)+(Dx*W) + W));

    TESTER(((dx*u)*(dx*w) + 2.0*w), ((Dx*U)*(Dx*W) + 2.0*W));

    TESTER(((dx*u)*(dx*w) + (dy*u)*(dy*w)), ((Dx*U)*(Dx*W) + (Dy*U)*(Dy*W)));


    TESTER((x*(dx*u)*(dx*w)), (X*(Dx*U)*(Dx*W))); 

    TESTER((u*(dx*u)*(dx*w)), (U*(Dx*U)*(Dx*W)));

    TESTER((x*(dx*u)*(dx*w)), (X*(Dx*U)*(Dx*W)));

    TESTER((u*(dx*u)*(dx*w)), (U*(Dx*U)*(Dx*W)));

    TESTER((x*(dx*u)*(dx*u)), (X*(Dx*U)*(Dx*U)));

    TESTER((u*(dx*u)*(dx*u)), (U*(Dx*U)*(Dx*U)));

    TESTER((u*(dx*u)), (U*(Dx*U)));

    TESTER((u*(dx*(u*x))), (U*(Dx*(U*X))));

    /* Unary operators */
    TESTER(sin(u), sin(U));

    TESTER(sin(u*u + w*w + w*u), sin(U*U + W*W + W*U));

    TESTER(sin(0.5*u), sin(0.5*U));

    TESTER(sin(u+w), sin(U+W));

    TESTER(sin(u*w), sin(U*W));

    TESTER(sin(0.5+u), sin(0.5+U));

    TESTER(sin(x), sin(X));

    TESTER(sin(u*x), sin(U*X));

    TESTER(sin(0.5*u*x), sin(0.5*U*X));

    TESTER(w*sin(u), W*sin(U));

    TESTER(w*sin(u-x), W*sin(U-X));

    TESTER(sin(dx*u), sin(Dx*U));

    TESTER(w*sin(dx*u), W*sin(Dx*U));


    TESTER(sin(u*cos(x)), sin(U*cos(X)));

    TESTER(sin(0.5*(w+u)), sin(0.5*(W+U)));

    TESTER(sin(0.5*w*u), sin(0.5*W*U));

    TESTER(sin(u*x+0.5*w*u), sin(U*X+0.5*W*U));

    TESTER(sin(cos(u*x)+0.5*w*u), sin(cos(U*X)+0.5*W*U));


    TESTER(sin(u)/u, sin(U)/U);

    TESTER(sin(cos(u)), sin((cos(U))));

    TESTER(sin(cos(x)), sin((cos(X))));

    TESTER((dx*sin(x)), (Dx*(sin(X))));

    TESTER((dx*sin(u)), (Dx*(sin(U))));

    TESTER((dx*exp(u)), (Dx*(exp(U))));

    TESTER((dx*exp(u+w)), (Dx*(exp(U+W))));

    TESTER((dx*exp(u*u)), (Dx*(exp(U*U))));

    TESTER((dx*exp(u*w)), (Dx*(exp(U*W))));

    TESTER((dx*exp(2.0*u)), (Dx*(exp(2.0*U))));

    TESTER((dx*exp(-u)), (Dx*(exp(-U))));

    TESTER((dx*exp(-1.0*u)), (Dx*(exp(-1.0*U))));

    TESTER((dx*(cos(x)*sin(x))), (Dx*(cos(X)*sin(X))));

    TESTER((dx*(cos(x)*sin(y))), (Dx*(cos(X)*sin(Y))));

    TESTER((dx*sin(cos(u))), (Dx*(sin(cos(U)))));

    TESTER((dx*sin(2.0*cos(u))), (Dx*(sin(2.0*cos(U)))));

    TESTER(w*(dx*sin(u)), W*(Dx*(sin(U))));

    TESTER(dx*(pow(u, 2.0)), Dx*(pow(U, 2.0)));

    TESTER(dx*(pow(u, 4.0)), Dx*(pow(U, 4.0)));

    TESTER(pow(u-x, 2.0), pow(U-X, 2.0));

    TESTER(log(u), log(U));

    TESTER(dx*(log(u)), Dx*(log(U)));

    TESTER(dx*(exp(log(u))), Dx*(exp(log(U))));

    TESTER(dx*(exp(-log(u))), Dx*(exp(-log(U))));

    TESTER(dx*(exp(2.0*log(u))), Dx*(exp(2.0*log(U))));

//    TESTER(w*((u-c_old)/dt) + (grad*w)*(grad*(u + c_old)/2.0),
//      W*((U-C_old)/dt) + (Dx*W)*(Dx*(U + C_old)/2.0) + (Dy*W)*(Dy*(U + C_old)/2.0));
             


    /* ------ user-defined operators ---------- */
    Expr p = 2.0*(dx*u) + w;
    Expr q = 2.0*u*w + x*w;

    Expr r = 3.0*u + x*w;
    Expr s = sin(x)*u + w;

    /* a singly-differentiable scalar-valued user-def function */
    Expr sfA = new UserDefOp(List(p,q,x), rcp(new MyScalarFunc()));
    Expr sfB = new UserDefOp(List(r,s,y), rcp(new MyScalarFunc()));

    /* a singly-differentiable vector-valued user-def function */
    Expr vfA = new UserDefOp(List(p,q,x), rcp(new MyVectorFunc1()));
    Expr vfB = new UserDefOp(List(r,s,y), rcp(new MyVectorFunc1()));

    /* a twice-differentiable vector-valued user-def function */
    Expr vfC = new UserDefOp(List(r,s,y), rcp(new MyVectorFunc2()));

    /* a twice-differentiable scalar-valued user-def function */
    Expr sf2 = new UserDefOp(List(r,s,y), rcp(new MyScalarFunc2()));

    ADReal P = 2.0*(Dx*U) + W;
    ADReal Q = 2.0*U*W + X*W;
    ADReal R = 3.0*U + X*W;
    ADReal S = sin(X)*U + W;
    ADReal adSF_PQ = P*sin(Q) + 2.0*P*X;
    ADReal adSF_RS = R*sin(S) + 2.0*R*Y;
    ADReal adVF_PQ = (P*sin(Q) + 2.0*P*X)*X + (P*Q*X)*Y;
    ADReal adVF_RS = (R*sin(S) + 2.0*R*Y)*X + (R*S*Y)*Y;

    ADReal adSF2 = R*S + S*S + Y*Y;

    TESTER1(sfA, adSF_PQ);

    TESTER1(sfB, adSF_RS);

    TESTER1(vfA*List(x,y), adVF_PQ);

    TESTER1(vfB*List(x,y), adVF_RS);

    TESTER1(vfC*List(x,y), adVF_RS);

    TESTER(sf2, adSF2);

    TESTER( FD(u*u, u), (2.0*U) );

    TESTER( FD(List(u*u, u), u)*List(1.0, 1.0), ((2.0*U)+1.0) );

    TESTER( 
      (List(1.0, 1.0) * (FD(List(sin(u)*v, v*v+u*u*u), List(u, v)) * List(1.0, 1.0))), 
      (V*cos(U) + 3.0*U*U + sin(U) + 2.0*V)
      );

    TESTER( FD(0.5*u*u, u), U );

    TESTER( FD(0.5*u*v, v), 0.5*U );

    TESTER( FD(sin(u), u), cos(U) );

    TESTER( cos(FD(sin(u), u)), cos(cos(U)) );

    TESTER1( dx*FD(0.5*u*u, u), Dx*U );

    Expr n = CellNormalExpr(2, "n");
    double Nx = 0.5;
    double Ny = ::sqrt(3.0)/2.0;

    TESTER( List(u, v)*n, U*Nx + V*Ny );

    TESTER(exp(log(u) * 2.3/(1.0 + 0.0*v/100.0)), 
      exp(log(U) * 2.3/(1.0 + 0.0*V/100.0)));


    goto finish; /* does nothing other than to avoid a label-not-used warning
                  * when debugging is off */

  finish:
    if (isOK)
    {
      Out::os() << "all tests PASSED!" << std::endl;
    }
    else
    {
      stat = -1;
      Out::os() << "overall test FAILED!" << std::endl;
      Out::os() << std::endl << "failed exprs: " << std::endl
                << std::endl;
      for (int i=0; i<failures.size(); i++)
      {
        Out::os() << failures[i] << std::endl;
      }
      Out::os() << std::endl;
    }
    TimeMonitor::summarize();
  }
	catch(std::exception& e)
  {
    stat = -1;
    Out::os() << "overall test FAILED!" << std::endl;
    Out::os() << "detected exception: " << e.what() << std::endl;
  }
  return stat;
}



void MyScalarFunc::eval1(const double* vars, double* f, double* df) const
{
  f[0] = vars[0]*sin(vars[1]) + 2.0*vars[2]*vars[0];
  df[0] = sin(vars[1]) + 2.0*vars[2];
  df[1] = vars[0]*cos(vars[1]);
  df[2] = 2.0*vars[0];
}

void MyScalarFunc::eval0(const double* vars, double* f) const
{
  f[0] = vars[0]*sin(vars[1]) + 2.0*vars[2]*vars[0];
}

void MyVectorFunc1::eval1(const double* vars, double* f, double* df) const
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

void MyVectorFunc1::eval0(const double* vars, double* f) const
{
  f[0] = vars[0]*sin(vars[1]) + 2.0*vars[2]*vars[0];
  f[1] = vars[0]*vars[1]*vars[2];
}

void MyVectorFunc2::eval2(const double* vars, double* f, double* df,
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

  d2f[6] = 0.0;                          // (0,0)
  d2f[7] = vars[2];                      // (1,0)  
  d2f[8] = 0.0;                          // (1,1)
  d2f[9] = vars[1];                     // (2,0)
  d2f[10] = vars[0];                     // (2,1)
  d2f[11] = 0.0;                         // (2,2)
}

void MyScalarFunc2::eval2(const double* vars, double* f, double* df,
  double* d2f) const
{
  f[0] = vars[0]*vars[1] + vars[1]*vars[1] + vars[2]*vars[2];

  df[0] = vars[1];
  df[1] = vars[0] + 2.0*vars[1];
  df[2] = 2.0*vars[2];

  d2f[0] = 0.0;                          // (0,0)
  d2f[1] = 1.0;                          // (1,0)  
  d2f[2] = 2.0;                          // (1,1)
  d2f[3] = 0.0;                          // (2,0)
  d2f[4] = 0.0;                          // (2,1)
  d2f[5] = 2.0;                          // (2,2)
}

