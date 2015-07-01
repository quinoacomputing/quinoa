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

#include "SundanceUnaryFunctor.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Teuchos;

void UnaryFunctor::eval1(const double* const x, 
  int nx, 
  double* f, 
  double* df_dx) const
{
  evalFDDerivs1(x, nx, f, df_dx);
}


void UnaryFunctor::eval2(const double* const x, 
  int nx, 
  double* f, 
  double* df_dx,
  double* d2f_dxx) const
{
  evalFDDerivs2(x, nx, f, df_dx, d2f_dxx);
}

void UnaryFunctor::eval3(const double* const x, 
  int nx, 
  double* f, 
  double* df_dx,
  double* d2f_dxx,
  double* d3f_dxxx) const
{
  evalFDDerivs3(x, nx, f, df_dx, d2f_dxx, d3f_dxxx);
}


void UnaryFunctor::evalFDDerivs1(const double* const x, 
  int nx, 
  double* f, 
  double* df_dx) const
{
  eval0(x, nx, f);

  double w1 = 1.0/h_/12.0;

  for (int i=0; i<nx; i++)
  {
    double xPlus1 = x[i] + h_;
    double xMinus1 = x[i] - h_;
    double fPlus1;
    double fMinus1;
    double xPlus2 = x[i] + 2.0*h_;
    double xMinus2 = x[i] - 2.0*h_;
    double fPlus2;
    double fMinus2;
    eval0(&xPlus1, 1, &fPlus1);
    eval0(&xMinus1, 1, &fMinus1);
    eval0(&xPlus2, 1, &fPlus2);
    eval0(&xMinus2, 1, &fMinus2);
    df_dx[i] = w1*( 8.0*(fPlus1-fMinus1) - (fPlus2 - fMinus2) );
  }
}

void UnaryFunctor::evalFDDerivs2(const double* const x, 
  int nx, 
  double* f, 
  double* df_dx,
  double* d2f_dxx) const
{
  eval0(x, nx, f);

  double w1 = 1.0/h_/12.0;
  double w2 = w1/h_;

  for (int i=0; i<nx; i++)
  {
    double xPlus1 = x[i] + h_;
    double xMinus1 = x[i] - h_;
    double fPlus1;
    double fMinus1;
    double xPlus2 = x[i] + 2.0*h_;
    double xMinus2 = x[i] - 2.0*h_;
    double fPlus2;
    double fMinus2;
    eval0(&xPlus1, 1, &fPlus1);
    eval0(&xMinus1, 1, &fMinus1);
    eval0(&xPlus2, 1, &fPlus2);
    eval0(&xMinus2, 1, &fMinus2);
    df_dx[i] = w1*( 8.0*(fPlus1-fMinus1) - (fPlus2 - fMinus2) );
    d2f_dxx[i] = w2*(16.0*(fPlus1 + fMinus1) 
      - 30.0*f[i] - (fPlus2+fMinus2));
  }
}

void UnaryFunctor::evalFDDerivs3(const double* const x, 
  int nx, 
  double* f, 
  double* df_dx,
  double* d2f_dxx,
  double* d3f_dxxx) const
{
  eval0(x, nx, f);

  double w1 = 1.0/h_/12.0;
  double w2 = w1/h_;
  double w3 = 0.5/h_/h_/h_;

  for (int i=0; i<nx; i++)
  {
    double xPlus1 = x[i] + h_;
    double xMinus1 = x[i] - h_;
    double fPlus1;
    double fMinus1;
    double xPlus2 = x[i] + 2.0*h_;
    double xMinus2 = x[i] - 2.0*h_;
    double fPlus2;
    double fMinus2;
    eval0(&xPlus1, 1, &fPlus1);
    eval0(&xMinus1, 1, &fMinus1);
    eval0(&xPlus2, 1, &fPlus2);
    eval0(&xMinus2, 1, &fMinus2);
    df_dx[i] = w1*( 8.0*(fPlus1-fMinus1) - (fPlus2 - fMinus2) );
    d2f_dxx[i] = w2*(16.0*(fPlus1 + fMinus1) 
      - 30.0*f[i] - (fPlus2+fMinus2));
    d3f_dxxx[i] = w3*(fPlus2 - 2.0*fPlus1 + 2.0*fMinus1 - fMinus2);
  }
}

bool UnaryFunctor::testDerivs(const double& x, const double& tol) const
{
  Tabs tabs;
  Out::os() << tabs << std::endl << tabs 
            << "comparing exact derivs to FD derivs for functor " 
            << name() << std::endl;

  double fExact;
  double fxExact;
  double fxxExact;

  double fFD;
  double fxFD;
  double fxxFD;

  bool isOK = true;

  /* test first differentiation */
  Out::os() << tabs << "computing first derivatives at x=" << x << std::endl;
  {
    Tabs tabs1;

    eval1(&x, 1, &fExact, &fxExact);
    Out::os() << tabs1 << "Exact: f=" << fExact << " df_dx=" << fxExact << std::endl; 
    evalFDDerivs1(&x, 1, &fFD, &fxFD);
    Out::os() << tabs1 << "FD:    f=" << fFD << " df_dx=" << fxFD << std::endl; 

    double fError = fabs(fFD - fExact)/(fabs(fExact) + h_);
    double fxError = fabs(fxFD - fxExact)/(fabs(fxExact)+h_);
    {
      Tabs tabs2;
      Out::os() << tabs2 << "| f-f_FD |=" << fError 
                << " | df-df_FD |=" << fxError << std::endl;
      if (fError > tol) 
      {
        Out::os() << tabs << "ERROR: function value mismatch!" << std::endl;
        isOK = false;
      }
      if (fxError > tol) 
      {
        Out::os() << tabs << "ERROR: first derivative mismatch!" << std::endl;
        isOK = false;
      }
    }
  }
  
  Out::os() << tabs << std::endl;

  
  /* test first and second differentiation */
  Out::os() << tabs << "computing first and second derivatives at x=" << x << std::endl;
  {
    Tabs tabs1;

    eval2(&x, 1, &fExact, &fxExact, &fxxExact);
    Out::os() << tabs1 << "Exact: f=" << fExact 
              << " df_dx=" << fxExact 
              << " d2f_dx2=" << fxxExact 
              << std::endl; 

    evalFDDerivs2(&x, 1, &fFD, &fxFD, &fxxFD);
    Out::os() << tabs1 << "FD:    f=" << fFD << " df_dx=" << fxFD  
              << " d2f_dx2=" << fxxFD << std::endl; 

    double fError = fabs(fFD - fExact)/(fabs(fExact) + h_);
    double fxError = fabs(fxFD - fxExact)/(fabs(fxExact)+h_);
    double fxxError = fabs(fxxFD - fxxExact)/(fabs(fxxExact)+h_);
    {
      Tabs tabs2;
      Out::os() << tabs2 << "| f-f_FD |=" << fError 
                << " | df-df_FD |=" << fxError 
                << " | d2f-d2f_FD |=" << fxxError << std::endl;
      if (fError > tol) 
      {
        Out::os() << tabs << "ERROR: function value mismatch!" << std::endl;
        isOK = false;
      }
      if (fxError > tol) 
      {
        Out::os() << tabs << "ERROR: first derivative mismatch!" << std::endl;
        isOK = false;
      }
      if (fxxError > tol) 
      {
        Out::os() << tabs << "ERROR: second derivative mismatch!" << std::endl;
        isOK = false;
      }
    }
  }

  Out::os() << tabs << std::endl;
  return isOK;
}

bool UnaryFunctor::testInvalidValue(const double& badValue) const 
{
  Tabs tabs;
  bool detectedError = false;

  Out::os() << std::endl << tabs 
            << "testing exception detection for bad value x=" << badValue
            << " for function  " << name() << std::endl;
  try
  {
    testDerivs(badValue, 1.0);
  }
  catch(std::exception& e)
  {
    detectedError = true;
  }

  if (!detectedError) 
  {
    Out::os() << tabs << "ERROR: missed detection of an exception at x=" << badValue
              << " for function " << name() << std::endl;
  }
  else
  {
    Out::os() << tabs << "the error was detected!" << std::endl;
  }
  return detectedError;
}


bool UnaryFunctor::test(int nx, const double& tol) const
{
  bool isOK = true;

  double a = -sqrt(1.5);
  double b = sqrt(2.0);

  /* first test a sample of points in the domain */
  if (domain()->hasLowerBound())
  {
    a = domain()->lowerBound();
  }
  if (domain()->hasUpperBound())
  {
    b = domain()->upperBound();
  }

  double c = (b-a)/((double) nx+1);
  for (int i=1; i<=nx; i++)
  {
    double x = a + i*c;
    Out::os() << "testing at " << x << std::endl;
    if (domain()->hasExcludedPoint() && fabs(x-domain()->excludedPoint())<1.0e-14)
    {
      continue;
    }
    isOK =  testDerivs(x, tol) && isOK ;
  }

  /* Test for exception detection at the excluded point, if any */
  if (domain()->hasExcludedPoint())
  {
    Out::os() << "testing detection of excluded point x=" 
              << domain()->excludedPoint() << std::endl;
    isOK =  testInvalidValue(domain()->excludedPoint()) && isOK ;
  }

  /* Test for exception detection at a point below the lower bound */
  if (domain()->hasLowerBound())
  {
    Out::os() << "testing exception below lower bound x=" 
              << domain()->lowerBound() << std::endl;
    isOK =  testInvalidValue(domain()->lowerBound() - 0.1) && isOK ;
  }
  
  /* Test for exception detection at a point above the upper bound */
  if (domain()->hasUpperBound())
  {
    Out::os() << "testing exception above lower bound x=" 
              << domain()->upperBound() << std::endl;

    isOK =  testInvalidValue(domain()->upperBound() + 0.1) && isOK ;
  }
  
  return isOK;
  

}
