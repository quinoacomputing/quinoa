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

#include "SundanceSpectralPreprocessor.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSumExpr.hpp"
#include "SundanceDiffOp.hpp"
#include "SundanceProductExpr.hpp"
#include "SundanceUnaryMinus.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceFuncElementBase.hpp"
#include "SundanceNonlinearUnaryOp.hpp"
#include "SundanceSpectralExpr.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_Time.hpp"


using namespace std;
using namespace Sundance;
using namespace Teuchos;

static Time& spectralExpansionTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("spectral expansion"); 
  return *rtn;
}


Expr SpectralPreprocessor::projectSpectral(const Expr& e)
{
  if (hasSpectral(e))
  {
    Array<Array<Expr> > terms;
    expandSpectral(e, terms);
    return projectSpectral(terms);
  }
  else
  {
    return e;
  }
}

Expr SpectralPreprocessor::projectSpectral(
  const Array<Array<Expr> >& terms)
{
  Expr rtn = 0.0;
  
  for (int t=0; t<terms.size(); t++)
  {
    Expr test;
    Expr deterministic;
    Array<Expr> specs;
    Array<Expr> testCoeffs;
    Array<Array<Expr> > specCoeffs;
    Array<SpectralBasis> specBas;

    parseProduct(terms[t], test, deterministic, specs);
    
    const SpectralExpr* specTest 
      = dynamic_cast<const SpectralExpr*>(test.ptr().get());
    TEUCHOS_TEST_FOR_EXCEPT(specTest==0);

    SpectralBasis testBasis = specTest->getSpectralBasis();
    for (int i=0; i<testBasis.nterms(); i++)
    {
      testCoeffs.append(specTest->getCoeff(i));
    } 

    for (int i=0; i<specs.size(); i++)
    {
      const SpectralExpr* s 
        = dynamic_cast<const SpectralExpr*>(specs[i].ptr().get());
      SpectralBasis bas = s->getSpectralBasis();
      specBas.append(bas);
      Array<Expr> c;
      for (int j=0; j<bas.nterms(); j++)
      {
        c.append(s->getCoeff(j));
      }
      specCoeffs.append(c);
    }

    if (specs.size()==0)
    {
      for (int i=0; i<testBasis.nterms(); i++)
      {
        rtn = rtn + testCoeffs[i]*deterministic;
      } 
    }
    else if (specs.size()==1)
    {
      for (int i=0; i<testBasis.nterms(); i++)
      {
        rtn = rtn + testCoeffs[i]*specCoeffs[0][i]*deterministic;
      } 
    }
    else if (specs.size()==2U)
    {
      for (int i=0; i<testBasis.nterms(); i++)
      {
        for (int j=0; j<specBas[0].nterms(); j++)
        {
          for (int k=0; k<specBas[1].nterms(); k++)
          {
            double w = testBasis.expectation(
              specBas[0].getElement(j), specBas[1].getElement(k), i) 
              / testBasis.expectation(0,i,i);
            if (w==0.0) continue;
            rtn = rtn + w*testCoeffs[i]
              *specCoeffs[0][j]*specCoeffs[1][k]*deterministic;
          }
        }
      } 
    }
    else
    {
      TEUCHOS_TEST_FOR_EXCEPTION(specs.size() > 2, std::runtime_error,
        "not ready to work with more than two spectral unks");
    }
  }

  return rtn;
}


void SpectralPreprocessor::parseProduct(const Array<Expr>& factors,
  Expr& test, Expr& deterministic, Array<Expr>& spec)
{
  bool hasTest = false;
  Expr det = 1.0;

  for (int i=0; i<factors.size(); i++)
  {
    const Expr& f = factors[i];
    if (isSpectralTest(f)) 
    {
      TEUCHOS_TEST_FOR_EXCEPTION(hasTest, std::runtime_error,
        "multiple test functions detected in product=" << factors);
      test = f;
      hasTest = true;
    }
    else if (isSpectral(f))
    {
      spec.append(f);
    }
    else 
    {
      det = det * f;
    }
  }

  deterministic = det;
}


bool SpectralPreprocessor::isSpectralTest(const Expr& f)
{
  const SpectralExpr* s 
    = dynamic_cast<const SpectralExpr*>(f.ptr().get());
  if (s)
    return s->hasTestFunctions();
  return false;
}

bool SpectralPreprocessor::isSpectral(const Expr& f)
{
  const SpectralExpr* s 
    = dynamic_cast<const SpectralExpr*>(f.ptr().get());
  return s!=0;
}


void SpectralPreprocessor::expandSpectral(const Expr& e, 
  Array<Array<Expr> >& terms)
{
  TimeMonitor timer(spectralExpansionTimer());

  if (!hasSpectral(e)) 
  {
    Array<Expr> z(1);
    z[0] = e;
    terms.append(z);
    return;
  }

  
  const NonlinearUnaryOp* nl 
    = dynamic_cast<const NonlinearUnaryOp*>(e.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(nl!=0, std::runtime_error,
    "expandSpectral() called on a nonlinear unary operator");

  const SumExpr* sum 
    = dynamic_cast<const SumExpr*>(e.ptr().get());

  if (sum) 
  {
    expandSpectralSum(sum, terms);
    return;
  }


  const ProductExpr* prod 
    = dynamic_cast<const ProductExpr*>(e.ptr().get());

  if (prod) 
  {
    expandSpectralProduct(prod, terms);
    return;
  }


  const UnaryMinus* u 
    = dynamic_cast<const UnaryMinus*>(e.ptr().get());

  if (u) 
  {
    expandSpectralUnaryMinus(u, terms);
    return;
  }


  const DiffOp* d 
    = dynamic_cast<const DiffOp*>(e.ptr().get());

  if (d) 
  {
    expandSpectralDiffOp(d, terms);
    return;
  }

  const FuncElementBase* f
    = dynamic_cast<const FuncElementBase*>(e.ptr().get());
  if (f)
  {
    Array<Expr> z(1);
    z[0] = e;
    terms.append(z);
    return;
  }

  const SpectralExpr* se 
    = dynamic_cast<const SpectralExpr*>(e.ptr().get());
  if (se)
  {
    Array<Expr> z(1);
    z[0] = e;
    terms.append(z);
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPT(true);
}


bool SpectralPreprocessor::hasSpectral(const Expr& e)
{
  const SpectralExpr* se 
    = dynamic_cast<const SpectralExpr*>(e.ptr().get());

  if (se) return true;

  const ExprWithChildren* ewc
    = dynamic_cast<const ExprWithChildren*>(e.ptr().get());

  if (ewc)
  {
    for (int i=0; i<ewc->numChildren(); i++)
    {
      if (hasSpectral(ewc->child(i))) return true;
    }
  }
  
  return false;
}

void SpectralPreprocessor::expandSpectralProduct(const ProductExpr* prod,
  Array<Array<Expr> >& terms)
{
  Array<Array<Expr> > L;
  Array<Array<Expr> > R;

  expandSpectral(prod->left(), L);
  expandSpectral(prod->right(), R);

  for (int i=0; i<L.size(); i++)
  {
    for (int j=0; j<R.size(); j++)
    {
      Array<Expr> t;
      for (int p=0; p<L[i].size(); p++) t.append(L[i][p]);
      for (int q=0; q<R[j].size(); q++) t.append(R[j][q]);
      terms.append(t);
    }
  }
}


void SpectralPreprocessor::expandSpectralSum(const SumExpr* sum,
  Array<Array<Expr> >& terms)
{
  Array<Array<Expr> > L;
  Array<Array<Expr> > R;

  expandSpectral(sum->left(), L);
  expandSpectral(sum->sign()*(sum->right()), R);

  for (int i=0; i<L.size(); i++)
  {
    terms.append(L[i]);
  }
  for (int i=0; i<R.size(); i++)
  {
    terms.append(R[i]);
  }
}

void SpectralPreprocessor::expandSpectralUnaryMinus(const UnaryMinus* u,
  Array<Array<Expr> >& terms)
{
  Array<Array<Expr> > A;
  expandSpectral(u->arg(), A);

  for (int i=0; i<A.size(); i++)
  {
    Array<Expr> t;
    t.append(Expr(-1.0));
    for (int j=0; j<A[i].size(); j++)
    {
      t.append(A[i][j]);
    }
    terms.append(t);
  }
}

void SpectralPreprocessor::expandSpectralDiffOp(const DiffOp* d,
  Array<Array<Expr> >& terms)
{
  Array<Array<Expr> > A;
  expandSpectral(d->arg(), A);

  for (int i=0; i<A.size(); i++)
  {
    Array<Array<Expr> > dA;
    expandDerivative(d->mi(), A[i], dA);
    for (int j=0; j<dA.size(); j++) terms.append(dA[j]);
  }
}


void SpectralPreprocessor::expandDerivative(const MultiIndex& mi,
  const Array<Expr>& factors,
  Array<Array<Expr> >& productRuleTerms)
{
  for (int i=0; i<factors.size(); i++)
  {
    Array<Expr> df;
    bool isZero = false;
    for (int j=0; j<factors.size(); j++) 
    {
      if (j!=i) df.append(factors[j]);
      else 
      {
        Expr g = takeDeriv(factors[j], mi);
        df.append(g);
        if (isZeroExpr(g)) 
        {
          isZero = true;
          break;
        }
      }
    }
    if (!isZero) productRuleTerms.append(df);
  }
}


Expr SpectralPreprocessor::takeDeriv(const Expr& f, const MultiIndex& mi)
{
  TEUCHOS_TEST_FOR_EXCEPT(mi.order() > 1);

  
  const SpectralExpr* se 
    = dynamic_cast<const SpectralExpr*>(f.ptr().get());

  Expr d = new Derivative(mi.firstOrderDirection());
  int n = se->getSpectralBasis().nterms();
  
  if (se)
  {
    Array<Expr> c(n);
    for (int i=0; i<n; i++)
    {
      c[i] = d*se->getCoeff(i);
    }
    return new SpectralExpr(se->getSpectralBasis(), c);
  }
  else
  {
    return d*f;
  }
}

bool SpectralPreprocessor::isZeroExpr(const Expr& f)
{
  const ZeroExpr* z = dynamic_cast<const ZeroExpr*>(f.ptr().get());
  return z != 0;
}








