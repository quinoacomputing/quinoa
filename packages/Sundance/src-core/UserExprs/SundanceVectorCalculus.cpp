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

#include "SundanceVectorCalculus.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceListExpr.hpp"


namespace Sundance
{

Expr gradient(int dim)
{
  Array<Expr> rtn(dim);
  for (int i=0; i<dim; i++)
  {
    rtn[i] = new Derivative(i);
  }
  return new ListExpr(rtn);
}

Expr div(const Expr& f)
{
  Expr rtn = 0.0;
  for (int i=0; i<f.size(); i++)
  {
    Expr d = new Derivative(i);
    rtn = rtn + d*f[i];
  }
  return rtn;
}


Expr cross(const Expr& a, const Expr& b)
{
  TEUCHOS_TEST_FOR_EXCEPTION(a.size() != b.size(), std::runtime_error,
    "mismatched vector sizes in cross(a,b): a.size()=" << a.size()
    << ", b.size()=" << b.size());

  TEUCHOS_TEST_FOR_EXCEPTION(a.size() < 2 || a.size() > 3, std::runtime_error,
    "cross(a,b) undefined for dim=" << a.size());

  if (a.size()==2)
  {
    return a[0]*b[1] - a[1]*b[0];
  }
  else
  {
    return List(
      cross(List(a[1],a[2]), List(b[1],b[2])),
      -cross(List(a[0],a[2]), List(b[0],b[2])),
      cross(List(a[0],a[1]), List(b[0],b[1]))
      );
  }
}

Expr curl(const Expr& f)
{
  Expr del = gradient(f.size());

  return cross(del, f);
}

Expr colonProduct(const Expr& A, const Expr& B)
{
  int nA = 0;
  int nB = 0;
  TEUCHOS_TEST_FOR_EXCEPTION(!isSquareMatrix(A, nA), std::runtime_error,
    "Colon product expected argument A=" << A << " to be a square matrix");
  TEUCHOS_TEST_FOR_EXCEPTION(!isSquareMatrix(B, nB), std::runtime_error,
    "Colon product expected argument B=" << B << " to be a square matrix");

  TEUCHOS_TEST_FOR_EXCEPTION(nA!=nB, std::runtime_error,
    "Colon product expected operands A=" << A << " and B=" << B 
    << " to have identical sizes");

  Expr rtn = 0.0;
  for (int i=0; i<nA; i++)
  {
    for (int j=0; j<nA; j++)
    {
      rtn = rtn + A[i][j] * B[i][j];
    }
  }

  return rtn;
}

Expr outerProduct(const Expr& A, const Expr& B)
{
  int nA = 0;
  int nB = 0;
  TEUCHOS_TEST_FOR_EXCEPTION(!isVector(A, nA), std::runtime_error,
    "Outer product expected argument A=" << A << " to be a vector");
  TEUCHOS_TEST_FOR_EXCEPTION(!isVector(B, nB), std::runtime_error,
    "Outer product expected argument B=" << B << " to be a vector");

  TEUCHOS_TEST_FOR_EXCEPTION(nA!=nB, std::runtime_error,
    "Colon product expected operands A=" << A << " and B=" << B 
    << " to have identical sizes");

  Array<Expr> rtn(nA);
  for (int i=0; i<nA; i++)
  {
    rtn[i] = A[i] * B;
  }

  return new ListExpr(rtn);
}

bool isVector(const Expr& x, int& N)
{
  N = 0;
  if (x.size() == x.totalSize()) 
  {
    N = x.size();
    return true;
  }
  else
  {
    return false;
  }
}


bool isSquareMatrix(const Expr& x, int& N)
{
  N = 0;
  /* do the simplest checks first */
  if (x.size()*x.size() != x.totalSize()) return false;
  N = x.size();
  for (int i=0; i<N; i++)
  {
    int M = 0;
    if (!isVector(x[i],M)) return false;
    if (M != N) return false;
  }
  return true;
}




}

