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

#include "SundanceDefs.hpp"
#include "SundanceSpectralExpr.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceComplexExpr.hpp"
#include "PlayaExceptions.hpp"
#include "Teuchos_Array.hpp"

using namespace Teuchos;
using namespace Sundance;

SpectralExpr::SpectralExpr(const SpectralBasis& sbasis, const Array<Expr>& coeffs)
  : ScalarExpr(), 
    coeffs_(coeffs),
    sbasis_(sbasis)

{
  TEUCHOS_TEST_FOR_EXCEPT(coeffs_.size() != sbasis_.nterms());
}


SpectralExpr::SpectralExpr(const SpectralBasis& sbasis, const Expr& coeffs)
  : ScalarExpr(), 
    coeffs_(),
    sbasis_(sbasis)

{
  int nterms = sbasis.nterms();
  coeffs_.resize(nterms);
  for (int i=0; i<nterms; i++)
    coeffs_[i] = coeffs[i];
}

void SpectralExpr::accumulateFuncSet(Set<int>& funcDofIDs, 
  const Set<int>& activeSet) const
{
  for (int i=0; i<coeffs_.size(); i++)
  {
    dynamic_cast<const ScalarExpr*>(coeffs_[i].ptr().get())->accumulateFuncSet(funcDofIDs, activeSet);
  }
}

SpectralBasis SpectralExpr::getSpectralBasis() const
{ 
  return sbasis_;
}


Expr SpectralExpr::getCoeff(int i) const
{
  return coeffs_[i];
}

Expr SpectralExpr::spectralDotProduct(const SpectralExpr* other) const
{
  Expr rtn = coeffs_[0] * other->coeffs_[0];
  for (int i=1; i<coeffs_.size(); i++) rtn = rtn + coeffs_[i]*other->coeffs_[i];
  return rtn;
}

bool SpectralExpr::hasTestFunctions() const
{
  bool rtn = coeffs_[0].ptr()->hasTestFunctions();
  for (int i=1; i<coeffs_.size(); i++)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(coeffs_[i].ptr()->hasTestFunctions() != rtn,
                         std::logic_error,
                         "expr " << toString() << " has a mix of test and "
                         "non-test coefficients");
    }
  return rtn;
}

bool SpectralExpr::hasUnkFunctions() const
{
  bool rtn = coeffs_[0].ptr()->hasUnkFunctions();
  for (int i=1; i<coeffs_.size(); i++)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(coeffs_[i].ptr()->hasUnkFunctions() != rtn,
                         std::logic_error,
                         "expr " << toString() << " has a mix of unk and "
                         "non-unk coefficients");
    }
  return rtn;
}

bool SpectralExpr::hasHungryDiffOp() const
{
  for (int i=0; i<coeffs_.size(); i++)
    {
      Expr re = coeffs_[i].real();
      Expr im = coeffs_[i].imag();
      const ScalarExpr* sr = dynamic_cast<const ScalarExpr*>(re.ptr().get());
      const ScalarExpr* si = dynamic_cast<const ScalarExpr*>(im.ptr().get());
      TEUCHOS_TEST_FOR_EXCEPTION(sr == 0 || si == 0, std::logic_error,
                         "spectral expr " << toString() << " contains a "
                         "non-scalar coefficient");
      if (sr->isHungryDiffOp() || si->isHungryDiffOp()) return true;
    }
  return false;
}


std::ostream& SpectralExpr::toText(std::ostream& os, bool paren) const
{
  os << "SpectralExpr{";
  for (int i=0; i<coeffs_.size(); i++)
    {
      coeffs_[i].ptr()->toText(os, paren);
      if (i < coeffs_.size()-1) os << ", ";
    }
  os << "}";
  return os;
}


XMLObject SpectralExpr::toXML() const 
{
  XMLObject rtn("SpectralExpr");
  for (int i=0; i<coeffs_.length(); i++)
    {
      rtn.addChild(coeffs_[i].toXML());
    }
  return rtn;
}


bool  SpectralExpr::lessThan(const ScalarExpr* other) const
{
  const SpectralExpr* s = dynamic_cast<const SpectralExpr*>(other);
  TEUCHOS_TEST_FOR_EXCEPTION(s==0, std::logic_error, "cast should never fail at this point");
  if (coeffs_.size() < s->coeffs_.size()) return true;
  if (coeffs_.size() > s->coeffs_.size()) return false;
  for (int i=0; i<coeffs_.size(); i++)
  {
    if (coeffs_[i].lessThan(s->coeffs_[i])) return true;
    if (s->coeffs_[i].lessThan(coeffs_[i])) return false;
  }
  return sbasis_.ptr()->lessThan(s->sbasis_.ptr().get());
}


namespace Sundance
{
  /** */
  Expr getSpectralCoeff(int i, const Expr& e)
  {
    const SpectralExpr* s 
      = dynamic_cast<const SpectralExpr*>(e[0].ptr().get());
    TEUCHOS_TEST_FOR_EXCEPTION(s!=0, std::runtime_error,
                       "getSpectralCoeff() called on non-spectral expr "
                       << e.toString());
    return s->getCoeff(i);
  }

  /** */
  SpectralBasis getSpectralBasis(const Expr& e) 
  {
    const SpectralExpr* s 
      = dynamic_cast<const SpectralExpr*>(e[0].ptr().get());
    TEUCHOS_TEST_FOR_EXCEPTION(s!=0, std::runtime_error,
                       "getSpectralBasis() called on non-spectral expr "
                       << e.toString());
    return s->getSpectralBasis();
  }
}
