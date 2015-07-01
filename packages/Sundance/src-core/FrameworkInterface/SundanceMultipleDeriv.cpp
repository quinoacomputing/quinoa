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

#include "SundanceMultipleDeriv.hpp"
#include "SundanceSymbolicFuncElement.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

MultipleDeriv::MultipleDeriv()
  : MultiSet<Deriv>()
{}

MultipleDeriv::MultipleDeriv(const Deriv& d)
  : MultiSet<Deriv>()
{
  put(d);
}
MultipleDeriv::MultipleDeriv(const Deriv& d1, const Deriv& d2)
  : MultiSet<Deriv>()
{
  put(d1);
  put(d2);
}

int MultipleDeriv::spatialOrder() const 
{
  int rtn = 0;
  for (MultipleDeriv::const_iterator i=this->begin(); i!=this->end(); i++)
  {
    if (i->isCoordDeriv())
    {
      rtn += 1;
    }
  }
  return rtn;
}

MultiIndex MultipleDeriv::spatialDeriv() const
{
  MultiIndex rtn;
  for (MultipleDeriv::const_iterator i=this->begin(); i!=this->end(); i++)
  {
    if (i->isCoordDeriv())
    {
      int d = i->coordDerivDir();
      rtn[d] += 1;
    }
  }
  return rtn;
}

MultiSet<FunctionIdentifier> MultipleDeriv::funcIDs() const
{
  MultiSet<FunctionIdentifier> rtn;
  for (MultipleDeriv::const_iterator i=this->begin(); i!=this->end(); i++)
  {
    if (i->isFunctionalDeriv())
    {
      rtn.put(i->fid());
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!i->isFunctionalDeriv(), std::runtime_error,
      "MultipleDeriv::funcIDs() found spatial deriv");
  }
  return rtn;
}

MultiSet<int> MultipleDeriv::dofIDs() const
{
  MultiSet<int> rtn;
  for (MultipleDeriv::const_iterator i=this->begin(); i!=this->end(); i++)
  {
    if (i->isFunctionalDeriv())
    {
      int f = i->dofID();
      rtn.put(f);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!i->isFunctionalDeriv(), std::runtime_error,
      "MultipleDeriv::sharedFuncIDs() found spatial deriv");
  }
  return rtn;
}

MultipleDeriv MultipleDeriv::product(const MultipleDeriv& other) const 
{
  MultipleDeriv rtn;
  
  for (MultipleDeriv::const_iterator i=this->begin(); i!=this->end(); i++)
  {
    rtn.put(*i);
  }
  for (MultipleDeriv::const_iterator i=other.begin(); i!=other.end(); i++)
  {
    rtn.put(*i);
  }
  return rtn;
}

MultipleDeriv MultipleDeriv::factorOutDeriv(const Deriv& x) const
{
  MultipleDeriv rtn = *this;

  MultipleDeriv::iterator i = rtn.find(x);

  /* remove a single copy of the given derivative */
  if (i != rtn.end()) rtn.erase(i);

  if (rtn.size() == this->size()) return MultipleDeriv();

  return rtn;
}


bool MultipleDeriv::containsDeriv(const MultipleDeriv& x) const
{
  for (MultipleDeriv::const_iterator i=x.begin(); i!=x.end(); i++)
  {
    if ( count(*i) <= x.count(*i) ) return false;
  }
  return true;
}

MultipleDeriv MultipleDeriv::factorOutDeriv(const MultipleDeriv& x) const
{
  MultipleDeriv rtn = *this;

  for (MultipleDeriv::const_iterator i=x.begin(); i!=x.end(); i++)
  {
    MultipleDeriv::iterator j = rtn.find(*i);

    /* remove a single copy of the given derivative */
    if (j != rtn.end()) rtn.erase(j);
  }

  if (rtn.size() == this->size()) return MultipleDeriv();
  return rtn;
}

bool MultipleDeriv
::isInRequiredSet(const Set<MultiSet<int> >& funcCombinations,
  const Set<MultiIndex>& multiIndices) const
{
  if (spatialOrder() == 0)
  {
    return funcCombinations.contains(dofIDs());
  }
  else
  {
    return multiIndices.contains(spatialDeriv());
  }
}


void MultipleDeriv
::productRulePermutations(ProductRulePerms& perms) const 
{
  int N = order();

  if (N==0)
  {
    MultipleDeriv md0;
    DerivPair p(md0, md0);
    perms.put(p, 1);
    return;
  }

  int p2 = pow2(N);

  for (int i=0; i<p2; i++)
  {
    MultipleDeriv left;
    MultipleDeriv right;
    Array<int> bits = bitsOfAnInteger(i, N);
    int j=0; 
    MultipleDeriv::const_iterator iter;
    for (iter=this->begin(); iter != this->end(); iter++, j++)
    {
      if (bits[j]==true)
      {
        left.put(*iter);
      }
      else
      {
        right.put(*iter);
      }
    }
    DerivPair p(left, right);
    if (!perms.containsKey(p))
    {
      perms.put(p, 1);
    }
    else
    {
      int count = perms.get(p);
      perms.put(p, count+1);
    }
  }
}

Array<int> MultipleDeriv::bitsOfAnInteger(int x, int n)
{
  TEUCHOS_TEST_FOR_EXCEPTION(x >= pow2(n), std::logic_error,
    "Invalid input to MultipleDeriv::bitsOfX");
                     
  Array<int> rtn(n);

  int r = x;
  for (int b=n-1; b>=0; b--)
  {
    rtn[b] = r/pow2(b);
    r = r - rtn[b]*pow2(b);
  }
  return rtn;
}

int MultipleDeriv::pow2(int n)
{
  static Array<int> p2(1,1);

  if (n >= p2.size())
  {
    int oldN = p2.size(); 
    for (int i=oldN; i<=n; i++) p2.push_back(p2[i-1]*2);
  }
  
  return p2[n];
}



namespace Sundance
{
Set<MultipleDeriv> applyTx(const Set<MultipleDeriv>& s,
  const MultiIndex& x)
{
  Set<MultipleDeriv> rtn;

  for (Set<MultipleDeriv>::const_iterator i=s.begin(); i!=s.end(); i++)
  {
    const MultipleDeriv& md = *i;
    for (MultipleDeriv::const_iterator j=md.begin(); j!=md.end(); j++)
    {
      const Deriv& d = *j;
      if (d.isFunctionalDeriv())
      {
        const MultiIndex& mi = d.opOnFunc().mi();
        MultiIndex miNew = mi+x;
        if (miNew.isValid())
        {
          Deriv dNew = d.derivWrtMultiIndex(miNew);
          MultipleDeriv mdNew = md;
          mdNew.erase(mdNew.find(d));
          mdNew.put(dNew);
          rtn.put(mdNew);
        }
      }
    }
  }
  return rtn;
}

Set<MultipleDeriv> Xx(const MultiIndex& x)
{
  Set<MultipleDeriv> rtn;

  TEUCHOS_TEST_FOR_EXCEPTION(x.order() < 0 || x.order() > 1, std::logic_error,
    "invalid multiindex " << x << " in this context");

  MultipleDeriv xmd = makeMultiDeriv(coordDeriv(x.firstOrderDirection()));
  rtn.put(xmd);
  return rtn;
}

Set<MultipleDeriv> applyZx(const Set<MultipleDeriv>& W,
  const MultiIndex& x)
{
  Set<MultipleDeriv> rtn;

  TEUCHOS_TEST_FOR_EXCEPTION(x.order() < 0 || x.order() > 1, std::logic_error,
    "invalid multiindex " << x << " in this context");

  for (Set<MultipleDeriv>::const_iterator i=W.begin(); i!=W.end(); i++)
  {
    const MultipleDeriv& md = *i;
    TEUCHOS_TEST_FOR_EXCEPTION(md.order() != 1, std::logic_error,
      "Only first-order multiple functional derivatives "
      "should appear in this function. The derivative "
      << md << " is not first-order.");

    const Deriv& d = *(md.begin());

    if (d.isFunctionalDeriv())
    {
      /* */
      TEUCHOS_TEST_FOR_EXCEPTION(!d.canBeDifferentiated(),
        std::logic_error, "function signature " << d << " cannot be "
        "differentiated further spatially");
      /* accept a functional derivative if the associated function 
       * is not identically zero */
      const SymbolicFuncElement* sfe = d.symbFuncElem();
      TEUCHOS_TEST_FOR_EXCEPTION(sfe==0, std::logic_error, 
        "can't cast function in "
        << d << " to a SymbolicFuncElement");
      if (sfe && !sfe->evalPtIsZero()) rtn.put(md);
    }
  }
  return rtn;
}


    
int factorial(const MultipleDeriv& ms)
{
  Sundance::Map<Deriv, int> counts;
    
  for (MultipleDeriv::const_iterator i=ms.begin(); i!=ms.end(); i++)
  {
    if (counts.containsKey(*i)) counts[*i]++;
    else counts.put(*i, 1);
  }

  int rtn = 1;
  for (Sundance::Map<Deriv, int>::const_iterator
         i=counts.begin(); i!=counts.end(); i++)
  {
    int f = 1;
    for (int j=1; j<=i->second; j++) f *= j;
    rtn *= f;
  }
  return rtn;
}

MultipleDeriv makeMultiDeriv(const Deriv& d)
{
  MultipleDeriv rtn;
  rtn.put(d);
  return rtn;
}

bool hasParameter(const MultipleDeriv& d)
{
  for (MultipleDeriv::const_iterator i=d.begin(); i!=d.end(); i++)
  {
    if (i->isParameter()) return true;
  }
  return false;
}



}
