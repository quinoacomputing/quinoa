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

#include "SundanceHermiteSpectralBasis.hpp"
#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMap.hpp"

using namespace std;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


HermiteSpectralBasis::HermiteSpectralBasis(int dim, int order)
  : SpectralBasisBase(),
    basis_(),
    dim_(dim),
    order_(order),
    maxterms_(-1),
    cijk_()
{
  Chaos Hermite(dim_, order_);
  maxterms_ = Hermite.tnterms();
  basis_.resize(maxterms_);
  for (int i=0; i<maxterms_; i++)
    {
      basis_[i] = i;
    }  
   cijk_ = rcp(new cijk(dim_, order_));
}

HermiteSpectralBasis::HermiteSpectralBasis(int dim, int order, int nterms)
  : SpectralBasisBase(),
    basis_(),
    dim_(dim),
    order_(order),
    maxterms_(nterms),
    cijk_()
{
  basis_.resize(maxterms_);
  for (int i=0; i<maxterms_; i++)
    {
      basis_[i] = i;
    }
  cijk_ = rcp(new cijk(dim_, order_));
}

HermiteSpectralBasis::HermiteSpectralBasis(int dim, int order, const Array<int>& basisarray)
  : SpectralBasisBase(),
    basis_(),
    dim_(dim),
    order_(order),
    maxterms_(-1),
    cijk_()
{
  maxterms_ = basisarray.size();
  basis_.resize(maxterms_);
  for (int i=0; i<maxterms_; i++)
    basis_[i] = basisarray[i];
  cijk_ = rcp(new cijk(dim_, order_));
}


int HermiteSpectralBasis::getDim() const 
{
  return dim_;
}

int HermiteSpectralBasis::getOrder() const
{
  return order_;
}

int HermiteSpectralBasis::nterms() const
{
  return maxterms_;
}

int HermiteSpectralBasis::getElement(int i) const
{
  return basis_[i];
}

double HermiteSpectralBasis::expectation(int i, int j, int k)
{
   return cijk_->expectation(i,j,k);
}


string HermiteSpectralBasis::toString() const
{
  return "HermiteSpectralBasis(" + Teuchos::toString(getDim())
    + ", " + Teuchos::toString(getOrder()) + ")";
}


bool HermiteSpectralBasis::lessThan(const SpectralBasisBase* other) const
{
  if (typeid(*this).before(typeid(*other))) return true;
  if (typeid(*other).before(typeid(*this))) return false;

  if (getDim() < other->getDim()) return true;
  if (other->getDim() < getDim()) return false;

  if (getOrder() < other->getOrder()) return true;
  if (other->getOrder() < getOrder()) return false;

  return false;
}
