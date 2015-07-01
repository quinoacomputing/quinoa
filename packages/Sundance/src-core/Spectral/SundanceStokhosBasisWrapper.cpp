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

#include "SundanceStokhosBasisWrapper.hpp"

#ifdef HAVE_SUNDANCE_STOKHOS

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMap.hpp"
#include "PlayaExceptions.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"


namespace Sundance
{
using namespace Teuchos;


StokhosBasisWrapper::StokhosBasisWrapper(
  const RCP<const PolyBasis>& basis)
  : basis_(basis), cijk_()
{
  fillCijk();
}

StokhosBasisWrapper::StokhosBasisWrapper(
  const RCP<const PolyBasis1D>& basis)
  : basis_(), cijk_()
{
  Array<RCP<const PolyBasis1D> > bases(1); 
  bases[0] = basis;
  basis_ = rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
  fillCijk();
}

void StokhosBasisWrapper::fillCijk()
{
  cijk_ = basis_->computeTripleProductTensor();
}


double StokhosBasisWrapper::expectation(int i, int j, int k)
{
  return cijk_->getValue(i,j,k);
}




bool StokhosBasisWrapper::lessThan(const SpectralBasisBase* other) const
{
  if (typeid(*this).before(typeid(*other))) return true;
  if (typeid(*other).before(typeid(*this))) return false;

  if (getDim() < other->getDim()) return true;
  if (other->getDim() < getDim()) return false;

  if (getOrder() < other->getOrder()) return true;
  if (other->getOrder() < getOrder()) return false;

  const PolyBasis* me = basis_.ptr().get();
  /* because the type comparisons have been neutral, this cast should not
   * fail */
  const StokhosBasisWrapper* otherGuy 
    = dynamic_cast<const StokhosBasisWrapper*>(other);
  TEUCHOS_TEST_FOR_EXCEPTION(otherGuy == 0, std::runtime_error,
    "unexpected cast failure");

  const PolyBasis* you = otherGuy->basis_.ptr().get();
  if (typeid(*me).before(typeid(*you))) return true;

  return false;
}


}

#endif
