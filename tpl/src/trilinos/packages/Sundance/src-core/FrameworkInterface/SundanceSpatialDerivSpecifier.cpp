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

#include "SundanceSpatialDerivSpecifier.hpp"

using namespace Sundance;
using namespace Sundance;



SpatialDerivSpecifier::SpatialDerivSpecifier()
  : EnumTypeField<SpatialDerivType>(IdentitySDT), 
    mi_(), normalDerivOrder_(-1)
{}

SpatialDerivSpecifier::SpatialDerivSpecifier(const MultiIndex& mi)
  : EnumTypeField<SpatialDerivType>(PartialSDT), 
    mi_(mi), normalDerivOrder_(-1)
{}

SpatialDerivSpecifier::SpatialDerivSpecifier(
  const SpatialDerivType& sdt,
  int order)
  : EnumTypeField<SpatialDerivType>(sdt),
    mi_(), normalDerivOrder_(order)
{
  assertNotType(PartialSDT);

  if (order > 0)
  {
    assertNotType(DivSDT);
  }
}

const MultiIndex& SpatialDerivSpecifier::mi() const
{
  assertNotType(DivSDT);
  assertNotType(NormalSDT);
  return mi_;
}

bool SpatialDerivSpecifier::isDivergence() const
{
  return isType(DivSDT);
}

bool SpatialDerivSpecifier::isNormal() const
{
  return isType(NormalSDT);
}

bool SpatialDerivSpecifier::isPartial() const
{
  return isType(PartialSDT);
}

bool SpatialDerivSpecifier::isIdentity() const
{
  return isType(IdentitySDT)
    || (isPartial() && mi().order()==0)
    || (isNormal() && normalDerivOrder()==0);
}

int SpatialDerivSpecifier::normalDerivOrder() const
{
  assertType(NormalSDT);
  return normalDerivOrder_;
}

int SpatialDerivSpecifier::derivOrder() const
{
  if (isDivergence()) return 1;
  if (isPartial()) return mi_.order();
  if (isNormal()) return normalDerivOrder_;
  return 0;
}


std::string SpatialDerivSpecifier::toString() const 
{
  TeuchosOStringStream os;
  os << *this;
  return os.str();
}



bool SpatialDerivSpecifier::operator<(const SpatialDerivSpecifier& other) const
{
  if (type() < other.type()) return true;
  if (type() > other.type()) return false;

  if (isPartial()) return mi() < other.mi();
  if (isNormal()) return normalDerivOrder() < other.normalDerivOrder();

  return false;
}


SpatialDerivSpecifier SpatialDerivSpecifier::derivWrtMultiIndex(const MultiIndex& mi) const 
{
  if (isPartial() || isIdentity()) 
  {
    return SpatialDerivSpecifier(mi_+mi);
  }
  else if (mi.order()==0)
  {
    return *this;
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "cannot take an arbitrary "
      "spatial derivative of SDS=" << *this);
    return *this; // -Wall
  }
  
}


namespace std
{
/** \relates SpatialDerivSpecifier */
ostream& operator<<(std::ostream& os, const Sundance::SpatialDerivSpecifier& sds)
{
  os << sds.type();
  if (sds.isPartial()) os << "(d=" << sds.mi() << ")";
  else os << "()";
  return os;
}


/** \relates SpatialDerivSpecifier */
ostream& operator<<(std::ostream& os, const Sundance::SpatialDerivType& sdt)
{
  static Array<string> names = tuple<string>("Identity", "Partial", "Normal", "Divergence");
  os << names[sdt];
  return os;
}

}
