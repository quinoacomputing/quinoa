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

#include "SundanceFunctionIdentifier.hpp"
#include "SundanceOrderedTuple.hpp"

using namespace Sundance;
using namespace Sundance;

FunctionIdentifier::FunctionIdentifier()
  : dofID_(-1), algSpec_(ScalarAT)
{}

FunctionIdentifier::FunctionIdentifier(
  const AlgebraSpecifier& algSpec)
  : dofID_(nextID()), algSpec_(algSpec)
{}

FunctionIdentifier::FunctionIdentifier(
  const FunctionIdentifier* parent,
  const AlgebraSpecifier& algSpec)
  : dofID_(-1), algSpec_(algSpec)
{
  /* make sure the parent exists */
  TEUCHOS_TEST_FOR_EXCEPT(parent==0);
  dofID_ = parent->dofID();

  /* check for various stupid cases that should never happen */
  TEUCHOS_TEST_FOR_EXCEPTION(!parent->algSpec_.isVector(), std::logic_error,
    "attempted to form a function ID for a component of a non-vector object:"
    "parent=" << parent->toString() << " component spec=" << algSpec);
  TEUCHOS_TEST_FOR_EXCEPTION(algSpec.isVector() || algSpec.isScalar(), std::runtime_error,
    "attempted to define a vector or scalar as a component of another object."
    "parent=" <<  parent->toString() << " component spec=" << algSpec);
}



int FunctionIdentifier::componentIndex() const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!(algSpec_.isCoordinateComponent() || algSpec_.isScalar()), std::logic_error,
    "attempted to find component index for a FID that is not a "
    "scalar or a coordinate component of a vector");
  if (algSpec_.isScalar()) return 0;
  return algSpec_.direction();
}

string FunctionIdentifier::toString() const 
{
  TeuchosOStringStream os;
  os << *this;
  return os.str();
}

FunctionIdentifier FunctionIdentifier::createNormal() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isVector(), std::logic_error,
    "attempted to find normal component of a FID that is not a vector");

  return FunctionIdentifier(this, normalAlgebraSpec());
}

bool FunctionIdentifier::operator<(const FunctionIdentifier& other) const 
{
  OrderedPair<int, AlgebraSpecifier> me(dofID_, algSpec_);
  OrderedPair<int, AlgebraSpecifier> you(other.dofID_, other.algSpec_);
  return me < you;
}

FunctionIdentifier FunctionIdentifier::createComponent(int d) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isVector(), std::logic_error,
    "attempted to find component of a FID that is not a vector");

  return FunctionIdentifier(this, coordAlgebraSpec(d));
}

namespace Sundance
{

FunctionIdentifier makeFuncID(int tensorOrder)
{
  if (tensorOrder==0)
    return FunctionIdentifier(scalarAlgebraSpec());
  else if (tensorOrder==1)
    return FunctionIdentifier(vectorAlgebraSpec());
  else
    TEUCHOS_TEST_FOR_EXCEPT(true);
  return FunctionIdentifier(scalarAlgebraSpec()); // -Wall
}

}


namespace std
{

ostream& operator<<(std::ostream& os, 
  const Sundance::FunctionIdentifier& fid)
{
  os << "FuncID(dofID=" << fid.dofID() << ", component type=" << fid.algSpec() << ")";
  return os;
}

}

