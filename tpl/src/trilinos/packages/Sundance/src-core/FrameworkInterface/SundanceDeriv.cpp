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

#include "SundanceDeriv.hpp"


#include "PlayaExceptions.hpp"
#include "SundanceParameter.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceSymbolicFunc.hpp"
#include "SundanceCommonFuncDataStub.hpp"
#include "SundanceTestFuncElement.hpp"
 

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;

Deriv::Deriv()
  : EnumTypeField<DerivType>(NullDT), myAlgSpec_(ScalarAT), fid_(), sds_(), 
    symbFuncElem_(0), symbFunc_(0), coordDerivDir_(-1)
{}

Deriv::Deriv(int coordDerivDir)
  : EnumTypeField<DerivType>(CoordDT), myAlgSpec_(coordDerivDir),
    fid_(FunctionIdentifier()), sds_(), 
    symbFuncElem_(0), symbFunc_(0),
    coordDerivDir_(coordDerivDir)
{}


Deriv::Deriv(
  const SymbolicFuncElement* func,
  const SpatialDerivSpecifier& d
  )
  : EnumTypeField<DerivType>(FunctionalDT), myAlgSpec_(),
    fid_(), sds_(d),
    symbFuncElem_(func),
    symbFunc_(0), coordDerivDir_(-1)
{
  TEUCHOS_TEST_FOR_EXCEPTION(func==0, std::logic_error, 
    "null function given to Deriv ctor");
  fid_ = func->fid();
  myAlgSpec_ = derivAlgSpec(fid_.algSpec(), d);

  checkConsistencyOfOperations();
}

Deriv::Deriv( 
  const SymbolicFunc* func,
  const SpatialDerivSpecifier& d
  )
  : EnumTypeField<DerivType>(FunctionalDT), myAlgSpec_(),
    fid_(), sds_(d),
    symbFuncElem_(0), symbFunc_(func), coordDerivDir_(-1)
{
  TEUCHOS_TEST_FOR_EXCEPTION(func==0, std::logic_error, 
    "null function given to Deriv ctor");
  fid_ = func->fid();
  myAlgSpec_ = derivAlgSpec(fid_.algSpec(), d);

  checkConsistencyOfOperations();
}


void Deriv::checkConsistencyOfOperations() const 
{
  TEUCHOS_TEST_FOR_EXCEPT(sds_.isDivergence() && !fid_.isVector());
  TEUCHOS_TEST_FOR_EXCEPT(sds_.isPartial() && fid_.isVector());
  TEUCHOS_TEST_FOR_EXCEPT(sds_.isNormal() && fid_.isVector());
}


bool Deriv::operator<(const Deriv& other) const 
{
  if (type() < other.type()) return true;
  if (type() > other.type()) return false;

  if (type() == CoordDT)
  {
    return coordDerivDir() < other.coordDerivDir();
  }
  else 
  {
    OrderedPair<FunctionIdentifier, SpatialDerivSpecifier> me(fid(),sds_);
    OrderedPair<FunctionIdentifier, SpatialDerivSpecifier> you(other.fid(),other.sds_);
    return me < you;
  }
}


bool Deriv::operator==(const Deriv& other) const
{
  return !(*this < other || other < *this);
}

std::string Deriv::toString(bool verb) const 
{
  static Array<string> coords = tuple<string>("x", "y", "z");
  switch(type())
  {
    case CoordDT:
      return coords[coordDerivDir()];
    case FunctionalDT:
    {
      string f;
      if (symbFunc_!=0) f = symbFunc_->toString();
      else if (symbFuncElem_!=0) f = symbFuncElem_->toString();
      else TEUCHOS_TEST_FOR_EXCEPT(true);

      if (opOnFunc().isPartial() && opOnFunc().derivOrder()>0)
      {
        return "D[" + f + ", " 
          + coords[opOnFunc().mi().firstOrderDirection()] + "]";
      }
      else if (opOnFunc().isDivergence())
      {
        return "div(" + f + ")";
      }
      else if (opOnFunc().isNormal() && opOnFunc().derivOrder()>0)
      {
        return "D[" + f + ", normal]";
      }
      else if (opOnFunc().isIdentity()) 
      {
        return f;
      }
      else TEUCHOS_TEST_FOR_EXCEPT(1);
    }
    case NullDT:
      return "NullDeriv()";
    default:
      TEUCHOS_TEST_FOR_EXCEPT(1);
      return "NullDeriv";
  }
}

const SymbolicFuncDescriptor* Deriv::sfdPtr() const 
{
  assertType(FunctionalDT);
  /* at this point, exactly one of the two function pointers should be
   * nonzero */
  TEUCHOS_TEST_FOR_EXCEPT(symbFuncElem_ == 0 && symbFunc_==0);
  TEUCHOS_TEST_FOR_EXCEPT(symbFuncElem_ != 0 && symbFunc_!=0);
  
  /* return the nonzero pointer */
  if (symbFuncElem_) return symbFuncElem_;
  if (symbFunc_) return symbFunc_;
  return symbFunc_; // -Wall
}

bool Deriv::isTestFunction() const 
{
  if( isFunctionalDeriv() )
  {
    return sfdPtr()->isTestFunction();
  }
  return false;
}

bool Deriv::isUnknownFunction() const 
{
  if( isFunctionalDeriv() )
  {
    return sfdPtr()->isUnknownFunction();
  }
  return false;
}


bool Deriv::isParameter() const 
{
  if( isFunctionalDeriv() )
  {
    return sfdPtr()->isParameter();
  }
  return false;
}

int Deriv::dofID() const 
{
  assertType(FunctionalDT);
  return fid_.dofID();
}

const AlgebraSpecifier& Deriv::funcAlgSpec() const 
{
  assertType(FunctionalDT);
  return fid_.algSpec();
}

bool Deriv::canBeDifferentiated() const 
{
  assertType(FunctionalDT);
  if (isParameter()) return false;
  return opOnFunc().isIdentity() || opOnFunc().isPartial();
}

int Deriv::coordDerivDir() const
{
  assertType(CoordDT);
  return coordDerivDir_;
}

RCP<const CommonFuncDataStub> Deriv::data() const
{
  assertType(FunctionalDT);
  TEUCHOS_TEST_FOR_EXCEPTION(symbFuncElem_==0 && symbFunc_==0, 
    std::logic_error,
    "Deriv::data() called, but deriv=" << *this << " does not contain a "
    "valid function");
  if (symbFuncElem_) return symbFuncElem_->commonData();
  if (symbFunc_) return symbFunc_->commonData();
  return symbFunc_->commonData(); // -Wall
}

const FunctionIdentifier& Deriv::fid() const 
{
  assertType(FunctionalDT);
  return fid_;
}

const SpatialDerivSpecifier& Deriv::opOnFunc() const 
{
  assertType(FunctionalDT);
  return sds_;
}

Deriv Deriv::derivWrtMultiIndex(const MultiIndex& mi) const
{
  assertType(FunctionalDT);
  TEUCHOS_TEST_FOR_EXCEPTION(mi.order()>0 && sds_.isDivergence(), std::logic_error,
    "cannot take spatial derivative of an atomic divergence operation");
  TEUCHOS_TEST_FOR_EXCEPTION(mi.order()>0 && isParameter(), std::logic_error,
    "cannot take spatial derivative of a parameter");

  SpatialDerivSpecifier d = sds_.derivWrtMultiIndex(mi);

  if (symbFuncElem_ != 0)
  {
    return Deriv(symbFuncElem_, d);
  }
  else if (symbFunc_ != 0)
  {
    return Deriv(symbFunc_, d);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(symbFuncElem_==0 && symbFunc_==0, 
    std::logic_error,
    "attempt to differentiate a null operative function");
  return *this; // -Wall
}


AlgebraSpecifier Deriv::derivAlgSpec(
  const AlgebraSpecifier& funcAlgSpec,
  const SpatialDerivSpecifier& d)
{
  if (d.derivOrder()==0) 
  {
    return funcAlgSpec;
  }
  else if (d.isPartial()) 
  {
    return AlgebraSpecifier(d.mi().firstOrderDirection());
  }
  else if (d.isDivergence())
  {
    return ScalarAT;
  }
  else if (d.isNormal())
  {
    return NormalAT;
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPT(true);
    return IdentitySDT;
  }
}




namespace Sundance
{

Deriv coordDeriv(int d) 
{
  return Deriv(d);
}

Deriv coordDeriv(const MultiIndex& mi) 
{
  return Deriv(mi.firstOrderDirection());
}


Deriv funcDeriv(const SymbolicFuncElement* symbFunc)
{
  return Deriv(symbFunc, SpatialDerivSpecifier());
}

Deriv funcDeriv(const SymbolicFuncElement* symbFunc,
  const MultiIndex& mi)
{
  return Deriv(symbFunc, SpatialDerivSpecifier(mi));
}

Deriv funcDeriv(const SymbolicFunc* symbFunc)
{
  return Deriv(symbFunc, SpatialDerivSpecifier());
}


Deriv funcDeriv(const SymbolicFunc* symbFunc,
  const MultiIndex& mi)
{
  return Deriv(symbFunc, SpatialDerivSpecifier(mi));
}


Deriv normalDeriv(const SymbolicFuncElement* symbFunc)
{
  return Deriv(symbFunc, SpatialDerivSpecifier(NormalSDT, 1));
}

Deriv divergenceDeriv(const SymbolicFunc* symbFunc)
{
  return Deriv(symbFunc, SpatialDerivSpecifier(DivSDT, 0));
}

}

