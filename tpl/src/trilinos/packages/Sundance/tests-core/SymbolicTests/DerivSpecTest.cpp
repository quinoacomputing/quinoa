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
#include "SundanceFunctionIdentifier.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceSymbolicFunc.hpp"
#include "SundanceUnknownParameter.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace std;

#define TEST_THROW(code, passFail) \
  TEUCHOS_TEST_THROW( code, std::exception, Out::os(), passFail)

#define TEST_NOTHROW(code, passFail) \
  TEUCHOS_TEST_NOTHROW( code, Out::os(), passFail)

Deriv funcDeriv(const Expr& u, const MultiIndex& mi = MultiIndex(0,0,0))
{
  if (u.size()==1)
  {
    const SymbolicFuncElement* sfe 
      = dynamic_cast<const SymbolicFuncElement*>(u[0].ptr().get());
    TEUCHOS_TEST_FOR_EXCEPT(sfe==0);
    return Sundance::funcDeriv(sfe, mi);
  }
  else
  {
    const SymbolicFunc* sf 
      = dynamic_cast<const SymbolicFunc*>(u.ptr().get());
    TEUCHOS_TEST_FOR_EXCEPT(sf==0);
    return Sundance::funcDeriv(sf);
  }
}

Deriv normalDeriv(const Expr& u)
{
  const SymbolicFuncElement* sfe 
    = dynamic_cast<const SymbolicFuncElement*>(u[0].ptr().get());
  TEUCHOS_TEST_FOR_EXCEPT(sfe==0);
  return normalDeriv(sfe);
}

Deriv divergenceDeriv(const Expr& u)
{
  const SymbolicFunc* sf 
    = dynamic_cast<const SymbolicFunc*>(u.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPT(sf==0);
  return divergenceDeriv(sf);
}

FunctionIdentifier fid(const Expr& u)
{
  const SymbolicFuncElement* sfe 
    = dynamic_cast<const SymbolicFuncElement*>(u[0].ptr().get());
  const SymbolicFunc* sf 
    = dynamic_cast<const SymbolicFunc*>(u.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPT(sf==0 && sfe==0);
  if (sf) return sf->fid();
  else return sfe->fid();

}


bool testVecFunction()
{
  bool pass = true;
  Tabs tab;
  int verb=1;

  /* make a vector function */
  int dim = 3;
  Expr u = new UnknownFunctionStub("u", 1, 3);
  Expr phi = new UnknownFunctionStub("u", 0, 1);
  Expr v = new TestFunctionStub("v", 1, 3);
  Expr alpha = new UnknownParameter("alpha");

  TEUCHOS_TEST_EQUALITY(u.size(), dim, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(v.size(), dim, Out::os(), pass);
  
  /* Verify that D_mi u fails for nonzero multiindex.  */
//  MultiIndex m1(0,1,0);
//  TEST_THROW(
//    funcDeriv(u, m1), pass
//    );
  
  
  Deriv d_du = funcDeriv(u);
  Deriv d_dAlpha = funcDeriv(alpha);
  Deriv d_dphi = funcDeriv(phi);
  Deriv d_dphi_x = d_dphi.derivWrtMultiIndex(MultiIndex(1,0,0));
  Deriv divU = divergenceDeriv(u);
  Deriv divV = divergenceDeriv(v);

  TEUCHOS_TEST_EQUALITY(divU.fid(), fid(u), Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dAlpha.fid(), fid(alpha), Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking algSpec for parameter");
  TEUCHOS_TEST_EQUALITY(d_dAlpha.algSpec().isScalar(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dAlpha.algSpec().isVector(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dAlpha.algSpec().isNormal(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dAlpha.algSpec().isCoordinateComponent(), false, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking algSpec for divergence");
  TEUCHOS_TEST_EQUALITY(divU.algSpec().isScalar(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.algSpec().isVector(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.algSpec().isNormal(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.algSpec().isCoordinateComponent(), false, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking algSpec for vector operative func");
  TEUCHOS_TEST_EQUALITY(d_du.algSpec().isScalar(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.algSpec().isVector(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.algSpec().isNormal(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.algSpec().isCoordinateComponent(), false, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking algSpec for scalar operative func");
  TEUCHOS_TEST_EQUALITY(d_dphi.algSpec().isScalar(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.algSpec().isVector(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.algSpec().isNormal(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.algSpec().isCoordinateComponent(), false, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking algSpec for d/d(Dx(phi))");
  TEUCHOS_TEST_EQUALITY(d_dphi_x.algSpec().isScalar(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.algSpec().isVector(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.algSpec().isNormal(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.algSpec().isCoordinateComponent(), true, Out::os(), pass);


  SUNDANCE_BANNER1(verb, tab, "checking distinction between functional and spatial derivs");
  TEUCHOS_TEST_EQUALITY(d_dAlpha.isFunctionalDeriv(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dAlpha.isCoordDeriv(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.isFunctionalDeriv(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.isCoordDeriv(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.isFunctionalDeriv(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.isCoordDeriv(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.isFunctionalDeriv(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.isCoordDeriv(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.isFunctionalDeriv(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.isCoordDeriv(), false, Out::os(), pass);


  SUNDANCE_BANNER1(verb, tab, "checking identification of differentiation order");
  TEUCHOS_TEST_EQUALITY(divU.opOnFunc().derivOrder(), 1, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.opOnFunc().derivOrder(), 0, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dAlpha.opOnFunc().derivOrder(), 0, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.opOnFunc().derivOrder(), 0, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.opOnFunc().derivOrder(), 1, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking identification of test and unk funcs");
  TEUCHOS_TEST_EQUALITY(d_du.isTestFunction(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.isUnknownFunction(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.isParameter(), false, Out::os(), pass);

  TEUCHOS_TEST_EQUALITY(d_dphi.isTestFunction(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.isUnknownFunction(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.isParameter(), false, Out::os(), pass);

  TEUCHOS_TEST_EQUALITY(d_dphi_x.isTestFunction(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.isUnknownFunction(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.isParameter(), false, Out::os(), pass);

  TEUCHOS_TEST_EQUALITY(divU.isTestFunction(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.isUnknownFunction(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.isParameter(), false, Out::os(), pass);

  TEUCHOS_TEST_EQUALITY(divV.isTestFunction(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divV.isUnknownFunction(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divV.isParameter(), false, Out::os(), pass);

  /* unknown parameters are both unknown and parameters */
  TEUCHOS_TEST_EQUALITY(d_dAlpha.isTestFunction(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dAlpha.isUnknownFunction(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dAlpha.isParameter(), true, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking identification of spatial operators acting on operative functions");

  TEUCHOS_TEST_EQUALITY(divU.opOnFunc().isDivergence(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.opOnFunc().isNormal(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.opOnFunc().isPartial(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.opOnFunc().isIdentity(), false, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking that asking for component information from a divergence throws an exception");

  TEST_THROW(divU.opOnFunc().normalDerivOrder(), pass);
  TEST_THROW(divU.opOnFunc().mi(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking that applying a partial derivative to a divergence throws an exception");
  TEST_THROW(divU.opOnFunc().derivWrtMultiIndex(MultiIndex(1,0,0)), pass);
  TEST_THROW(divU.opOnFunc().derivWrtMultiIndex(MultiIndex(0,1,0)), pass);
  TEST_THROW(divU.opOnFunc().derivWrtMultiIndex(MultiIndex(0,0,1)), pass);
  SUNDANCE_BANNER1(verb, tab, "checking that applying a zero-order partial derivative to a divergence does not throw an exception");
  TEST_NOTHROW(divU.opOnFunc().derivWrtMultiIndex(MultiIndex(0,0,0)), pass);

  SUNDANCE_BANNER1(verb, tab, "checking that applying a partial derivative to a divergence throws an exception");
  TEST_THROW(divU.derivWrtMultiIndex(MultiIndex(1,0,0)), pass);
  TEST_THROW(divU.derivWrtMultiIndex(MultiIndex(0,1,0)), pass);
  TEST_THROW(divU.derivWrtMultiIndex(MultiIndex(0,0,1)), pass);
  SUNDANCE_BANNER1(verb, tab, "checking that applying a zero-order partial derivative to a divergence does not throw an exception");
  TEST_NOTHROW(divU.derivWrtMultiIndex(MultiIndex(0,0,0)), pass);


  SUNDANCE_BANNER1(verb, tab, "checking that applying a partial derivative to a parameter throws an exception");
  TEST_THROW(d_dAlpha.derivWrtMultiIndex(MultiIndex(1,0,0)), pass);
  TEST_THROW(d_dAlpha.derivWrtMultiIndex(MultiIndex(0,1,0)), pass);
  TEST_THROW(d_dAlpha.derivWrtMultiIndex(MultiIndex(0,0,1)), pass);
  SUNDANCE_BANNER1(verb, tab, "checking that applying a zero-order partial derivative to a divergence does not throw an exception");
  TEST_NOTHROW(d_dAlpha.derivWrtMultiIndex(MultiIndex(0,0,0)), pass);

  SUNDANCE_BANNER1(verb, tab, "all done!");

  return pass;
}

int main(int argc, char** argv)
{
  bool pass = true;
  try
    {
      pass = pass && testVecFunction();
    }
  catch(std::exception& e)
    {
      pass = false;
      Out::os() << "unexpected exception: " << e.what() << std::endl;
    }
  if (pass) 
  {
    std::cerr << "test PASSED" << std::endl;
    return 0;
  }
  else 
  {
    std::cerr << "test FAILED" << std::endl;
    return -1;
  }
}

