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

#include "SundanceFunctionSupportResolver.hpp"
#include "PlayaTabs.hpp"
#include "SundanceSumOfBCs.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceOut.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceSpectralExpr.hpp"
#include "SundanceUnknownParameterElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "PlayaExceptions.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace std;

FunctionSupportResolver::FunctionSupportResolver(
  const Expr& eqns,
  const Expr& bcs,
  const Array<Expr>& vars,
  const Array<Expr>& unks,
  const Expr& unkParams,
  const Expr& fixedParams,
  const Array<Expr>& fixedFields,
  bool isVariational)
  : 
    eqns_(eqns),
    bcs_(bcs),
    integralSum_(dynamic_cast<const SumOfIntegrals*>(eqns.ptr().get())),
    bcSum_(0),
    varFuncSet_(),
    unkFuncSet_(),
    unkParamSet_(),
    fixedParamSet_(),
    regions_(),
    regionToIndexMap_(),
    varsOnRegions_(),
    unksOnRegions_(),
    bcVarsOnRegions_(),
    bcUnksOnRegions_(),
    testToRegionsMap_(),
    unkToRegionsMap_(),
    varFuncData_(),
    unkFuncData_(),
    varFuncs_(flattenSpectral(vars)),
    unkFuncs_(flattenSpectral(unks)),
    fixedFields_(flattenSpectral(fixedFields)),
    unkParams_(unkParams),
    fixedParams_(fixedParams),
    varIDToReducedIDMap_(),
    unkIDToReducedIDMap_(),
    unkParamIDToReducedUnkParamIDMap_(),
    fixedParamIDToReducedFixedParamIDMap_(),
    varIDToBlockMap_(),
    unkIDToBlockMap_(),
    unreducedVarID_(),
    unreducedUnkID_(),
    unreducedUnkParamID_(),
    unreducedFixedParamID_(),
    isVariationalProblem_(isVariational)
{
  Tabs tab0(0);
  Tabs tab1;

  /* begin with a sanity check to ensure that the input equation set 
   * exists and is integral form */

  TEUCHOS_TEST_FOR_EXCEPTION(eqns.ptr().get()==0, std::runtime_error,
    "FunctionSupportResolver ctor detected empty equation set input");

  TEUCHOS_TEST_FOR_EXCEPTION(integralSum_==0, std::runtime_error,
    "FunctionSupportResolver ctor detected an input equation set that "
    "is not in integral form");

  bool hasBCs = false;
  if (bcs.ptr().get() != 0)
  {
    bcSum_ = dynamic_cast<const SumOfBCs*>(bcs.ptr().get());
    TEUCHOS_TEST_FOR_EXCEPTION(bcSum_==0, std::runtime_error,
      "FunctionSupport ctor detected an input Essential "
      "BC that is not an EssentialBC object. "
      "Object is " << bcs);
    hasBCs = true;
  }

  /* upgrade base verbosity level if one of the terms is being watched */
  int verb = 0;
  if (integralSum_->hasWatchedTerm() || (hasBCs && bcSum_->hasWatchedTerm()))
  {
    int v1 = integralSum_->eqnSetSetupVerb();
    int v2 = 0;
    if (hasBCs) v2 = bcSum_->eqnSetSetupVerb();
    verb = std::max(verb, max(v1, v2));
  }
  SUNDANCE_BANNER1(verb, tab0, "FunctionSupportResolver setup");

  if (hasBCs)
  {
    SUNDANCE_MSG1(verb, tab1 << "Problem has EssentialBCs");
  }
  else
  {
    SUNDANCE_MSG1(verb, tab1 << "Problem has no EssentialBCs");
  }

  SUNDANCE_MSG1(verb, tab1 << "verbosity is " << verb);

  /* 
   * See whether the variational functions are TestFunction objects
   * (as in a problem where we've already taken variations, or in 
   * a Galerkin-like formulation of a non-variational problem) 
   * or UnknownFunction objects, as in a variational problem. 
   */
  bool varsAreTestFunctions = false;
  for (int b=0; b<vars.size(); b++)
  {
    for (int i=0; i<vars[b].size(); i++)
    {
      const TestFuncElement* t 
        = dynamic_cast<const TestFuncElement*>(vars[b][i].ptr().get());
      TEUCHOS_TEST_FOR_EXCEPTION(t == 0 && varsAreTestFunctions==true,
        std::runtime_error,
        "variational function list " << vars
        << " contains a mix of test and "
        "non-test functions");
      TEUCHOS_TEST_FOR_EXCEPTION(t != 0 && i>0 && varsAreTestFunctions==false,
        std::runtime_error,
        "variational function list " << vars
        << " contains a mix of test and "
        "non-test functions");
      if (t!=0) varsAreTestFunctions=true;
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION(varsAreTestFunctions == true
    && isVariationalProblem_,
    std::runtime_error,
    "test functions given to a variational equation set");

  TEUCHOS_TEST_FOR_EXCEPTION(varsAreTestFunctions == false
    && !isVariationalProblem_,
    std::runtime_error,
    "variational functions are unknowns, but equation "
    "set is not variational");

  if (isVariationalProblem_)
  {
    SUNDANCE_MSG1(verb, tab1 << "is variational problem");
  }
  else
  {
    SUNDANCE_MSG1(verb, tab1 << "is not in variational form");
  }


  /*
   * Now we collect arrays of the function data (e.g., basis functions).
   * Several function components may share common data. For example, the
   * components of a function discretized with Nedelec elements will
   * all point back to the same data. The symbolic engine needs to know 
   * about components, but the DOF map system only needs bases. 
   */
  varFuncData_.resize(vars.size());
  unkFuncData_.resize(unks.size());
  varIDToReducedIDMap_.resize(vars.size());
  unreducedVarID_.resize(vars.size());
  

  /* map each var and unknown function's ID numbers to its
   * position in the input function lists */
  for (int b=0; b<vars.size(); b++)
  {
    Tabs tab2;
    unreducedVarID_[b].resize(vars[b].size());
    int k=0;
    for (int i=0; i<vars[b].size(); i++)
    {
      const FuncElementBase* t 
        = dynamic_cast<const FuncElementBase*>(vars[b][i].ptr().get());
      int fid = t->fid().dofID();
      if (varFuncSet_.contains(fid)) continue;
      varFuncData_[b].append(getSharedFunctionData(t));
      varFuncSet_.put(fid);
      varIDToBlockMap_.put(fid, b);
      varIDToReducedIDMap_[b].put(fid, k);
      unreducedVarID_[b][k] = fid;
      k++;
    }
    SUNDANCE_MSG2(verb, tab2 << "block=" << b << " var functions are " 
      << unreducedVarID_[b]);
  }

  /* set up func ID maps for unks */
  unkIDToReducedIDMap_.resize(unks.size());
  unreducedUnkID_.resize(unks.size());
  for (int b=0; b<unks.size(); b++)
  {
    Tabs tab2;
    unreducedUnkID_[b].resize(unks[b].size());
    int k=0;
    for (int i=0; i<unks[b].size(); i++)
    {
      const UnknownFuncElement* u 
        = dynamic_cast<const UnknownFuncElement*>(unks[b][i].ptr().get());
      TEUCHOS_TEST_FOR_EXCEPTION(u==0, std::runtime_error, 
        "EquationSet ctor input unk function "
        << unks[b][i] 
        << " does not appear to be a unk function");
      int fid = u->fid().dofID();
      if (unkFuncSet_.contains(fid)) continue;
      unkFuncData_[b].append(getSharedFunctionData(u));
      unkFuncSet_.put(fid);
      unkIDToBlockMap_.put(fid, b);
      unkIDToReducedIDMap_[b].put(fid, k);
      unreducedUnkID_[b][k] = fid;
      k++;
    }
    SUNDANCE_MSG2(verb, tab2 << "block=" << b << " unk functions are " 
      << unreducedUnkID_[b]);
  }


  
  
  /* set up func ID maps for unk parameters */
  unreducedUnkParamID_.resize(unkParams.size());
  for (int i=0; i<unkParams.size(); i++)
  {
    const UnknownParameterElement* u 
      = dynamic_cast<const UnknownParameterElement*>(unkParams[i].ptr().get());
    TEUCHOS_TEST_FOR_EXCEPTION(u==0, std::runtime_error, 
      "EquationSet ctor input unk parameter "
      << unkParams[i] 
      << " does not appear to be a unk parameter");
    int fid = u->fid().dofID();
    unkParamSet_.put(fid);
    unkParamIDToReducedUnkParamIDMap_.put(fid, i);
    unreducedUnkParamID_[i] = fid;
  }
  SUNDANCE_MSG2(verb, tab1 << "unk parameters are " 
    << unreducedUnkParamID_);

  
  
  /* set up func ID maps for fixed parameters */
  unreducedFixedParamID_.resize(fixedParams.size());
  for (int i=0; i<fixedParams.size(); i++)
  {
    const UnknownParameterElement* u 
      = dynamic_cast<const UnknownParameterElement*>(fixedParams[i].ptr().get());
    TEUCHOS_TEST_FOR_EXCEPTION(u==0, std::runtime_error, 
      "EquationSet ctor input fixed parameter "
      << fixedParams[i] 
      << " does not appear to be a fixed parameter");
    int fid = u->fid().dofID();
    fixedParamSet_.put(fid);
    fixedParamIDToReducedFixedParamIDMap_.put(fid, i);
    unreducedFixedParamID_[i] = fid;
  }
  SUNDANCE_MSG2(verb, tab1 << "fixed parameters are " 
    << unreducedFixedParamID_);

  Set<OrderedHandle<CellFilterStub> > regionSet;

  /* Do the non-bc eqns first */
  SUNDANCE_MSG1(verb, tab1 << "processing integral terms");
  {
    Tabs tab2;
    SUNDANCE_MSG3(verb, tab2 << integralSum_->rqcToExprMap());
  }
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         r=integralSum_->rqcToExprMap().begin(); 
       r!=integralSum_->rqcToExprMap().end(); r++)
  {
    Tabs tab15;
    Tabs tab2;
    RegionQuadCombo rqc = r->first;
    int rqcVerb = verb;
    if (rqc.watch().isActive()) 
    {
      rqcVerb=rqc.watch().param("equation set setup");
      SUNDANCE_MSG1(max(verb, rqcVerb), tab15 << "processing RQC = " << rqc);
    }

    Expr term = r->second;

    SUNDANCE_MSG2(rqcVerb, tab2 << "expr = " << term);
    OrderedHandle<CellFilterStub> reg = rqc.domain();

    if (!regionSet.contains(reg)) 
    {
      regionSet.put(reg);
      Set<int> vf = integralSum_->funcsOnRegion(reg, varFuncSet_);
      Set<int> uf = integralSum_->funcsOnRegion(reg, unkFuncSet_);
      SUNDANCE_MSG2(rqcVerb, tab2 << "vf = " << vf);
      SUNDANCE_MSG2(rqcVerb, tab2 << "uf = " << uf);
      varsOnRegions_.put(reg, vf);
      unksOnRegions_.put(reg, uf);
    }
    else
    {
      Set<int>& t = varsOnRegions_.get(reg);
      t.merge(integralSum_->funcsOnRegion(reg, varFuncSet_));
      Set<int>& u = unksOnRegions_.get(reg);
      u.merge(integralSum_->funcsOnRegion(reg, unkFuncSet_));
    }
    SUNDANCE_MSG2(rqcVerb, tab2 << "varsOnRegion = " << varsOnRegions_.get(reg));
    SUNDANCE_MSG2(rqcVerb, tab2 << "unksOnRegion = " << unksOnRegions_.get(reg));
  }
  

  /* now do the BCs */
  if (hasBCs)
  {
    /* functions found in the BCs both in the overall lists and 
     * also in the bc-specific lists */
    for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
           r=bcSum_->rqcToExprMap().begin(); 
         r!=bcSum_->rqcToExprMap().end(); r++)
    {
      Tabs tab15;
      RegionQuadCombo rqc = r->first;
      int rqcVerb = verb;
      if (rqc.watch().isActive()) 
      {
        rqcVerb=rqc.watch().param("equation set setup");
        SUNDANCE_MSG1(max(verb, rqcVerb), tab15 << "processing RQC = " << rqc);
      }

      Expr term = r->second;
      OrderedHandle<CellFilterStub> reg = rqc.domain();

      if (!regionSet.contains(reg)) 
      {
        regionSet.put(reg);
        varsOnRegions_.put(reg, bcSum_->funcsOnRegion(reg, varFuncSet_));
        unksOnRegions_.put(reg, bcSum_->funcsOnRegion(reg, unkFuncSet_));
        bcVarsOnRegions_.put(reg, bcSum_->funcsOnRegion(reg, varFuncSet_));
        bcUnksOnRegions_.put(reg, bcSum_->funcsOnRegion(reg, unkFuncSet_));
      }
      else
      {
        if (!bcVarsOnRegions_.containsKey(reg))
        {
          bcVarsOnRegions_.put(reg, bcSum_->funcsOnRegion(reg, varFuncSet_));
        }
        if (!bcUnksOnRegions_.containsKey(reg))
        {
          bcUnksOnRegions_.put(reg, bcSum_->funcsOnRegion(reg, unkFuncSet_));
        }
        Set<int>& t = varsOnRegions_.get(reg);
        t.merge(bcSum_->funcsOnRegion(reg, varFuncSet_));
        Set<int>& u = unksOnRegions_.get(reg);
        u.merge(bcSum_->funcsOnRegion(reg, unkFuncSet_));
      }
    }
  }

  /* convert sets to arrays */
  regions_ = regionSet.elements();

  reducedVarsOnRegions_.resize(regions_.size());
  reducedUnksOnRegions_.resize(regions_.size());
  for (int r=0; r<regions_.size(); r++)
  {
    Tabs tab1;
    regionToIndexMap_.put(regions_[r], r);
    reducedVarsOnRegions_[r].resize(numVarBlocks());
    reducedUnksOnRegions_[r].resize(numUnkBlocks());
    OrderedHandle<CellFilterStub> cf = regions_[r];
    const Set<int>& v = this->varsOnRegion(r);
    const Set<int>& u = this->unksOnRegion(r);
    const Set<int>& bv = this->varsOnRegion(r);
    const Set<int>& bu = this->unksOnRegion(r);
    Set<int> vf = v;
    Set<int> uf = u;
    vf.merge(bv);
    uf.merge(bu);
    for (Set<int>::const_iterator i=vf.begin(); i!=vf.end(); i++)
    {
      int fid = *i;
      SUNDANCE_MSG2(verb, tab1 << "adding VF=" << fid << " to testToRegionMap");
      if (testToRegionsMap_.containsKey(fid))
      {
        testToRegionsMap_[fid].put(cf);
      }
      else
      {
        Set<OrderedHandle<CellFilterStub> > s;
        s.put(cf);
        testToRegionsMap_.put(fid, s);
      }
      int rv = reducedVarID(fid);
      int br = blockForVarID(fid);
      reducedVarsOnRegions_[r][br].put(rv);
    }
    for (Set<int>::const_iterator i=uf.begin(); i!=uf.end(); i++)
    {
      int fid = *i;
      if (unkToRegionsMap_.containsKey(fid))
      {
        unkToRegionsMap_[fid].put(cf);
      }
      else
      {
        Set<OrderedHandle<CellFilterStub> > s;
        s.put(cf);
        unkToRegionsMap_.put(fid, s);
      }
      int ru = reducedUnkID(fid);
      int bc = blockForUnkID(fid);
      reducedUnksOnRegions_[r][bc].put(ru);
    }
  }
}




Array<Expr> FunctionSupportResolver::flattenSpectral(const Array<Expr>& expr) const
{
  Array<Expr> rtn(expr.size());
  for (int i=0; i<expr.size(); i++)
  {
    const Expr& e = expr[i];
    rtn[i] = flattenSpectral(e);
  }
  return rtn;
}

Expr FunctionSupportResolver::flattenSpectral(const Expr& expr) const
{
  Array<Expr> rtn(expr.size());
  for (int i=0; i<expr.size(); i++)
  {
    if (expr[i].size() == 1)
    {
      const SpectralExpr* se 
        = dynamic_cast<const SpectralExpr*>(expr[i][0].ptr().get());
      if (se != 0)
      {
        int nt = se->getSpectralBasis().nterms();
        Array<Expr> e(nt);
        for (int j=0; j<nt; j++)
        {
          e[j] = se->getCoeff(j);
        }
        rtn[i] = new ListExpr(e);
      }
      else
      {
        rtn[i] = expr[i];
      }
    }
    else
    {
      rtn[i] = flattenSpectral(expr[i]);
    }
  }
  Expr r = new ListExpr(rtn);
  return r.flatten();
                  
}

int FunctionSupportResolver::reducedVarID(int varID) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!hasVarID(varID), std::runtime_error, 
    "varID " << varID << " not found in equation set");

  int b = blockForVarID(varID);
  return varIDToReducedIDMap_[b].get(varID);
}

int FunctionSupportResolver::reducedUnkID(int unkID) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!hasUnkID(unkID), std::runtime_error, 
    "unkID " << unkID << " not found in equation set");

  int b = blockForUnkID(unkID);
  return unkIDToReducedIDMap_[b].get(unkID);
}


int FunctionSupportResolver::reducedUnkParamID(int unkID) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!hasUnkParamID(unkID), std::runtime_error, 
    "unkParamID " << unkID << " not found in equation set");

  return unkParamIDToReducedUnkParamIDMap_.get(unkID);
}

int FunctionSupportResolver::reducedFixedParamID(int parID) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!hasFixedParamID(parID), std::runtime_error, 
    "fixedParamID " << parID << " not found in equation set");

  return fixedParamIDToReducedFixedParamIDMap_.get(parID);
}

int FunctionSupportResolver::blockForVarID(int varID) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!varIDToBlockMap_.containsKey(varID), std::runtime_error,
    "key " << varID << " not found in map "
    << varIDToBlockMap_);
  return varIDToBlockMap_.get(varID);
}

int FunctionSupportResolver::blockForUnkID(int unkID) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!unkIDToBlockMap_.containsKey(unkID), std::runtime_error,
    "key " << unkID << " not found in map "
    << unkIDToBlockMap_);
  return unkIDToBlockMap_.get(unkID);
}

const Set<OrderedHandle<CellFilterStub> >&  FunctionSupportResolver::regionsForTestFunc(int testID) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!testToRegionsMap_.containsKey(testID), std::runtime_error,
    "key " << testID << " not found in map "
    << testToRegionsMap_);
  return testToRegionsMap_.get(testID);
}

const Set<OrderedHandle<CellFilterStub> >&  FunctionSupportResolver::regionsForUnkFunc(int unkID) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!unkToRegionsMap_.containsKey(unkID), std::runtime_error,
    "key " << unkID << " not found in map "
    << testToRegionsMap_);
  return unkToRegionsMap_.get(unkID);
}

int FunctionSupportResolver::indexForRegion(const OrderedHandle<CellFilterStub>& region) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!regionToIndexMap_.containsKey(region), std::runtime_error,
    "key " << region << " not found in map "
    << regionToIndexMap_);
  return regionToIndexMap_.get(region);
}

bool FunctionSupportResolver::hasBCs() const 
{
  return bcSum_ != 0;
}


namespace Sundance
{

RCP<const CommonFuncDataStub> getSharedFunctionData(const FuncElementBase* f)
{
  const UnknownFuncElement* u = dynamic_cast<const UnknownFuncElement*>(f);
  const TestFuncElement* t = dynamic_cast<const TestFuncElement*>(f);

  if (u != 0) 
  {
    return rcp_dynamic_cast<const CommonFuncDataStub>(u->commonData());
  }
  if (t != 0)
  {
    return rcp_dynamic_cast<const CommonFuncDataStub>(t->commonData());
  }
  
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, 
    "unrecognized function type: " << typeid(*f).name());
  return u->commonData(); // -Wall, will never be called;
}
}

