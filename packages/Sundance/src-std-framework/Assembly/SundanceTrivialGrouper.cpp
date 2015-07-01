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

#include "SundanceTrivialGrouper.hpp"
#include "SundanceRefIntegral.hpp"
#include "SundanceQuadratureIntegral.hpp"
#include "SundanceReducedIntegral.hpp"
#include "SundanceMaximalQuadratureIntegral.hpp"
#include "SundanceCurveQuadratureIntegral.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceIntegralGroup.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceReducedQuadrature.hpp"
#include "SundanceMap.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;
using namespace Teuchos;



void TrivialGrouper::findGroups(const EquationSet& eqn,
  const CellType& maxCellType,
  int spatialDim,
  const CellType& cellType,
  int cellDim,
  const QuadratureFamily& quad,
  const RCP<SparsitySuperset>& sparsity,
  bool isInternalBdry,
  Array<RCP<IntegralGroup> >& groups,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh) const
{
  Tabs tab(0);

  SUNDANCE_MSG1(setupVerb(),
    tab << "in TrivialGrouper::findGroups(), num derivs = " 
    << sparsity->numDerivs());
  SUNDANCE_MSG1(setupVerb(), 
    tab << "cell type = " << cellType);
  SUNDANCE_MSG1(setupVerb(), 
    tab << "sparsity = " << std::endl << *sparsity << std::endl);

  const ReducedQuadrature* rq = dynamic_cast<const ReducedQuadrature*>(quad.ptr(
).get());
  bool useReducedQuad = (rq != 0);
  SUNDANCE_MSG1(setupVerb(), tab << "using reduced quadrature: " 
    << useReducedQuad);


  int vecCount=0;
  int constCount=0;

  bool isMaximal = cellType == maxCellType;
  bool useMaxIntegral = isMaximal;
  
  /* turn off grouping for submaximal cells. This works around 
   * a bug detected by Rob Kirby that
   * shows up with Nitsche BCs in mixed-element discretizations */
  bool doGroups = true;
  if (!isMaximal) doGroups = false;


  typedef Sundance::Map<OrderedQuartet<int, BasisFamily, int, BasisFamily>, Array<RCP<ElementIntegral> > > twoFormMap;
  typedef Sundance::Map<OrderedTriple<int,int,BasisFamily>, Array<RCP<ElementIntegral> > > oneFormMap;
  Sundance::Map<OrderedQuartet<int, BasisFamily, int, BasisFamily>, Array<RCP<ElementIntegral> > > twoForms;
  Sundance::Map<OrderedQuartet<int, BasisFamily, int, BasisFamily>, Array<int> > twoFormResultIndices;
  Sundance::Map<OrderedTriple<int,int,BasisFamily>, Array<RCP<ElementIntegral> > > oneForms;
  Sundance::Map<OrderedTriple<int,int,BasisFamily>, Array<int> > oneFormResultIndices;

  for (int i=0; i<sparsity->numDerivs(); i++)
  {
    const MultipleDeriv& d = sparsity->deriv(i);
    SUNDANCE_MSG3(setupVerb(),
      tab << "--------------------------------------------------");
    SUNDANCE_MSG3(setupVerb(),
      tab << "defining integration policy for " << d);
    SUNDANCE_MSG3(setupVerb(),
      tab << "--------------------------------------------------");
      
    if (d.order()==0) 
    {
      RCP<ElementIntegral> integral ;
      int resultIndex;
      if (sparsity->isConstant(i))
      {
        if (globalCurve.isCurveValid() && globalCurve.isCurveIntegral() && isMaximal )
        { // ----- curve Integral ------
          integral = rcp(new CurveQuadratureIntegral( maxCellType, true ,
              quad, globalCurve , mesh , setupVerb() ) );
        }
        else
        { // --- no curve integral ---
          integral = rcp(new RefIntegral(spatialDim, maxCellType,
              cellDim, cellType, quad , isInternalBdry, globalCurve , mesh , setupVerb()));
        }
        resultIndex = constCount++;
      }
      else
      {
        if (useReducedQuad)
        { 
          integral = rcp(new ReducedIntegral(spatialDim, maxCellType,
              cellDim, cellType, quad , isInternalBdry, globalCurve , mesh , setupVerb()));
        }
        else if (useMaxIntegral)
        {
          if (globalCurve.isCurveValid() && globalCurve.isCurveIntegral() && isMaximal )
          { // ----- curve Integral ------
            integral = rcp(new CurveQuadratureIntegral(maxCellType, false ,
                quad, globalCurve , mesh , setupVerb()));
          }
          else
          { // --- no curve integral ---
            integral = rcp(new MaximalQuadratureIntegral(maxCellType,
                quad, globalCurve , mesh , setupVerb()));
          }
        }
        else // no maxCell Integral
        {
          integral = rcp(new QuadratureIntegral(spatialDim, maxCellType, 
              cellDim, cellType, quad, isInternalBdry, globalCurve , mesh ,
              setupVerb()));
        }
        resultIndex = vecCount++;
      }
      integral->setVerb(integrationVerb(), transformVerb());
      SUNDANCE_MSG3(setupVerb(), tab << "is zero-form");
      groups.append(rcp(new IntegralGroup(tuple(integral),
            tuple(resultIndex), setupVerb())));
    }
    else
    {
      BasisFamily testBasis;
      BasisFamily unkBasis;
      MultiIndex miTest;
      MultiIndex miUnk;
      int rawTestID = -1;
      int rawUnkID = -1;
      int rawParamID = -1;
      int testID = -1;
      int unkID = -1;
      int paramID = -1;
      int testBlock = -1;
      int unkBlock = -1;
      bool isOneForm;
      bool hasParam;
      extractWeakForm(eqn, d, testBasis, unkBasis, miTest, miUnk, 
        rawTestID, rawUnkID,
        testID, unkID,
        testBlock, unkBlock,
        rawParamID, paramID,
        isOneForm, hasParam);
      
      TEUCHOS_TEST_FOR_EXCEPT(hasParam && !isOneForm);

      /* The parameter index acts as an index into a multivector. If
       * this one-form is not a parametric derivative, we use zero as
       * the multivector index */
      int mvIndex = 0;
      if (hasParam) mvIndex = paramID; 
      
      /* In variational problems we might have (u,v) and (v,u). Because 
       * the derivative is stored as an unordered multiset it can't 
       * distinguish between the two cases. We need to check the equation 
       * set to see if the two functions show up as variations and 
       * unknowns. If so, then we need to produce the transposed integral.
       */
      bool transposeNeeded = false;
      if (!isOneForm && rawTestID!=rawUnkID 
        && eqn.hasVarID(rawUnkID) && eqn.hasUnkID(rawTestID))
      {
        transposeNeeded = true;
      }


      if (isOneForm)
      {
        SUNDANCE_MSG3(setupVerb(), tab << "is one-form");
      }
      else
      {
        SUNDANCE_MSG3(setupVerb(), tab << "is two-form");
      }


      if (hasParam)
      {
        SUNDANCE_MSG3(setupVerb(), tab << "is a parametric derivative");
      }
      else
      {
        SUNDANCE_MSG3(setupVerb(), tab << "is not a parametric derivative");
      }

      SUNDANCE_MSG3(setupVerb(), 
        tab << "test ID: " << testID << " block=" << testBlock);

      if (!isOneForm)
      {
        SUNDANCE_MSG3(setupVerb(), tab << "unk funcID: " << unkID << " block=" << unkBlock);
      }

      if (hasParam)
      {
        SUNDANCE_MSG3(setupVerb(), tab << "paramID=" << paramID);
      }
                   
      SUNDANCE_MSG3(setupVerb(), tab << "deriv = " << d);
      if (sparsity->isConstant(i))
      {
        SUNDANCE_MSG3(setupVerb(), tab << "coeff is constant");
      }
      else
      {
        SUNDANCE_MSG3(setupVerb(), tab << "coeff is non-constant");
      }

      RCP<ElementIntegral> integral;
      RCP<ElementIntegral> transposedIntegral;
      int resultIndex;
      if (sparsity->isConstant(i))
      {
        if (isOneForm)
        {
          int alpha=0;
          if (miTest.order()==1)
          {
            alpha = miTest.firstOrderDirection();
          }
          // test if we need to make curve Integral
          if (globalCurve.isCurveValid() && globalCurve.isCurveIntegral() && isMaximal )
          { // ----- curve Integral ------
            integral = rcp(new CurveQuadratureIntegral(maxCellType, true ,
                testBasis, alpha,
                miTest.order(), quad, globalCurve , mesh ,setupVerb()));
          }
          else
          { // --- no curve integral ---
            integral = rcp(new RefIntegral(spatialDim, maxCellType,
                cellDim, cellType,
                testBasis, alpha,
                miTest.order(), quad , isInternalBdry, globalCurve , mesh ,setupVerb()));
          }
        }
        else
        {
          int alpha=0;
          int beta=0;
          if (miTest.order()==1)
          {
            alpha = miTest.firstOrderDirection();
          }
          if (miUnk.order()==1)
          {
            beta = miUnk.firstOrderDirection();
          }
          // test if we need to make curve Integral
          if (globalCurve.isCurveValid() && globalCurve.isCurveIntegral() && isMaximal )
          { // ----- curve Integral ------
            integral = rcp(new CurveQuadratureIntegral(maxCellType, true ,
                testBasis, alpha,
                miTest.order(),
                unkBasis, beta,
                miUnk.order(), quad, globalCurve , mesh ,setupVerb()));
          }
          else // --- no curve integral ---
          {
            integral = rcp(new RefIntegral(spatialDim, maxCellType,
                cellDim, cellType,
                testBasis, alpha, miTest.order(),
                unkBasis, beta, miUnk.order(), quad , isInternalBdry, globalCurve , mesh ,setupVerb()));
          }
          if (transposeNeeded)
          {
            if (globalCurve.isCurveValid() && globalCurve.isCurveIntegral() && isMaximal )
            { // ----- curve Integral ------
              transposedIntegral = rcp(new CurveQuadratureIntegral(maxCellType, true ,
                  unkBasis, beta,
                  miUnk.order(),
                  testBasis, alpha,
                  miTest.order(),
                  quad, globalCurve , mesh ,setupVerb()));
            }
            else // --- no curve integral ---
            {
              transposedIntegral = rcp(new RefIntegral(spatialDim, maxCellType,
                  cellDim, cellType,
                  unkBasis, beta,
                  miUnk.order(),
                  testBasis, alpha,
                  miTest.order(), quad , isInternalBdry, globalCurve , mesh ,setupVerb()));
            }
          }
        }
        resultIndex = constCount++;
      }
      else /* sparsity->isVector(i) */
      {
        if (isOneForm)
        {
          int alpha=0;
          if (miTest.order()==1)
          {
            alpha = miTest.firstOrderDirection();
          }
          if (useReducedQuad)
          {
            integral = rcp(new ReducedIntegral(spatialDim, maxCellType,
                cellDim, cellType,
                testBasis, alpha,
                miTest.order(), quad , isInternalBdry, globalCurve , mesh ,setupVerb()));
          }
          else if (useMaxIntegral)
          {
            if ( globalCurve.isCurveValid() && globalCurve.isCurveIntegral() )
            { // ----- curve Integral ------
              integral = rcp(new CurveQuadratureIntegral(maxCellType, false ,
                  testBasis, alpha,
                  miTest.order(), quad, globalCurve , mesh ,setupVerb()));
            }
            else
            {// --- no curve integral ---
              integral = rcp(new MaximalQuadratureIntegral(maxCellType,
                  testBasis, alpha,
                  miTest.order(), quad, globalCurve , mesh ,setupVerb()));
            }
          }
          else // no maxCell Integral
          {
            integral = rcp(new QuadratureIntegral(spatialDim, maxCellType,
                cellDim, cellType,
                testBasis, alpha, 
                miTest.order(), quad, isInternalBdry, globalCurve , mesh ,setupVerb()));
          }
        }
        else /* two-form */
        {
          int alpha=0;
          int beta=0;
          if (miTest.order()==1)
          {
            alpha = miTest.firstOrderDirection();
          }
          if (miUnk.order()==1)
          {
            beta = miUnk.firstOrderDirection();
          }
          if (useReducedQuad)
          {
            integral = rcp(new ReducedIntegral(spatialDim, maxCellType,
                cellDim, cellType,
                testBasis, alpha, miTest.order(),
                unkBasis, beta, miUnk.order(), quad , isInternalBdry, globalCurve , mesh ,setupVerb()));
            if (transposeNeeded)
            {
              transposedIntegral = rcp(new ReducedIntegral(spatialDim, maxCellType,
                  cellDim, cellType,
                  unkBasis, beta,
                  miUnk.order(),
                  testBasis, alpha,
                  miTest.order(), quad , isInternalBdry, globalCurve , mesh ,setupVerb()));
            }
          }
          else if (useMaxIntegral)
          {
            if ( globalCurve.isCurveValid() && globalCurve.isCurveIntegral() )
            { // ----- curve Integral ------
              integral = rcp(new CurveQuadratureIntegral(maxCellType, false ,
                  testBasis, alpha,
                  miTest.order(),
                  unkBasis, beta,
                  miUnk.order(), quad, globalCurve , mesh ,setupVerb()));
            }
            else
            {// --- no curve integral ---
              integral = rcp(new MaximalQuadratureIntegral(maxCellType,
                  testBasis, alpha,
                  miTest.order(),
                  unkBasis, beta,
                  miUnk.order(), quad, globalCurve , mesh ,setupVerb()));
            }
            if (transposeNeeded)
            {
              if ( globalCurve.isCurveValid() && globalCurve.isCurveIntegral() )
              { // ----- curve Integral ------
                transposedIntegral = rcp(new CurveQuadratureIntegral(maxCellType, false ,
                    unkBasis, beta, miUnk.order(),
                    testBasis, alpha, miTest.order(), quad,
                    globalCurve , mesh ,setupVerb()));
              }
              else
              { // --- no curve integral ---
                transposedIntegral = rcp(new MaximalQuadratureIntegral(maxCellType,
                    unkBasis, beta, miUnk.order(),
                    testBasis, alpha, miTest.order(), quad,
                    globalCurve , mesh ,setupVerb()));
              }
            }
          }
          else // no MaxCell integral
          {
            integral = rcp(new QuadratureIntegral(spatialDim, maxCellType,
                cellDim, cellType,
                testBasis, alpha, 
                miTest.order(),
                unkBasis, beta, 
                miUnk.order(), quad, isInternalBdry, globalCurve , mesh ,setupVerb()));
            if (transposeNeeded)
            {
              transposedIntegral = rcp(new QuadratureIntegral(spatialDim, maxCellType,
                  cellDim, cellType,
                  unkBasis, beta, miUnk.order(),
                  testBasis, alpha, miTest.order(), quad, 
                  isInternalBdry, globalCurve , mesh ,setupVerb()));
            }
          }
        }
        resultIndex = vecCount++;
      }

      /* Set the verbosity for the integrals */
      integral->setVerb(integrationVerb(), transformVerb());
      if (transposeNeeded)
      {
        transposedIntegral->setVerb(integrationVerb(), transformVerb());
      }
          
      
      if (isOneForm)
      {
        if (doGroups)
        {
          OrderedTriple<int,int,BasisFamily> testKey(rawTestID, mvIndex, testBasis);
          if (!oneForms.containsKey(testKey))
          {
            oneForms.put(testKey, tuple(integral));
            oneFormResultIndices.put(testKey, tuple(resultIndex));
          }
          else
          {
            oneForms[testKey].append(integral);
            oneFormResultIndices[testKey].append(resultIndex);
          }
        }
        else
        {
          groups.append(rcp(new IntegralGroup(tuple(testID), tuple(testBlock),
                tuple(mvIndex),
                tuple(integral),
                tuple(resultIndex), tuple(d), setupVerb())));
        }
      }
      else
      {
        if (!doGroups)
        {
          groups.append(rcp(new IntegralGroup(tuple(testID), tuple(testBlock),
                tuple(unkID), tuple(unkBlock),
                tuple(integral),
                tuple(resultIndex), tuple(d), setupVerb())));
          if (transposeNeeded)
          {
            groups.append(rcp(new IntegralGroup(tuple(unkID), tuple(unkBlock),
                  tuple(testID), tuple(testBlock),
                  tuple(transposedIntegral),
                  tuple(resultIndex), tuple(d), setupVerb())));
          }
        }
        else
        {
          Tabs tab3;
          OrderedQuartet<int, BasisFamily, int, BasisFamily> testUnkKey(rawTestID, testBasis, rawUnkID, unkBasis);


          SUNDANCE_MSG2(setupVerb(), tab3 << "key=" << testUnkKey);
          if (!twoForms.containsKey(testUnkKey))
          {
            Tabs tab4;
            SUNDANCE_MSG2(setupVerb(), tab4 << "key not found");
            twoForms.put(testUnkKey, tuple(integral));
            twoFormResultIndices.put(testUnkKey, tuple(resultIndex));
          }
          else
          {
            Tabs tab4;
            SUNDANCE_MSG2(setupVerb(), tab4 << "key found");
            twoForms[testUnkKey].append(integral);
            twoFormResultIndices[testUnkKey].append(resultIndex);
          }
          if (transposeNeeded)
          {
            OrderedQuartet<int, BasisFamily, int, BasisFamily> unkTestKey(rawUnkID, unkBasis, rawTestID, testBasis);
            
            if (!twoForms.containsKey(unkTestKey))
            {
              Tabs tab4;
              SUNDANCE_MSG2(setupVerb(), tab4 << "key not found");
              twoForms.put(unkTestKey, tuple(transposedIntegral));
              twoFormResultIndices.put(unkTestKey, tuple(resultIndex));
            }
            else
            {
              Tabs tab4;
              SUNDANCE_MSG2(setupVerb(), tab4 << "key found");
              twoForms[unkTestKey].append(transposedIntegral);
              twoFormResultIndices[unkTestKey].append(resultIndex);
            }
          }
        }
      }
    }
  }

  if (doGroups)
  {
    Tabs tab;
    SUNDANCE_MSG2(setupVerb(), tab << "creating integral groups");
    for (twoFormMap::const_iterator i=twoForms.begin(); i!=twoForms.end(); i++)
    {
      Tabs tab3;
      SUNDANCE_MSG2(setupVerb(), tab3 << "integral group number="
        << groups.size());
      int rawTestID = i->first.a();
      BasisFamily testBasis = i->first.b();
      int rawUnkID = i->first.c();
      BasisFamily unkBasis = i->first.d();
      int testID = eqn.reducedVarID(rawTestID);
      int unkID = eqn.reducedUnkID(rawUnkID);
      int testBlock = eqn.blockForVarID(rawTestID);
      int unkBlock = eqn.blockForUnkID(rawUnkID);
      const Array<RCP<ElementIntegral> >& integrals = i->second;
      const Array<int>& resultIndices 
        = twoFormResultIndices.get(i->first);
      SUNDANCE_MSG2(setupVerb(), tab3 << "creating two-form integral group" << std::endl
        << tab3 << "testID=" << rawTestID << std::endl
        << tab3 << "unkID=" << rawUnkID << std::endl
        << tab3 << "testBlock=" << testBlock << std::endl
        << tab3 << "unkBlock=" << unkBlock << std::endl
        << tab3 << "testBasis=" << testBasis << std::endl
        << tab3 << "unkBasis=" << unkBasis << std::endl
        << tab3 << "resultIndices=" << resultIndices);
      Array<MultipleDeriv> grpDerivs;
      for (int j=0; j<resultIndices.size(); j++)
      {
        MultipleDeriv d = sparsity->deriv(resultIndices[j]);
        SUNDANCE_MSG2(setupVerb(), tab3 << "deriv " << j << " " 
          << d);
        grpDerivs.append(d);
      }
      groups.append(rcp(new IntegralGroup(tuple(testID), tuple(testBlock), 
            tuple(unkID), tuple(unkBlock),
            integrals, resultIndices, grpDerivs, setupVerb())));
    }

    for (oneFormMap::const_iterator i=oneForms.begin(); i!=oneForms.end(); i++)
    {
      Tabs tab3;
      SUNDANCE_MSG2(setupVerb(), tab3 << "integral group number="
        << groups.size());
      int rawTestID = i->first.a();
      int mvIndex = i->first.b();
      int testID = eqn.reducedVarID(rawTestID);
      int testBlock = eqn.blockForVarID(rawTestID);
      const Array<RCP<ElementIntegral> >& integrals = i->second;
      const Array<int>& resultIndices 
        = oneFormResultIndices.get(i->first);
      SUNDANCE_MSG2(setupVerb(), tab3 << "creating one-form integral group" << std::endl
        << tab3 << "testID=" << testID << std::endl
        << tab3 << "resultIndices=" << resultIndices);
      Array<MultipleDeriv> grpDerivs;
      for (int j=0; j<resultIndices.size(); j++)
      {
        MultipleDeriv d = sparsity->deriv(resultIndices[j]);
        SUNDANCE_MSG2(setupVerb(), tab3 << "deriv " << j << " " 
          << d);
        grpDerivs.append(d);
      }
      groups.append(rcp(new IntegralGroup(tuple(testID), tuple(testBlock),
            tuple(mvIndex),
            integrals, resultIndices, grpDerivs, setupVerb())));
    }
  }
  
  
}

