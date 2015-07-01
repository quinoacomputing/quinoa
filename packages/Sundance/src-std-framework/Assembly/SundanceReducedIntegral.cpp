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

#include "SundanceReducedIntegral.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceGaussianQuadratureType.hpp"
#include "SundanceQuadratureType.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Teuchos;

using std::ios_base;
using std::setw;
using std::endl;

extern "C" 
{
int dgemm_(const char* transA, const char* transB,
  const int* M, const int *N, const int* K,
  const double* alpha, 
  const double* A, const int* ldA,
  const double* B, const int* ldB,
  const double* beta,
  double* C, const int* ldC);
}

static Time& reduced0IntegrationTimer() 
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("reduced 0-form integration"); 
  return *rtn;
}


static Time& reduced1IntegrationTimer() 
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("reduced 1-form integration"); 
  return *rtn;
}


static Time& reduced2IntegrationTimer() 
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("reduced 2-form integration"); 
  return *rtn;
}



ReducedIntegral::ReducedIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim, 
  const CellType& cellType,
  const QuadratureFamily& quad_in,
  bool isInternalBdry,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh,
  int verb)
  : ElementIntegral(spatialDim, maxCellType, dim, cellType, isInternalBdry, globalCurve , mesh,
    verb), W_()
{
  Tabs tab0(0);

  SUNDANCE_MSG1(setupVerb(),
    tab0 << "************* creating reduced 0-form integrals ********");
  if (setupVerb()) describe(Out::os());
  
  /* we need to sum the quadrature weights 
     to compute the volume of the reference cell */
  QuadratureFamily quad = new GaussianQuadrature(1);

  /* not supporting ACI for now */
  TEUCHOS_TEST_FOR_EXCEPT(globalCurve.isCurveValid());

  Array<Point> quadPts;
  Array<double> quadWeights;

  W_.resize(1);
  W_[0].resize(1);

  quad.getPoints(cellType, quadPts, quadWeights);  

  for (int q=0; q<quadWeights.size(); q++) {
	  W_[0][0] += quadWeights[q];
  }
}

ReducedIntegral::ReducedIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim, 
  const CellType& cellType,
  const BasisFamily& testBasis,
  int alpha,
  int testDerivOrder,
  const QuadratureFamily& quad_in,
  bool isInternalBdry,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh,
  int verb)
  : ElementIntegral(spatialDim, maxCellType, dim, cellType, 
    testBasis, alpha, testDerivOrder, isInternalBdry, globalCurve , mesh , verb), W_()
{
  Tabs tab0(0);
  SUNDANCE_MSG1(setupVerb(),
    tab0 << "************* creating reduced 1-form integrals ********");
  if (setupVerb()) describe(Out::os());
  assertLinearForm();

  W_.resize(nFacetCases());

  /* Determine the quadrature order needed for exact integrations */
  QuadratureType qType = new GaussianQuadratureType();
  int reqOrder = qType.findValidOrder(cellType, 
    std::max(1, testBasis.order()));

  SUNDANCE_MSG2(setupVerb(),
    tab0 << "using quadrature order=" << reqOrder);
   
  /* Create a quadrature family of the required order */
  QuadratureFamily quad = qType.createQuadFamily(reqOrder);

  /* not supporting ACI for now */
  TEUCHOS_TEST_FOR_EXCEPT(globalCurve.isCurveValid());

  /* We now loop over the different evaluation cases, integrating the
   * basis functions for each. Because this is a reference integral,
   * we can actually do the untransformed integrals here. */
  for (int fc=0; fc<nFacetCases(); fc++)
  {
    Tabs tab1;
    SUNDANCE_MSG2(setupVerb(),
      tab1 << "evaluation case=" << fc << " of " << nFacetCases());
    /* initialize size of untransformed integral results array */
    W_[fc].resize(nRefDerivTest() * nNodesTest());

    /* initialize values of integrals to zero */
    for (int i=0; i<W_[fc].size(); i++) { W_[fc][i]=0.0; }

    Array<Array<Array<Array<double> > > > testBasisVals(nRefDerivTest());
  
    /* get quadrature points */

    Array<Point> quadPts;
    Array<double> quadWeights;
    getQuad(quad, fc, quadPts, quadWeights);

    int nQuad = quadPts.size();

    /* compute the basis functions */
    for (int r=0; r<nRefDerivTest(); r++)
    {
      Tabs tab2;
      SUNDANCE_MSG2(setupVerb(), tab2 << "evaluating basis derivative " 
        << r << " of " << nRefDerivTest());

      testBasisVals[r].resize(testBasis.dim());
      MultiIndex mi;
      if (testDerivOrder==1) mi[r] = 1;
      SpatialDerivSpecifier deriv(mi);
      testBasis.refEval(evalCellType(), quadPts, deriv,
        testBasisVals[r], setupVerb());
    }

    /* do the quadrature */
    SUNDANCE_MSG2(setupVerb(), tab1 << "doing quadrature");
    int vecComp = 0;
    for (int q=0; q<nQuad; q++)
    {
      for (int t=0; t<nRefDerivTest(); t++)
      {
        for (int nt=0; nt<nNodesTest(); nt++)
        {
          value(fc, t, nt) 
            += chop(quadWeights[q] * testBasisVals[t][vecComp][q][nt]) ;
        }
      }
    }    

    for (int i=0; i<W_[fc].size(); i++) W_[fc][i] = chop(W_[fc][i]);

    addFlops(3*nQuad*nRefDerivTest()*nNodesTest() + W_[fc].size());
  }

  /* print the result */
  SUNDANCE_MSG4(setupVerb(), tab0 << "--------------------------------------");
  SUNDANCE_MSG4(setupVerb(), tab0 << "reduced linear form integral results");
  if (setupVerb() >= 4)
  {
    for (int fc=0; fc<nFacetCases(); fc++)
    {
      Tabs tab1;
      SUNDANCE_MSG4(setupVerb(), tab1 << "------ evaluation case " << fc << " of "
        << nFacetCases() << "-------");
      
      for (int r=0; r<nRefDerivTest(); r++)
      {
        Tabs tab2;

        MultiIndex mi;
        if (testDerivOrder==1) mi[r] = 1;
        SUNDANCE_MSG1(setupVerb(), tab2 << "multiindex=" << mi);

        ios_base::fmtflags oldFlags = Out::os().flags();
        Out::os().setf(ios_base::right);    
        Out::os().setf(ios_base::showpoint);
        for (int nt=0; nt<nNodesTest(); nt++)
        {
          Tabs tab3;
          Out::os() << tab3 << setw(10) << nt 
                    << setw(12) << std::setprecision(5) << value(fc, r, nt) 
                    << std::endl;
        }
        Out::os().flags(oldFlags);
      }
    }
  }

  SUNDANCE_MSG1(setupVerb(), tab0 << "done reduced linear form ctor");
}




ReducedIntegral::ReducedIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim,
  const CellType& cellType,
  const BasisFamily& testBasis,
  int alpha,
  int testDerivOrder,
  const BasisFamily& unkBasis,
  int beta,
  int unkDerivOrder, 
  const QuadratureFamily& quad_in,
  bool isInternalBdry,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh,
  int verb)
  : ElementIntegral(spatialDim, maxCellType,  dim, cellType,
    testBasis, alpha, testDerivOrder, 
    unkBasis, beta, unkDerivOrder, isInternalBdry, globalCurve , mesh ,verb), W_()

{
  Tabs tab0(0);
  SUNDANCE_MSG1(setupVerb(),
    tab0 << "************* creating reduced 2-form integrals ********");
  if (setupVerb()) describe(Out::os());

  assertBilinearForm();

  W_.resize(nFacetCases());

  QuadratureType qType = new GaussianQuadratureType();
  int reqOrder = qType.findValidOrder(cellType,
    std::max(1, unkBasis.order() + testBasis.order()));

  SUNDANCE_MSG2(setupVerb(),
    tab0 << "using quadrature order=" << reqOrder);
  QuadratureFamily quad = qType.createQuadFamily(reqOrder);

  /* not supporting ACI for now */
  TEUCHOS_TEST_FOR_EXCEPT(globalCurve.isCurveValid());


  SUNDANCE_MSG2(setupVerb(),
    tab0 << "processing evaluation cases");

  for (int fc=0; fc<nFacetCases(); fc++)
  {
    Tabs tab1;
    SUNDANCE_MSG1(setupVerb(), tab1 << "------ evaluation case " << fc << " of "
      << nFacetCases() << "-------");
    
    W_[fc].resize(nRefDerivTest() * nNodesTest()  * nRefDerivUnk() * nNodesUnk());
    for (int i=0; i<W_[fc].size(); i++) W_[fc][i]=0.0;

    Array<Array<Array<Array<double> > > > testBasisVals(nRefDerivTest());
    Array<Array<Array<Array<double> > > > unkBasisVals(nRefDerivUnk());

    Array<Point> quadPts;
    Array<double> quadWeights;
    getQuad(quad, fc, quadPts, quadWeights);

    int nQuad = quadPts.size();

    for (int r=0; r<nRefDerivTest(); r++)
    {
      Tabs tab2;
      SUNDANCE_MSG2(setupVerb(), tab2 
        << "evaluating test function basis derivative " 
        << r << " of " << nRefDerivTest());
      testBasisVals[r].resize(testBasis.dim());
      MultiIndex mi;
      if (testDerivOrder==1) mi[r] = 1;
      SpatialDerivSpecifier deriv(mi);
      testBasis.refEval(evalCellType(), quadPts, deriv,
        testBasisVals[r], setupVerb());
    }

    for (int r=0; r<nRefDerivUnk(); r++)
    {
      Tabs tab2;
      SUNDANCE_MSG2(setupVerb(), tab2 
        << "evaluating unknown function basis derivative " 
        << r << " of " << nRefDerivUnk());
      unkBasisVals[r].resize(unkBasis.dim());
      MultiIndex mi;
      if (unkDerivOrder==1) mi[r] = 1;
      SpatialDerivSpecifier deriv(mi);
      unkBasis.refEval(evalCellType(), 
        quadPts, deriv, unkBasisVals[r], setupVerb());
    }

    SUNDANCE_MSG2(setupVerb(), tab1 << "doing quadrature...");
    int vecComp = 0;
    for (int q=0; q<nQuad; q++)
    {
      for (int t=0; t<nRefDerivTest(); t++)
      {
        for (int nt=0; nt<nNodesTest(); nt++)
        {
          for (int u=0; u<nRefDerivUnk(); u++)
          {
            for (int nu=0; nu<nNodesUnk(); nu++)
            {
              value(fc, t, nt, u, nu) 
                += chop(quadWeights[q] * testBasisVals[t][vecComp][q][nt]
                  * unkBasisVals[u][vecComp][q][nu]);
            }
          }
        }
      }
    }
    SUNDANCE_MSG2(setupVerb(), tab1 << "...done");
    addFlops(4*nQuad*nRefDerivTest()*nNodesTest()*nRefDerivUnk()*nNodesUnk()
      + W_[fc].size());
    for (int i=0; i<W_[fc].size(); i++) W_[fc][i] = chop(W_[fc][i]);
  }

  SUNDANCE_MSG1(setupVerb(), tab0 
    << "----------------------------------------");
  SUNDANCE_MSG4(setupVerb(), tab0 
    << "reduced bilinear form integral results");
  if (setupVerb() >= 4)
  {
    for (int fc=0; fc<nFacetCases(); fc++)
    {
      Tabs tab1;
      SUNDANCE_MSG4(setupVerb(), tab1 << "evaluation case " << fc << " of "
        << nFacetCases());
      
      for (int rt=0; rt<nRefDerivTest(); rt++)
      {
        for (int ru=0; ru<nRefDerivUnk(); ru++)
        {
          Tabs tab2;
          MultiIndex miTest;
          if (testDerivOrder==1) miTest[rt] = 1;
          MultiIndex miUnk;
          if (unkDerivOrder==1) miUnk[ru] = 1;
          SUNDANCE_MSG1(setupVerb(), tab2 << "test multiindex=" << miTest
            << " unk multiindex=" << miUnk);
          
          ios_base::fmtflags oldFlags = Out::os().flags();
          Out::os().setf(ios_base::right);    
          Out::os().setf(ios_base::showpoint);
          for (int nt=0; nt<nNodesTest(); nt++)
          {
            Tabs tab3;
            Out::os() << tab3 << setw(10) << nt;
            for (int nu=0; nu<nNodesUnk(); nu++)
            {
              Out::os() << setw(12) << std::setprecision(5)
                        << value(fc, rt, nt, ru, nu) ;
            }
            Out::os() << std::endl;
          }
          Out::os().flags(oldFlags);
        }
      }
    }
  }

  SUNDANCE_MSG1(setupVerb(), tab0 << "done reduced bilinear form ctor");
}




void ReducedIntegral::transformZeroForm(const CellJacobianBatch& JTrans,
  const CellJacobianBatch& JVol,
  const Array<int>& isLocalFlag,  
  const Array<int>& facetIndex,
  const RCP<Array<int> >& cellLIDs,
  const double* const coeffs,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(reduced0IntegrationTimer());

  TEUCHOS_TEST_FOR_EXCEPTION(order() != 0, std::logic_error,
    "ReducedIntegral::transformZeroForm() called "
    "for form of order " << order());

  Tabs tabs;  
  SUNDANCE_MSG1(integrationVerb(), tabs << "doing zero form by reduced integration");

  double& a = (*A)[0];
  int flops = 0;

  /* if we don't need to check whether elements are local, we
   * can streamline the loop. This will be the case when
   * we are evaluating a functional but not its gradient */
  double w = W_[0][0];
  if ((int) isLocalFlag.size()==0)
  {
    /* not supporting ACI for now */
    if (globalCurve().isCurveValid())
    {
      TEUCHOS_TEST_FOR_EXCEPT(globalCurve().isCurveValid());
    }
    else
    {
      for (int c=0; c<JVol.numCells(); c++)
      {
 				flops+=3;
 				a += w * coeffs[c] * fabs(JVol.detJ()[c]);
      }
    }

  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION( (int) isLocalFlag.size() != JVol.numCells(),
      std::runtime_error,
      "mismatch between isLocalFlag.size()=" 
      << isLocalFlag.size()
      << " and JVol.numCells()=" << JVol.numCells());

    if (globalCurve().isCurveValid())
    {  
      TEUCHOS_TEST_FOR_EXCEPT(globalCurve().isCurveValid());
    }
    else         /* ---------- NO ACI logic ----------- */
    {
      for (int c=0; c<JVol.numCells(); c++)
      {
        if (isLocalFlag[c])
        {
    			flops+=3;
    			a += w * coeffs[c] * fabs(JVol.detJ()[c]);
        }
      }
    }
  }
  addFlops(flops);
}

void ReducedIntegral::transformOneForm(const CellJacobianBatch& JTrans,  
  const CellJacobianBatch& JVol,
  const Array<int>& facetIndex,
  const RCP<Array<int> >& cellLIDs,
  const double* const coeffs,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(reduced1IntegrationTimer());
  TEUCHOS_TEST_FOR_EXCEPTION(order() != 1, std::logic_error,
    "ReducedIntegral::transformOneForm() called for form "
    "of order " << order());
  Tabs tabs;  
  SUNDANCE_MSG1(integrationVerb(),
    tabs << "doing one form by reduced integration");

	int nfc = nFacetCases();

  /* If the derivative order is zero, the only transformation to be done 
   * is to multiply by the cell's Jacobian determinant.  Each entry also needs
   * to be multiplied by the coefficient for the current cell. */
  if (testDerivOrder() == 0)
  {
    double* aPtr = &((*A)[0]);
    int count = 0;
    if (globalCurve().isCurveValid())
    {    /* ---------- ACI logic ----------- */
      TEUCHOS_TEST_FOR_EXCEPT(globalCurve().isCurveValid());      
    }
    else         /* ---------- NO ACI logic ----------- */
    {
      for (int c=0; c<JVol.numCells(); c++)
      {
        double w_cell = coeffs[c] * fabs(JVol.detJ()[c]);
        int fc = 0;
        if (nfc != 1) fc = facetIndex[c];

     		const Array<double>& w = W_[fc];
     		for (int n=0; n<nNodes(); n++, count++)
     		{
    		  aPtr[count] += w_cell*w[n];
     		}
      }
    }
    addFlops(JVol.numCells() * (nNodes() + 2));
  }
  else
  {
    /* If the derivative order is nonzero, then we have to do a transformation. 
     * If we're also on a cell of dimension lower than maximal, we need to refer
     * to the facet index of the facet being integrated. */
    int nCells = JVol.numCells();
    double one = 1.0;
    int nTransRows = nRefDerivTest();

    createOneFormTransformationMatrix(JTrans, JVol);

    SUNDANCE_MSG3(transformVerb(),
      Tabs() << "transformation matrix=" << G(alpha()));
    int nNodes0 = nNodes();
      
    if (nFacetCases()==1)
    {
      /* 
       * on maximal cell
       */
      if (globalCurve().isCurveValid())
      {   /* ---------- ACI logic ----------- */
        TEUCHOS_TEST_FOR_EXCEPT(globalCurve().isCurveValid()); 
      }
      else           /* ---------- NO ACI logic----------- */
      {
        int N = 1;
        for (int c=0; c<JVol.numCells(); c++)
        {
          double* aPtr = &((*A)[c*nNodes0]);
          ::dgemm_("N", "N", &nNodes0, &N, &nTransRows, &(coeffs[c]), &(W_[0][0]),
            &nNodes0, &(G(alpha())[c*nTransRows]), &nTransRows, &one,
            aPtr, &nNodes0);
        }
      }
    }
    else
    {
      /* If we're on a lower-dimensional cell and have to transform, 
       * we've got to do each transformation using a different facet case */
      if (globalCurve().isCurveValid())
      {                 /* ---------- ACI logic ----------- */
        TEUCHOS_TEST_FOR_EXCEPT(globalCurve().isCurveValid()); 
      }
      else  /* ---------- NO ACI logic ----------- */
      {
        int N = 1;
        for (int c=0; c<JVol.numCells(); c++)
        {
          int fc = 0;
          if (nfc != 1) fc = facetIndex[c];
          double* aPtr = &((*A)[c*nNodes0]);
          ::dgemm_("N", "N", &nNodes0, &N, &nTransRows, &(coeffs[c]), &(W_[fc][0]),
            &nNodes0, &(G(alpha())[c*nTransRows]), &nTransRows, &one,
            aPtr, &nNodes0);
        }
      }
    }
      
    addFlops(2 * nNodes0 * nCells * nTransRows);
  }
}

void ReducedIntegral::transformTwoForm(const CellJacobianBatch& JTrans,
  const CellJacobianBatch& JVol,
  const Array<int>& facetIndex, 
  const RCP<Array<int> >& cellLIDs,
  const double* const coeffs,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(reduced2IntegrationTimer());
  TEUCHOS_TEST_FOR_EXCEPTION(order() != 2, std::logic_error,
    "ReducedIntegral::transformTwoForm() called for form "
    "of order " << order());
  
  Tabs tabs;  
  SUNDANCE_MSG1(transformVerb(), tabs << "doing two form by reduced integration");

  int nfc = nFacetCases();

  /* If the derivative orders are zero, the only transformation to be done 
   * is to multiply by the cell's Jacobian determinant */
  if (testDerivOrder() == 0 && unkDerivOrder() == 0)
  {
    if (globalCurve().isCurveValid())
    {     /* ----------- ACI logic ------------ */
      TEUCHOS_TEST_FOR_EXCEPT(globalCurve().isCurveValid()); 
    }
    else        /* ---------- NO ACI logic----------- */
    {
      double* aPtr = &((*A)[0]);
      int count = 0;
      for (int c=0; c<JVol.numCells(); c++)
      {
        int fc = 0;
        if (nFacetCases() != 1) fc = facetIndex[c];

        const Array<double>& w = W_[fc];
        double w_cell = coeffs[c] * fabs(JVol.detJ()[c]);
        for (int n=0; n<nNodes(); n++, count++)
        {
          aPtr[count] += w_cell*w[n];
        }
      }
    }
    addFlops(JVol.numCells() * (nNodes() + 1));
  }
  else
  {
    /* If the derivative order is nonzero, then we have to do a transformation. 
     * If we're also on a cell of dimension lower than maximal, we need to refer
     * to the facet index of the facet being integrated. */
    int nCells = JVol.numCells();
    double one = 1.0;
    int nTransRows = nRefDerivUnk()*nRefDerivTest();

    createTwoFormTransformationMatrix(JTrans, JVol);
      
    double* GPtr;
    if (testDerivOrder() == 0)
    {
      GPtr = &(G(beta())[0]);
      SUNDANCE_MSG2(transformVerb(),
        Tabs() << "transformation matrix=" << G(beta()));
    }
    else if (unkDerivOrder() == 0)
    {
      GPtr = &(G(alpha())[0]);
      SUNDANCE_MSG2(transformVerb(),
        Tabs() << "transformation matrix=" << G(alpha()));
    }
    else
    {
      GPtr = &(G(alpha(), beta())[0]);
      SUNDANCE_MSG2(transformVerb(),
        Tabs() << "transformation matrix=" 
        << G(alpha(),beta()));
    }
      
    int nNodes0 = nNodes();

    if (nFacetCases()==1)
    {
      /* 
       * on maximal cell
       */
      if (globalCurve().isCurveValid())
      {          /* ---------- ACI logic ----------- */
        TEUCHOS_TEST_FOR_EXCEPT(globalCurve().isCurveValid());         
      }
      else /* ---------- NO ACI ----------- */
      {
        int N = 1;
        for (int c=0; c<JVol.numCells(); c++)
        {
          double* aPtr = &((*A)[c*nNodes0]);
          double* gPtr = &(GPtr[c*nTransRows]);
          SUNDANCE_MSG2(integrationVerb(),
            tabs << "transforming c=" << c << ", W=" << W_[0]);
          
          ::dgemm_("N", "N", &nNodes0, &N, &nTransRows, &(coeffs[c]), &(W_[0][0]),
            &nNodes0, gPtr, &nTransRows, &one,
            aPtr, &nNodes0);
        }
      }
    }
    else
    {
      /* If we're on a lower-dimensional cell and have to transform, 
       * we've got to do each transformation using a different facet case */
      if (globalCurve().isCurveValid())
      {   /* ---------- ACI logic ----------- */
        TEUCHOS_TEST_FOR_EXCEPT(globalCurve().isCurveValid());         
      }
      else         /* ---------- NO ACI ----------- */
      {
        int N = 1;
        for (int c=0; c<JVol.numCells(); c++)
        {
          int fc = 0;
          if (nfc != 1) fc = facetIndex[c];
          double* aPtr = &((*A)[c*nNodes0]);
          double* gPtr = &(GPtr[c*nTransRows]);
          SUNDANCE_MSG2(integrationVerb(),
            tabs << "c=" << c << ", facet case=" << fc
            << " W=" << W_[fc]);

          ::dgemm_("N", "N", &nNodes0, &N, &nTransRows, &(coeffs[c]), &(W_[fc][0]),
            &nNodes0, gPtr, &nTransRows, &one,
            aPtr, &nNodes0);
        }
      }
    }// from else of (nFacetCases()==1)
      
    addFlops(2 * nNodes0 * nCells * nTransRows);
  }
}

