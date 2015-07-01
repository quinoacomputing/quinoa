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

#include "SundanceRefIntegral.hpp"
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

static Time& ref0IntegrationTimer() 
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("ref 0-form integration"); 
  return *rtn;
}

static Time& ref1IntegrationTimer() 
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("ref 1-form integration"); 
  return *rtn;
}


static Time& ref2IntegrationTimer() 
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("ref 2-form integration"); 
  return *rtn;
}


RefIntegral::RefIntegral(int spatialDim,
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
    tab0 << "************* creating reference 0-form integrals ********");
  if (setupVerb()) describe(Out::os());
  
  /* we need to sum the quadrature weights 
     to compute the volume of the reference cell */
  QuadratureFamily quad = new GaussianQuadrature(2);
  /* If we have a valid curve (in case of Adaptive Cell Integration)
   * then we have to choose the quadrature which the user specified*/
  if (globalCurve.isCurveValid()){
	 quad = quad_in;
	 Tabs tab1;
	 SUNDANCE_MSG1(setupVerb(),tab1 << "ACI change quadrature to Quadrature of order: "<<quad.order());
  }
  quad_ = quad;

  Array<Point> quadPts;
  Array<double> quadWeights;

  W_.resize(1);
  W_[0].resize(1);

  quad.getPoints(cellType, quadPts, quadWeights);  

  quadWeights_ = quadWeights;

  for (int q=0; q<quadWeights.size(); q++) {
	  W_[0][0] += quadWeights[q];
  }
}

RefIntegral::RefIntegral(int spatialDim,
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
    tab0 << "************* creating reference 1-form integrals ********");
  if (setupVerb()) describe(Out::os());
  assertLinearForm();

  W_.resize(nFacetCases());
  W_ACI_F1_.resize(nFacetCases());

  /* Determine the quadrature order needed for exact integrations */
  QuadratureType qType = new GaussianQuadratureType();
  int reqOrder = qType.findValidOrder(cellType, 
    std::max(1, testBasis.order()));

  SUNDANCE_MSG2(setupVerb(),
    tab0 << "using quadrature order=" << reqOrder);
   
  /* Create a quadrature family of the required order */
  QuadratureFamily quad = qType.createQuadFamily(reqOrder);
  
  /* If we have a valid curve (in case of Adaptive Cell Integration)
   * then we have to choose the quadrature which the user specified*/
  if (globalCurve.isCurveValid()){
	 quad = quad_in;
	 Tabs tab1;
	 SUNDANCE_MSG1(setupVerb(),tab1 << "ACI change quadrature to Quadrature of order: "<<quad.order());
  }
  quad_ = quad;

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

    getQuad(quad, fc, quadPts_, quadWeights_);

    int nQuad = quadPts_.size();

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
      testBasis.refEval(evalCellType(), quadPts_, deriv,
        testBasisVals[r], setupVerb());
    }

    /* do the quadrature */
    SUNDANCE_MSG2(setupVerb(), tab1 << "doing quadrature");
    int vecComp = 0;
    W_ACI_F1_[fc].resize(nQuad);
    for (int q=0; q<nQuad; q++)
    {
      W_ACI_F1_[fc][q].resize(nRefDerivTest());
      for (int t=0; t<nRefDerivTest(); t++)
      {
    	W_ACI_F1_[fc][q][t].resize(nNodesTest());
        for (int nt=0; nt<nNodesTest(); nt++)
        {
          value(fc, t, nt) 
            += chop(quadWeights_[q] * testBasisVals[t][vecComp][q][nt]) ;
          W_ACI_F1_[fc][q][t][nt] = chop(testBasisVals[t][vecComp][q][nt]);
        }
      }
    }    

    for (int i=0; i<W_[fc].size(); i++) W_[fc][i] = chop(W_[fc][i]);

    addFlops(3*nQuad*nRefDerivTest()*nNodesTest() + W_[fc].size());
  }

  /* print the result */
  SUNDANCE_MSG4(setupVerb(), tab0 << "--------------------------------------");
  SUNDANCE_MSG4(setupVerb(), tab0 << "reference linear form integral results");
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

  SUNDANCE_MSG1(setupVerb(), tab0 << "done reference linear form ctor");
}




RefIntegral::RefIntegral(int spatialDim,
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
    tab0 << "************* creating reference 2-form integrals ********");
  if (setupVerb()) describe(Out::os());

  assertBilinearForm();

  W_.resize(nFacetCases());
  W_ACI_F2_.resize(nFacetCases());

  QuadratureType qType = new GaussianQuadratureType();
  int reqOrder = qType.findValidOrder(cellType,
      std::max(1, unkBasis.order() + testBasis.order()));

  SUNDANCE_MSG2(setupVerb(),
      tab0 << "using quadrature order=" << reqOrder);
  QuadratureFamily quad = qType.createQuadFamily(reqOrder);

  /* If we have a valid curve (in case of Adaptive Cell Integration)
   * then we have to choose the quadrature which the user specified*/
  if (globalCurve.isCurveValid()){
	 quad = quad_in;
	 Tabs tab1;
	 SUNDANCE_MSG1(setupVerb(),tab1 << "ACI change quadrature to Quadrature of order: "<<quad.order());
  }
  quad_ = quad;

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
        
    getQuad(quad, fc, quadPts_, quadWeights_);
    int nQuad = quadPts_.size();

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
      testBasis.refEval(evalCellType(), quadPts_, deriv,
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
    	quadPts_, deriv, unkBasisVals[r], setupVerb());
    }

    SUNDANCE_MSG2(setupVerb(), tab1 << "doing quadrature...");
    int vecComp = 0;
    W_ACI_F2_[fc].resize(nQuad);
    for (int q=0; q<nQuad; q++)
    {
      W_ACI_F2_[fc][q].resize(nRefDerivTest());
      for (int t=0; t<nRefDerivTest(); t++)
      {
        W_ACI_F2_[fc][q][t].resize(nNodesTest());
        for (int nt=0; nt<nNodesTest(); nt++)
        {
          W_ACI_F2_[fc][q][t][nt].resize(nRefDerivUnk());
          for (int u=0; u<nRefDerivUnk(); u++)
          {
            W_ACI_F2_[fc][q][t][nt][u].resize(nNodesUnk());
            for (int nu=0; nu<nNodesUnk(); nu++)
            {
              value(fc, t, nt, u, nu) 
                += chop(quadWeights_[q] * testBasisVals[t][vecComp][q][nt]
                  * unkBasisVals[u][vecComp][q][nu]);
              W_ACI_F2_[fc][q][t][nt][u][nu] = chop( testBasisVals[t][vecComp][q][nt]
                                               * unkBasisVals[u][vecComp][q][nu] );
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
    << "reference bilinear form integral results");
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

  SUNDANCE_MSG1(setupVerb(), tab0 << "done reference bilinear form ctor");
}




void RefIntegral::transformZeroForm(const CellJacobianBatch& JVol,
  const Array<int>& isLocalFlag,  
  const RCP<Array<int> >& cellLIDs,
  const double& coeff,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(ref0IntegrationTimer());

  TEUCHOS_TEST_FOR_EXCEPTION(order() != 0, std::logic_error,
    "RefIntegral::transformZeroForm() called "
    "for form of order " << order());

  Tabs tabs;  
  SUNDANCE_MSG1(integrationVerb(), tabs << "doing zero form by reference");

  double& a = (*A)[0];
  int flops = 0;
  const Array<int>* cellLID = cellLIDs.get();

  /* if we don't need to check whether elements are local, we
   * can streamline the loop. This will be the case when
   * we are evaluating a functional but not its gradient */
  double w = coeff * W_[0][0];
  if ((int) isLocalFlag.size()==0)
  {
     	if (globalCurve().isCurveValid())
     	{     /* ---------- ACI logic ----------- */

     		Array<double> quadWeightsTmp = quadWeights_;
     		Array<Point> quadPointsTmp = quadPts_;
     		bool isCutByCurve;

     		for (int c=0; c<JVol.numCells(); c++)
     		{
     			int fc = 0;
   				/* call the special integration routine */
   				quadWeightsTmp = quadWeights_;
   				quadPointsTmp = quadPts_;
   				quad_.getAdaptedWeights(cellType(), dim(), (*cellLID)[c], fc ,mesh(),
   						globalCurve(), quadPointsTmp, quadWeightsTmp, isCutByCurve);
   				/* if we have special weights then do the same as before */
   				if (isCutByCurve){
   					double sumweights = 0;
   					for (int j=0; j < quadWeightsTmp.size(); j++) sumweights += chop(quadWeightsTmp[j]);
   					flops+=3+quadWeightsTmp.size();  //Todo: the curve stuff not counted
   					a += coeff * sumweights * fabs(JVol.detJ()[c]);
   				} else {
   					flops+=2;  //Todo: the curve stuff not counted
   					a += w * fabs(JVol.detJ()[c]);
   				}
     		}
     	}
     	else /* -------- NO ACI logic ------- */
     	{
     		for (int c=0; c<JVol.numCells(); c++)
     		{
 				flops+=2;
 				a += w * fabs(JVol.detJ()[c]);
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

      int fc = 0;
      if (globalCurve().isCurveValid())
      {   /* ---------- ACI logic ----------- */
    		Array<double> quadWeightsTmp = quadWeights_;
     		Array<Point> quadPointsTmp = quadPts_;
     		bool isCutByCurve;

    		for (int c=0; c<JVol.numCells(); c++)
    		{
    		  if (isLocalFlag[c])
    		  {

    			/* call the special integration routine */
    			quadWeightsTmp = quadWeights_;
    			quadPointsTmp = quadPts_;
    			quad_.getAdaptedWeights(cellType(), dim(), (*cellLID)[c], fc , mesh(),
    					globalCurve(), quadPointsTmp, quadWeightsTmp, isCutByCurve);
    			/* if we do not have special weights then do the same as before */
    			if (isCutByCurve){
    				double sumweights = 0;
    				for (int j=0; j < quadWeightsTmp.size(); j++) sumweights += chop(quadWeightsTmp[j]);
    				flops+=3+quadWeightsTmp.size();  //Todo: the curve stuff not counted
    				a += coeff * sumweights * fabs(JVol.detJ()[c]);
    			} else {
    				flops+=2;  //Todo: the curve stuff not counted
    				a += w * fabs(JVol.detJ()[c]);
    			}
    		  }
    		}
    	}
        else         /* ---------- NO ACI logic ----------- */
    	{
    		for (int c=0; c<JVol.numCells(); c++)
    		{
      		  if (isLocalFlag[c])
      		  {
    			flops+=2;
    			a += w * fabs(JVol.detJ()[c]);
      		  }
    		}
    	}
  }
  addFlops(flops);
}

void RefIntegral::transformOneForm(const CellJacobianBatch& JTrans,  
  const CellJacobianBatch& JVol,
  const Array<int>& facetIndex,
  const RCP<Array<int> >& cellLIDs,
  const double& coeff,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(ref1IntegrationTimer());
  TEUCHOS_TEST_FOR_EXCEPTION(order() != 1, std::logic_error,
    "RefIntegral::transformOneForm() called for form "
    "of order " << order());
  Tabs tabs;  
  SUNDANCE_MSG1(integrationVerb(),
    tabs << "doing one form by reference");

	int nfc = nFacetCases();

  /* If the derivative order is zero, the only transformation to be done 
   * is to multiply by the cell's Jacobian determinant */
  if (testDerivOrder() == 0)
  {
    double* aPtr = &((*A)[0]);
    int count = 0;
    if (globalCurve().isCurveValid())
    {    /* ---------- ACI logic ----------- */

    	 Array<double> quadWeightsTmp = quadWeights_;
    	 Array<Point> quadPointsTmp = quadPts_;
    	 bool isCutByCurve = false;

    	 for (int c=0; c<JVol.numCells(); c++)
    	 {
    	   double detJ = coeff * fabs(JVol.detJ()[c]);
    	   int fc = 0;
    	   if (nfc != 1) fc = facetIndex[c];

    	   quadWeightsTmp = quadWeights_;
    	   quadPointsTmp = quadPts_;
    	   /* call the special integration routine */
    	   quad_.getAdaptedWeights(cellType(), dim(), (*cellLIDs)[c] , fc ,
    			   mesh() , globalCurve() , quadPointsTmp , quadWeightsTmp , isCutByCurve );
    	   if (isCutByCurve)
    	   {
    		   Array<double> w;
    		   w.resize(nNodesTest()); //recalculate the special weights
    		   for (int nt = 0 ; nt < nNodesTest() ; nt++){
    			   w[nt] = 0.0;
    			   for (int q=0 ; q < quadWeightsTmp.size() ; q++)
    				   w[nt] += chop(quadWeightsTmp[q] * W_ACI_F1_[fc][q][0][nt]);
    		   }
    		   // do the integration here
    		   for (int n=0; n<nNodes(); n++, count++)
    		   {
    			   aPtr[count] += detJ*w[n];
    		   }
    	    }
    	    else
    	    {
         		const Array<double>& w = W_[fc];
    	    	for (int n=0; n<nNodes(); n++, count++)
    	    	{
    	    	  aPtr[count] += detJ*w[n];
    	    	}
    	    }
    	 }
    }
    else         /* ---------- NO ACI logic ----------- */
    {
     	 for (int c=0; c<JVol.numCells(); c++)
     	 {
     	    double detJ = coeff * fabs(JVol.detJ()[c]);
     	    int fc = 0;
     	    if (nfc != 1) fc = facetIndex[c];

     		const Array<double>& w = W_[fc];
     		for (int n=0; n<nNodes(); n++, count++)
     		{
    		  aPtr[count] += detJ*w[n];
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
    int nTransRows = nRefDerivTest();

    createOneFormTransformationMatrix(JTrans, JVol);

    SUNDANCE_MSG3(transformVerb(),
      Tabs() << "transformation matrix=" << G(alpha()));
    int nNodes0 = nNodes();
      
    if (nFacetCases()==1)
    {
      /* if we're on a maximal cell, we can do transformations 
       * for all cells in a single batch. 
       */
      if (globalCurve().isCurveValid())
      {   /* ---------- ACI logic ----------- */

     	 Array<double> quadWeightsTmp = quadWeights_;
     	 Array<Point> quadPointsTmp = quadPts_;
     	 bool isCutByCurve = false;

    	 double* aPtr_tmp = &((*A)[0]);
    	 int count = 0;
    	 const int oneI = 1;
    	 for (int c=0; c<JVol.numCells(); c++)
    	 {
    	     int fc = 0;
    	     if (nfc != 1) fc = facetIndex[c];
    		 /* call the special integration routine */
	         quadWeightsTmp = quadWeights_;
	         quadPointsTmp = quadPts_;
    		 quad_.getAdaptedWeights(cellType(), dim(), (*cellLIDs)[c], fc ,
           mesh() , globalCurve() , quadPointsTmp , quadWeightsTmp , isCutByCurve);
    		 if (isCutByCurve){
    			 Array<double> w;
    			 w.resize(nRefDerivTest()*nNodes()); //recalculate the special weights
				 for (int i = 0 ; i < nRefDerivTest()*nNodes() ; w[i] = 0.0 , i++);
    		     for (int t=0; t<nRefDerivTest(); t++)
    			   for (int nt = 0 ; nt < nNodesTest() ; nt++){
    				 for (int q=0 ; q < quadWeightsTmp.size() ; q++)
    					 //Index formula: nNodesTest()*testDerivDir + testNode
    					 w[nNodesTest()*t + nt] += chop(quadWeightsTmp[q] * W_ACI_F1_[0][q][t][nt]);
    			   }
   		         ::dgemm_("N", "N", &nNodes0, &oneI , &nTransRows, &coeff, &(w[0]),
   		            &nNodes0, &(G(alpha())[0]), &nTransRows, &one,
   		            &(aPtr_tmp[count]), &nNodes0);
    		 }else{
    		     ::dgemm_("N", "N", &nNodes0, &oneI , &nTransRows, &coeff, &(W_[0][0]),
    		        &nNodes0, &(G(alpha())[0]), &nTransRows, &one,
    		        &(aPtr_tmp[count]), &nNodes0);
    		 }
             count += nNodes();
    	 } // end from the for loop over the cells
      }
      else           /* ---------- NO ACI logic----------- */
      {
      ::dgemm_("N", "N", &nNodes0, &nCells, &nTransRows, &coeff, &(W_[0][0]),
        &nNodes0, &(G(alpha())[0]), &nTransRows, &one, 
        &((*A)[0]), &nNodes0);
      }
    }
    else
    {
      /* If we're on a lower-dimensional cell and have to transform, 
       * we've got to do each transformation using a different facet case */
        if (globalCurve().isCurveValid())
        {                 /* ---------- ACI logic ----------- */

        	Array<double> quadWeightsTmp = quadWeights_;
        	Array<Point> quadPointsTmp = quadPts_;
        	bool isCutByCurve = false;

            int N = 1;
            for (int c=0; c<JVol.numCells(); c++)
            {
                int fc = 0;
                if (nfc != 1) fc = facetIndex[c];
                double* aPtr = &((*A)[c*nNodes0]);
            	/* call the special integration routine */
            	quadWeightsTmp = quadWeights_;
            	quadPointsTmp = quadPts_;
            	quad_.getAdaptedWeights(cellType(), dim(), (*cellLIDs)[c], fc ,
            			mesh() , globalCurve(), quadPointsTmp, quadWeightsTmp, isCutByCurve);
            	if (isCutByCurve){
            		Array<double> w;
            		w.resize(nRefDerivTest()*nNodes()); //recalculate the special weights
            		for (int i = 0 ; i < nRefDerivTest()*nNodes() ; w[i] = 0.0 , i++);
            		for (int t=0; t<nRefDerivTest(); t++)
            			for (int nt = 0 ; nt < nNodesTest() ; nt++){
            				for (int q=0 ; q < quadWeightsTmp.size() ; q++)
            					//Index formula: nNodesTest()*testDerivDir + testNode
            					w[nNodesTest()*t + nt] += chop(quadWeightsTmp[q] * W_ACI_F1_[fc][q][t][nt]);
            			}
            		::dgemm_("N", "N", &nNodes0, &N , &nTransRows, &coeff, &(w[0]),
            				&nNodes0, &(G(alpha())[0]), &nTransRows, &one,
            				aPtr, &nNodes0);
            	}else{
            		::dgemm_("N", "N", &nNodes0, &N, &nTransRows, &coeff, &(W_[fc][0]),
            				&nNodes0, &(G(alpha())[c*nTransRows]), &nTransRows, &one,
            				aPtr, &nNodes0);
            	}
            }
        }
        else  /* ---------- NO ACI logic ----------- */
        {
            int N = 1;
            for (int c=0; c<JVol.numCells(); c++)
            {
                int fc = 0;
                if (nfc != 1) fc = facetIndex[c];
                double* aPtr = &((*A)[c*nNodes0]);
                ::dgemm_("N", "N", &nNodes0, &N, &nTransRows, &coeff, &(W_[fc][0]),
                		&nNodes0, &(G(alpha())[c*nTransRows]), &nTransRows, &one,
                		aPtr, &nNodes0);
            }
        }
    }
      
    addFlops(2 * nNodes0 * nCells * nTransRows);
  }
}

void RefIntegral::transformTwoForm(const CellJacobianBatch& JTrans,
  const CellJacobianBatch& JVol,
  const Array<int>& facetIndex, 
  const RCP<Array<int> >& cellLIDs,
  const double& coeff,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(ref2IntegrationTimer());
  TEUCHOS_TEST_FOR_EXCEPTION(order() != 2, std::logic_error,
    "RefIntegral::transformTwoForm() called for form "
    "of order " << order());
  
  Tabs tabs;  
  SUNDANCE_MSG1(transformVerb(), tabs << "doing two form by reference");

  int nfc = nFacetCases();

	  SUNDANCE_MSG1(transformVerb(), tabs << "doing two form by reference ... ");
  /* If the derivative orders are zero, the only transformation to be done 
   * is to multiply by the cell's Jacobian determinant */
  if (testDerivOrder() == 0 && unkDerivOrder() == 0)
  {
      if (globalCurve().isCurveValid())
      {     /* ----------- ACI logic ------------ */

    	   Array<double> quadWeightsTmp = quadWeights_;
    	   Array<Point> quadPointsTmp = quadPts_;
    	   bool isCutByCurve = false;

    	   double* aPtr = &((*A)[0]);
    	   int count = 0;
    	   for (int c=0; c<JVol.numCells(); c++)
    	   {
    	     int fc = 0;
    	     if (nFacetCases() != 1) fc = facetIndex[c];

    	     /* ---------- ACI ----------- */
    	     /* call the special integration routine */
    	     quadWeightsTmp = quadWeights_;
    	     quadPointsTmp = quadPts_;
    	     quad_.getAdaptedWeights(cellType(), dim(), (*cellLIDs)[c] , fc ,
    	    		 mesh(), globalCurve(), quadPointsTmp, quadWeightsTmp, isCutByCurve);
    	     if (isCutByCurve){
    	    	 Array<double> w;
    	    	 int ci = 0;
    	    	 w.resize(nNodesTest()*nNodesUnk()); //recalculate the special weights
    	    	 for (int nt = 0 ; nt < nNodesTest() ; nt++)
    	    		 for(int nu=0 ; nu < nNodesUnk() ; nu++ , ci++){
    	    			 w[ci] = 0.0;
    	    			 for (int q=0 ; q < quadWeightsTmp.size() ; q++)
    	    				 w[ci] += chop( quadWeightsTmp[q] * W_ACI_F2_[fc][q][0][nt][0][nu] );
    	    		 }
    	    	 // do the integration here
    	    	 double detJ = coeff * fabs(JVol.detJ()[c]);
    	    	 for (int n=0; n<nNodes(); n++, count++)
    	    	 {
    	    		 aPtr[count] += detJ*w[n];
    	    	 }
    	     }
    	     else
    	     {
    	    	  const Array<double>& w = W_[fc];
    	    	  double detJ = coeff * fabs(JVol.detJ()[c]);
    	    	  for (int n=0; n<nNodes(); n++, count++)
    	    	  {
    	    		  aPtr[count] += detJ*w[n];
    	    	  }
    	     }
    	   }
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
    		  double detJ = coeff * fabs(JVol.detJ()[c]);
    		  for (int n=0; n<nNodes(); n++, count++)
    		  {
    			  aPtr[count] += detJ*w[n];
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
      /* if we're on a maximal cell, we can do transformations 
       * for all cells in a single batch. 
       */
      if (globalCurve().isCurveValid())
      {          /* ---------- ACI logic ----------- */

    	 Array<double> quadWeightsTmp = quadWeights_;
    	 Array<Point> quadPointsTmp = quadPts_;
    	 bool isCutByCurve = false;

       	 for (int c=0; c<JVol.numCells(); c++)
       	 {
             int fc = 0;
             if (nfc != 1) fc = facetIndex[c];

             double* aPtr = &((*A)[c*nNodes0]);
             double* gPtr = &(GPtr[c*nTransRows]);
             int oneI = 1;
       		 /* call the special integration routine */
         	//SUNDANCE_MSG1(transformVerb(), tabs << "before quad_.getAdaptedWeights");
         	 quadWeightsTmp = quadWeights_;
             quadPointsTmp = quadPts_;
       		 quad_.getAdaptedWeights(cellType(), dim(), (*cellLIDs)[c], fc ,
             mesh(),globalCurve(), quadPointsTmp, quadWeightsTmp, isCutByCurve);
         	//SUNDANCE_MSG1(transformVerb(), tabs << "after quad_.getAdaptedWeights");
       		 if (isCutByCurve){
       			 Array<double> w;
       			 w.resize(nNodesUnk()*nNodesTest()*nRefDerivUnk()*nRefDerivTest());
       			 for ( int i = 0 ; i < w.size() ; i++) w[i] = 0.0;
       			 //recalculate the special weights
       		     for (int t=0; t<nRefDerivTest(); t++){
       		        for (int nt=0; nt<nNodesTest(); nt++)
       		          for (int u=0; u<nRefDerivUnk(); u++)
       		            for (int nu=0; nu<nNodesUnk(); nu++)
       		            	for (int q=0 ; q < quadWeightsTmp.size() ; q++)
       		                // unkNode + nNodesUnk()*testNode  + nNodes()*(unkDerivDir + nRefDerivUnk()*testDerivDir)
       		                    w[nu + nNodesUnk()*nt  + nNodes()*(u + nRefDerivUnk()*t)] +=
       		                    		chop(quadWeightsTmp[q]*W_ACI_F2_[0][q][t][nt][u][nu]);
       		     }
      		      ::dgemm_("N", "N", &nNodes0, &oneI , &nTransRows, &coeff, &(w[0]),
      		        &nNodes0, &(gPtr[0]), &nTransRows, &one,
      		        &(aPtr[0]), &nNodes0);
       		  }else{
       		     ::dgemm_("N", "N", &nNodes0, &oneI , &nTransRows, &coeff, &(W_[0][0]),
       		        &nNodes0, &(gPtr[0]), &nTransRows, &one,
       		        &(aPtr[0]), &nNodes0);
       		  }
       	 } // end from the for loop over the cells
      }
      else /* ---------- NO ACI ----------- */
      {
        	 ::dgemm_("N", "N", &nNodes0, &nCells, &nTransRows, &coeff, &(W_[0][0]),
        		 &nNodes0, GPtr, &nTransRows, &one,
        		 &((*A)[0]), &nNodes0);
      }
    }
    else
    {
      /* If we're on a lower-dimensional cell and have to transform, 
       * we've got to do each transformation using a different facet case */
        if (globalCurve().isCurveValid())
        {   /* ---------- ACI logic ----------- */
            int oneI = 1;
            Array<double> quadWeightsTmp = quadWeights_;
            Array<Point> quadPointsTmp = quadPts_;
            bool isCutByCurve = false;

            for (int c=0; c<JVol.numCells(); c++)
            {
              int fc = 0;
              if (nfc != 1) fc = facetIndex[c];
              double* aPtr = &((*A)[c*nNodes0]);
              double* gPtr = &(GPtr[c*nTransRows]);
              SUNDANCE_MSG2(integrationVerb(),
                tabs << "c=" << c << ", facet case=" << fc
                << " W=" << W_[fc]);

              /* call the special integration routine */
              quadWeightsTmp = quadWeights_;
              quadPointsTmp = quadPts_;
              quad_.getAdaptedWeights(cellType(), dim(), (*cellLIDs)[c], fc ,
            		  mesh(), globalCurve(), quadPointsTmp, quadWeightsTmp, isCutByCurve);
              if (isCutByCurve){
            	  Array<double> w;
            	  w.resize(nNodesUnk()*nNodesTest()*nRefDerivUnk()*nRefDerivTest());
            	  for ( int i = 0 ; i < w.size() ; i++) w[i] = 0.0;
            	  //recalculate the special weights
            	  for (int t=0; t<nRefDerivTest(); t++){
            		  for (int nt=0; nt<nNodesTest(); nt++)
            			  for (int u=0; u<nRefDerivUnk(); u++)
            				  for (int nu=0; nu<nNodesUnk(); nu++)
            					  for (int q=0 ; q < quadWeightsTmp.size() ; q++)
            						  // unkNode + nNodesUnk()*testNode  + nNodes()*(unkDerivDir + nRefDerivUnk()*testDerivDir)
            						  w[nu + nNodesUnk()*nt  + nNodes()*(u + nRefDerivUnk()*t)] +=
            								  chop( quadWeightsTmp[q]*W_ACI_F2_[fc][q][t][nt][u][nu] );
            	  }
            	  ::dgemm_("N", "N", &nNodes0, &oneI , &nTransRows, &coeff, &(w[0]),
            			  &nNodes0, &(gPtr[0]), &nTransRows, &one,
            			  &(aPtr[0]), &nNodes0);
				  }else{
					  ::dgemm_("N", "N", &nNodes0, &oneI , &nTransRows, &coeff, &(W_[fc][0]),
							  &nNodes0, &(gPtr[0]), &nTransRows, &one,
							  &(aPtr[0]), &nNodes0);
				  }
            }
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

              ::dgemm_("N", "N", &nNodes0, &N, &nTransRows, &coeff, &(W_[fc][0]),
            		  &nNodes0, gPtr, &nTransRows, &one,
            		  aPtr, &nNodes0);
            }
        }
    }// from else of (nFacetCases()==1)
      
    addFlops(2 * nNodes0 * nCells * nTransRows);
  }
}
