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

#include "SundanceQuadratureIntegral.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Teuchos;

using std::endl;
using std::setw;
using std::setprecision;

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

static Time& quadrature0Timer() 
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("0-form quadrature"); 
  return *rtn;
}

static Time& quadrature1Timer() 
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("1-form quadrature"); 
  return *rtn;
}

static Time& quadrature2Timer() 
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("2-form quadrature"); 
  return *rtn;
}


QuadratureIntegral::QuadratureIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim, 
  const CellType& cellType,
  const QuadratureFamily& quad,
  bool isInternalBdry,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh,
  int verb)
  : QuadratureIntegralBase(spatialDim, maxCellType, dim, cellType, quad, 
    isInternalBdry, globalCurve, mesh, verb),
    W_(),
    useSumFirstMethod_(true)
{
  Tabs tab0(0);
  
  SUNDANCE_MSG1(setupVerb(), tab0 << "QuadratureIntegral ctor for 0-form");
  if (setupVerb()) describe(Out::os());
  
  SUNDANCE_MSG1(setupVerb(), tab0 << "quadrature family=" << quad);  


  W_.resize(nFacetCases());
  quadPts_.resize(nFacetCases());
  quadWeights_.resize(nFacetCases());
  quad_ = quad;

  for (int fc=0; fc<nFacetCases(); fc++)
  {
    Tabs tab2;
    SUNDANCE_MSG1(setupVerb(), tab2 << "facet case=" << fc
      << " of " << nFacetCases());        
    Tabs tab3;
    /* create the quad points and weights */
    Array<double> quadWeights;
    Array<Point> quadPts;
    getQuad(quad, fc, quadPts, quadWeights);
    nQuad_ = quadPts.size();
    quadPts_[fc] = quadPts;
    quadWeights_[fc] = quadWeights;
      
    W_[fc].resize(nQuad());
      
    for (int q=0; q<nQuad(); q++)
    {
      W_[fc][q] = quadWeights[q];
    }
  }    
}


QuadratureIntegral::QuadratureIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim, 
  const CellType& cellType,
  const BasisFamily& testBasis,
  int alpha,
  int testDerivOrder,
  const QuadratureFamily& quad,
  bool isInternalBdry,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh,
  int verb)
  : QuadratureIntegralBase(spatialDim, maxCellType, dim, cellType, 
    testBasis, alpha, testDerivOrder, quad , isInternalBdry, globalCurve , mesh, verb),
    W_(),
    useSumFirstMethod_(true)
{
  Tabs tab0;
  
  SUNDANCE_MSG1(setupVerb(), tab0 << "QuadratureIntegral ctor for 1-form");
  if (setupVerb()) describe(Out::os());
  assertLinearForm();

  
  SUNDANCE_MSG1(setupVerb(), tab0 << "quadrature family=" << quad);  
  
  quad_ = quad;
  W_.resize(nFacetCases());
  // store the quadrature points and weights
  quadPts_.resize(nFacetCases());
  quadWeights_.resize(nFacetCases());
  W_ACI_F1_.resize(nFacetCases());

  for (int fc=0; fc<nFacetCases(); fc++)
  {
    Tabs tab2;

    SUNDANCE_MSG1(setupVerb(), tab2 << "facet case=" << fc
      << " of " << nFacetCases());        

    /* create the quad points and weights */
    Array<double> quadWeights;
    Array<Point> quadPts;
    getQuad(quad, fc, quadPts, quadWeights);
    nQuad_ = quadPts.size();
    quadPts_[fc] = quadPts;
    quadWeights_[fc] = quadWeights;
      
    W_[fc].resize(nQuad() * nRefDerivTest() * nNodesTest());

    SUNDANCE_MSG1(setupVerb(), tab2 << "num nodes for test function " << nNodesTest());

    Array<Array<Array<Array<double> > > > testBasisVals(nRefDerivTest());

    for (int r=0; r<nRefDerivTest(); r++)
    {
      Tabs tab3;
      SUNDANCE_MSG1(setupVerb(), tab3 
        << "evaluating basis functions for ref deriv direction " << r);
      MultiIndex mi;
      testBasisVals[r].resize(testBasis.dim());
      if (testDerivOrder==1) mi[r] = 1;
      SpatialDerivSpecifier deriv(mi);
      testBasis.refEval(evalCellType(), quadPts, deriv, 
        testBasisVals[r], setupVerb());
    }

      

    int vecComp = 0;
    W_ACI_F1_[fc].resize(nQuad());
    for (int q=0; q<nQuad(); q++)
    {
      W_ACI_F1_[fc][q].resize(nRefDerivTest());
      for (int t=0; t<nRefDerivTest(); t++)
      {
        W_ACI_F1_[fc][q][t].resize(nNodesTest());
        for (int nt=0; nt<nNodesTest(); nt++)
        {
          wValue(fc, q, t, nt) 
            = chop(quadWeights[q] * testBasisVals[t][vecComp][q][nt]) ;
          W_ACI_F1_[fc][q][t][nt] = chop(testBasisVals[t][vecComp][q][nt]);
        }
      }
    }

    addFlops(2*nQuad()*nRefDerivTest()*nNodesTest());
  }
}




QuadratureIntegral::QuadratureIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim,
  const CellType& cellType,
  const BasisFamily& testBasis,
  int alpha,
  int testDerivOrder,
  const BasisFamily& unkBasis,
  int beta,
  int unkDerivOrder,
  const QuadratureFamily& quad,
  bool isInternalBdry,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh,
  int verb)
  : QuadratureIntegralBase(spatialDim, maxCellType, dim, cellType, 
    testBasis, alpha, testDerivOrder, 
    unkBasis, beta, unkDerivOrder, quad , isInternalBdry, globalCurve , mesh , verb),
    W_(),
    useSumFirstMethod_(true)
{
  Tabs tab0;
  
  SUNDANCE_MSG1(setupVerb(), tab0 << "QuadratureIntegral ctor for 2-form");
  if (setupVerb()) describe(Out::os());
  assertBilinearForm();
  quad_ = quad;

  W_.resize(nFacetCases());
  W_ACI_F2_.resize(nFacetCases());
  // store the quadrature points and weights
  quadPts_.resize(nFacetCases());
  quadWeights_.resize(nFacetCases());

  for (int fc=0; fc<nFacetCases(); fc++)
  {
    /* get the quad pts and weights */
    Array<double> quadWeights;
    Array<Point> quadPts;
    getQuad(quad, fc, quadPts, quadWeights);
    // store the points for each facet case
    quadPts_[fc] = quadPts;
    quadWeights_[fc] = quadWeights;
    nQuad_ = quadPts.size();

    W_[fc].resize(nQuad() * nRefDerivTest() * nNodesTest()  
      * nRefDerivUnk() * nNodesUnk());


    /* compute the basis functions */
    Array<Array<Array<Array<double> > > > testBasisVals(nRefDerivTest());
    Array<Array<Array<Array<double> > > > unkBasisVals(nRefDerivUnk());


    for (int r=0; r<nRefDerivTest(); r++)
    {
      testBasisVals[r].resize(testBasis.dim());
      MultiIndex mi;
      if (testDerivOrder==1) mi[r] = 1;
      SpatialDerivSpecifier deriv(mi);
      testBasis.refEval(evalCellType(), quadPts, deriv, 
        testBasisVals[r], setupVerb());
    }

    for (int r=0; r<nRefDerivUnk(); r++)
    {
      unkBasisVals[r].resize(unkBasis.dim());
      MultiIndex mi;
      if (unkDerivOrder==1) mi[r] = 1;
      SpatialDerivSpecifier deriv(mi);
      unkBasis.refEval(evalCellType(), 
        quadPts, deriv, unkBasisVals[r], setupVerb());
    }


    int vecComp = 0;
    /* form the products of basis functions at each quad pt */
    W_ACI_F2_[fc].resize(nQuad());
    for (int q=0; q<nQuad(); q++)
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
              wValue(fc, q, t, nt, u, nu)
                = chop(quadWeights[q] * testBasisVals[t][vecComp][q][nt] 
                  * unkBasisVals[u][vecComp][q][nu]);
              W_ACI_F2_[fc][q][t][nt][u][nu] =
                chop(testBasisVals[t][vecComp][q][nt] * unkBasisVals[u][vecComp][q][nu]);
            }
          }
        }
      }
    }

    addFlops(3*nQuad()*nRefDerivTest()*nNodesTest()*nRefDerivUnk()*nNodesUnk()
      + W_[fc].size());
    for (int i=0; i<W_[fc].size(); i++) W_[fc][i] = chop(W_[fc][i]);
    
  }

}


void QuadratureIntegral::transformZeroForm(const CellJacobianBatch& JTrans,  
  const CellJacobianBatch& JVol,
  const Array<int>& isLocalFlag,
  const Array<int>& facetIndex,
  const RCP<Array<int> >& cellLIDs,
  const double* const coeff,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(quadrature0Timer());
  Tabs tabs;
  TEUCHOS_TEST_FOR_EXCEPTION(order() != 0, std::logic_error,
    "QuadratureIntegral::transformZeroForm() called "
    "for form of order " << order());

  TEUCHOS_TEST_FOR_EXCEPTION( (int) isLocalFlag.size() != 0 
    && (int) isLocalFlag.size() != JVol.numCells(),
    std::runtime_error,
    "mismatch between isLocalFlag.size()=" << isLocalFlag.size()
    << " and JVol.numCells()=" << JVol.numCells());

  bool checkLocalFlag = (int) isLocalFlag.size() != 0;

  const Array<int>* cellLID = cellLIDs.get();
  
  SUNDANCE_MSG1(integrationVerb(), tabs << "doing zero form by quadrature");
  double& a = (*A)[0];
  SUNDANCE_MSG5(integrationVerb(), tabs << "input A=");
  if (integrationVerb() >= 5) writeTable(Out::os(), tabs, *A, 6);
  double* coeffPtr = (double*) coeff;

  int flops = 0;
  if (nFacetCases()==1)
  {
      if (globalCurve().isCurveValid())
      {       /* ---------- ACI logic ----------- */

    	Array<double> quadWeightsTmp = quadWeights_[0];
    	Array<Point> quadPointsTmp = quadPts_[0];
    	bool isCutByCurve;

   	    const Array<double>& w = W_[0];
   	    for (int c=0; c<JVol.numCells(); c++)
   	    {
   	      if (checkLocalFlag && !isLocalFlag[c])
   	      {
   	        coeffPtr += nQuad();
   	        continue;
   	      }
   	      double detJ = fabs(JVol.detJ()[c]);
   	      int fc = 0;

          quadWeightsTmp = quadWeights_[fc];
          quadPointsTmp = quadPts_[fc];
          quad_.getAdaptedWeights(cellType(), dim(), (*cellLID)[c], fc ,mesh(),
               globalCurve(), quadPointsTmp, quadWeightsTmp, isCutByCurve);
      	   /* if we have special weights then do the same as before */
      	  if (isCutByCurve)
      	  {
               for (int q=0; q<nQuad(); q++, coeffPtr++)
               {
                 a += quadWeightsTmp[q]*(*coeffPtr)*detJ;
               }
               flops += 3*nQuad();
      	   }
      	   else
      	   {
               for (int q=0; q<nQuad(); q++, coeffPtr++)
               {
                 a += w[q]*(*coeffPtr)*detJ;
               }
               flops += 3*nQuad();
    	   }

         }
      }
      else        /* ---------- NO ACI logic ----------- */
      {
     	    const Array<double>& w = W_[0];
     	    for (int c=0; c<JVol.numCells(); c++)
     	    {
     	      if (checkLocalFlag && !isLocalFlag[c])
     	      {
     	        coeffPtr += nQuad();
     	        continue;
     	      }
     	      double detJ = fabs(JVol.detJ()[c]);

     	      for (int q=0; q<nQuad(); q++, coeffPtr++)
     	      {
     	    	  a += w[q]*(*coeffPtr)*detJ;
     	      }
     	      flops += 3*nQuad();
     	    }
      }
  }
  else
  {
      if (globalCurve().isCurveValid())
      {     /* ---------- ACI logic ----------- */
			Array<double> quadWeightsTmp = quadWeights_[0];
			Array<Point> quadPointsTmp = quadPts_[0];
			bool isCutByCurve;

    	    for (int c=0; c<JVol.numCells(); c++)
    	    {
    	      if (checkLocalFlag && !isLocalFlag[c])
    	      {
    	        coeffPtr += nQuad();
    	        continue;
    	      }
    	      double detJ = fabs(JVol.detJ()[c]);
    	      int fc = facetIndex[c];
    	      const Array<double>& w = W_[facetIndex[c]];


    	      quadWeightsTmp = quadWeights_[fc];
    	      quadPointsTmp = quadPts_[fc];
    	      quad_.getAdaptedWeights(cellType(), dim(), (*cellLID)[c], fc ,mesh(),
    	    		  globalCurve(), quadPointsTmp, quadWeightsTmp, isCutByCurve);
    	      /* if we have special weights then do the same as before */
    	      if (isCutByCurve){
    	    	  for (int q=0; q<nQuad(); q++, coeffPtr++)
    	    		  a += quadWeightsTmp[q]*(*coeffPtr)*detJ;
    	    	  flops += 3*nQuad();
    	      } // end cut by curve
    	      else{
    	    	  for (int q=0; q<nQuad(); q++, coeffPtr++)
    	    		  a += w[q]*(*coeffPtr)*detJ;
    	    	  flops += 3*nQuad();
    	      }
    	    }
      }
      else         /* ---------- NO ACI logic----------- */
      {
    	    for (int c=0; c<JVol.numCells(); c++)
    	    {
	      Tabs tab1;
	      SUNDANCE_MSG4(integrationVerb(), tab1 << "cell #" << c
			    << " of " << JVol.numCells());
    	      if (checkLocalFlag && !isLocalFlag[c])
    	      {
		Tabs tab2;
		SUNDANCE_MSG4(integrationVerb(), tab2 
			      << "skipped -- is not local cell");
    	        coeffPtr += nQuad();
    	        continue;
    	      }
    	      double detJ = fabs(JVol.detJ()[c]);
    	      const Array<double>& w = W_[facetIndex[c]];

    	      for (int q=0; q<nQuad(); q++, coeffPtr++)
    	      {
		Tabs tab3;
		SUNDANCE_MSG4(integrationVerb(), 
			      tab3 << "q=" << q << "\t w=" << w[q]
			      << "\t\t coeff=" << *coeffPtr << "\t\t |J|=" << detJ);
		a += w[q]*(*coeffPtr)*detJ;
    	      }
    	      flops += 3*nQuad();
    	    }
      }
  }
  SUNDANCE_MSG5(integrationVerb(), tabs << "output A = ");
  if (integrationVerb() >= 5) writeTable(Out::os(), tabs, *A, 6);

  addFlops(flops);
}


void QuadratureIntegral::transformOneForm(const CellJacobianBatch& JTrans,  
  const CellJacobianBatch& JVol,
  const Array<int>& facetIndex,
  const RCP<Array<int> >& cellLIDs,
  const double* const coeff,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(quadrature1Timer());
  Tabs tabs;
  TEUCHOS_TEST_FOR_EXCEPTION(order() != 1, std::logic_error,
    "QuadratureIntegral::transformOneForm() called for form "
    "of order " << order());
  SUNDANCE_MSG2(integrationVerb(), tabs << "doing one form by quadrature");
  int flops = 0;

  const Array<int>* cellLID = cellLIDs.get();

  /* If the derivative order is zero, the only thing to be done 
   * is to multiply by the cell's Jacobian determinant and sum over the
   * quad points */
  if (testDerivOrder() == 0)
  {
    double* aPtr = &((*A)[0]);
    SUNDANCE_MSG5(integrationVerb(), tabs << "input A = ");
    if (integrationVerb() >= 5) writeTable(Out::os(), tabs, *A, 6);
  
    double* coeffPtr = (double*) coeff;
    int offset = 0 ;

      if (globalCurve().isCurveValid())
      {         /* ---------- ACI logic ----------- */

    	    Array<double> quadWeightsTmp = quadWeights_[0];
    	    Array<Point> quadPointsTmp = quadPts_[0];
    	    bool isCutByCurve;

    	    for (int c=0; c<JVol.numCells(); c++, offset+=nNodes())
    	    {
    	      Tabs tab2;
    	      double detJ = fabs(JVol.detJ()[c]);
    	      int fc = 0;
    	      if (nFacetCases() != 1) fc = facetIndex[c];
    	      const Array<double>& w = W_[fc];
    	      SUNDANCE_MSG4(integrationVerb(), tab2 << "c=" << c << " detJ=" << detJ);

    	      quadWeightsTmp = quadWeights_[fc];
    	      quadPointsTmp = quadPts_[fc];
    	      /* call the special integration routine */
    	      quad_.getAdaptedWeights(cellType(), dim(), (*cellLID)[c] , fc ,
    	    		  mesh() , globalCurve() , quadPointsTmp , quadWeightsTmp , isCutByCurve );
    	      if (isCutByCurve){
    	    	  Array<double> wi;
    	    	  wi.resize(nQuad() * nNodes()); //recalculate the special weights
    	    	  for (int ii = 0 ; ii < wi.size() ; ii++ ) { wi[ii] = 0.0; }
    	    	  for (int n = 0 ; n < nNodes() ; n++)
    	    		  for (int q=0 ; q < quadWeightsTmp.size() ; q++)
    	    			  //Indexing: testNode + nNodesTest()*(testDerivDir + nRefDerivTest()*q)
    	    			  wi[n + nNodes()*q] +=
    	    					  chop(quadWeightsTmp[q] * W_ACI_F1_[fc][q][0][n]);
    	    	  // if it is cut by curve then use this vector
    	    	  for (int q=0; q<nQuad(); q++, coeffPtr++){
    	    		  double f = (*coeffPtr)*detJ;
    	    		  for (int n=0; n<nNodes(); n++){
    	    			  aPtr[offset+n] += f*wi[n + nNodes()*q];
    	    		  }
    	    	  }
    	      } // end isCutByCurve
    	      else
    	      {
    	    	  for (int q=0; q<nQuad(); q++, coeffPtr++){
    	    		  double f = (*coeffPtr)*detJ;
    	    		  for (int n=0; n<nNodes(); n++){
    	    			  aPtr[offset+n] += f*w[n + nNodes()*q];
    	    		  }
    	    	  }
    	      }
    	   } // from cell loop
      }
      else    /* ---------- NO ACI logic ----------- */
      {
  	    for (int c=0; c<JVol.numCells(); c++, offset+=nNodes())
  	    {
  	      Tabs tab2;
  	      double detJ = fabs(JVol.detJ()[c]);
  	      int fc = 0;
  	      if (nFacetCases() != 1) fc = facetIndex[c];
  	      const Array<double>& w = W_[fc];
  	      SUNDANCE_MSG4(integrationVerb(), tab2 << "c=" << c << " detJ=" << detJ);

    	  for (int q=0; q<nQuad(); q++, coeffPtr++)
    	  {
    		  Tabs tab3;
    		  double f = (*coeffPtr)*detJ;
    		  SUNDANCE_MSG4(integrationVerb(), tab3 << "q=" << q << " coeff=" <<
            *coeffPtr << " coeff*detJ=" << f);
    		  for (int n=0; n<nNodes(); n++)
    		  {
    			  Tabs tab4;
    			  SUNDANCE_MSG4(integrationVerb(), tab4 << "n=" << n << " w=" <<
              w[n + nNodes()*q]);
    			  aPtr[offset+n] += f*w[n + nNodes()*q];
    		  }
    	  }
  	    } // from for loop over cells
      }
      if (integrationVerb() >= 4)
      {
	    Tabs tab2;
        Out::os() << tab2 << "integration results on cell:" << std::endl;
        Out::os() << tab2 << setw(10) << "n" << setw(12) << "I_n" << std::endl;
        for (int n=0; n<nNodes(); n++) 
        {
          Out::os() << tab2 << setw(10) << n 
                    << setw(12) << setprecision(6) << aPtr[offset+n] << std::endl;
        }
      }
      
    SUNDANCE_MSG5(integrationVerb(), tabs << "output A = ");
    if (integrationVerb() >= 5) writeTable(Out::os(), tabs, *A, 6);
    addFlops( JVol.numCells() * (1 + nQuad() * (1 + 2*nNodes())) );
  }
  else
  {
    /* If the derivative order is nonzero, then we have to do a transformation. 
     * If we're also on a cell of dimension lower than maximal, we need to refer
     * to the facet index of the facet being integrated. */

    createOneFormTransformationMatrix(JTrans, JVol);

    SUNDANCE_MSG4(transformVerb(), 
      Tabs() << "transformation matrix=" << G(alpha()));

    double* GPtr = &(G(alpha())[0]);      

    if (useSumFirstMethod())
    {
      transformSummingFirst(JVol.numCells(), facetIndex, cellLIDs, GPtr, coeff, A);
    }
    else
    {
      transformSummingLast(JVol.numCells(), facetIndex, cellLIDs, GPtr, coeff, A);
    }
  }
  addFlops(flops);
}


void QuadratureIntegral::transformTwoForm(const CellJacobianBatch& JTrans,
  const CellJacobianBatch& JVol,
  const Array<int>& facetIndex,
  const RCP<Array<int> >& cellLIDs,
  const double* const coeff,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(quadrature2Timer());
  Tabs tabs;
  TEUCHOS_TEST_FOR_EXCEPTION(order() != 2, std::logic_error,
    "QuadratureIntegral::transformTwoForm() called for form "
    "of order " << order());
  SUNDANCE_MSG2(integrationVerb(), tabs << "doing two form by quadrature");

  const Array<int>* cellLID = cellLIDs.get();

  /* If the derivative orders are zero, the only thing to be done 
   * is to multiply by the cell's Jacobian determinant and sum over the
   * quad points */
  if (testDerivOrder() == 0 && unkDerivOrder() == 0)
  {
    double* aPtr = &((*A)[0]);
    double* coeffPtr = (double*) coeff;
    int offset = 0 ;

    if (globalCurve().isCurveValid())
    {        /* ---------- ACI logic ----------- */

    	  Array<double> quadWeightsTmp = quadWeights_[0];
    	  Array<Point> quadPointsTmp = quadPts_[0];
    	  bool isCutByCurve;

    	  for (int c=0; c<JVol.numCells(); c++, offset+=nNodes())
    	  {
    	      double detJ = fabs(JVol.detJ()[c]);
    	      int fc = 0;
    	      if (nFacetCases() != 1) fc = facetIndex[c];
    	      const Array<double>& w = W_[fc];

    	      quadWeightsTmp = quadWeights_[fc];
    	      quadPointsTmp = quadPts_[fc];
    	      /* call the special integration routine */
    	      quad_.getAdaptedWeights(cellType(), dim(), (*cellLID)[c] , fc ,
    	    		  mesh() , globalCurve() , quadPointsTmp , quadWeightsTmp , isCutByCurve );
    	      if (isCutByCurve){
    	    	  Array<double> wi;
    	    	  wi.resize(nQuad() * nNodesTest() *nNodesUnk() ); //recalculate the special weights
    	    	  for (int ii = 0 ; ii < wi.size() ; ii++ ) wi[ii] = 0.0;
    	    	  for (int nt = 0 ; nt < nNodesTest() ; nt++){
    	    		  for (int nu=0; nu<nNodesUnk(); nu++)
    	    			  for (int q=0 ; q < quadWeightsTmp.size() ; q++)
    	    				  //Indexing: unkNode + nNodesUnk()*(testNode + nNodesTest()*(unkDerivDir + nRefDerivUnk()*(testDerivDir + nRefDerivTest()*q)))
    	    				  wi[nu + nNodesUnk()*(nt + nNodesTest()*q)] +=
    	    						  chop(quadWeightsTmp[q] * W_ACI_F2_[fc][q][0][nt][0][nu]);
    	    	  }
    	    	  for (int q=0; q<nQuad(); q++, coeffPtr++){
    	    		  double f = (*coeffPtr)*detJ;
    	    		  for (int n=0; n<nNodes(); n++){
    	    			  aPtr[offset+n] += f*wi[n + nNodes()*q];
    	    		  }
    	    	  }
    	      }// end isCutByCurve
    	      else
    	      {
    	    	  for (int q=0; q<nQuad(); q++, coeffPtr++){
    	    		  double f = (*coeffPtr)*detJ;
    	    		  for (int n=0; n<nNodes(); n++){
    	    			  aPtr[offset+n] += f*w[n + nNodes()*q];
    	    		  }
    	    	  }
    	      }
    	  }
      }
      else         /* ---------- NO ACI logic ----------- */
      {

    	  for (int c=0; c<JVol.numCells(); c++, offset+=nNodes())
    	  {
    	      double detJ = fabs(JVol.detJ()[c]);
    	      int fc = 0;
    	      if (nFacetCases() != 1) fc = facetIndex[c];
    	      const Array<double>& w = W_[fc];
    	      for (int q=0; q<nQuad(); q++, coeffPtr++)
    	      {
    	    	  double f = (*coeffPtr)*detJ;
    	    	  for (int n=0; n<nNodes(); n++)
    	    	  {
    	    		  aPtr[offset+n] += f*w[n + nNodes()*q];
    	    	  }
    	      }
    	  }
      }
    addFlops( JVol.numCells() * (1 + nQuad() * (1 + 2*nNodes())) );
  }
  else
  {
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
        
      
    if (useSumFirstMethod())
    {
      transformSummingFirst(JTrans.numCells(), facetIndex, cellLIDs, GPtr, coeff, A);
    }
    else
    {
      transformSummingLast(JTrans.numCells(), facetIndex, cellLIDs, GPtr, coeff, A);
    }
  }
}

void QuadratureIntegral
::transformSummingFirst(int nCells,
  const Array<int>& facetIndex,
  const RCP<Array<int> >& cellLIDs,
  const double* const GPtr,
  const double* const coeff,
  RCP<Array<double> >& A) const
{
  double* aPtr = &((*A)[0]);
  double* coeffPtr = (double*) coeff;
  const Array<int>* cellLID = cellLIDs.get();
  
  int transSize = 0; 
  if (order()==2)
  {
    transSize = nRefDerivTest() * nRefDerivUnk();
  }
  else
  {
    transSize = nRefDerivTest();
  }

  /* The sum workspace is used to store the sum of untransformed quantities */
  static Array<double> sumWorkspace;

  int swSize = transSize * nNodes();
  sumWorkspace.resize(swSize);
  
  
  /*
   * The number of operations for the sum-first method is 
   * 
   * Adds: (N_c * nNodes * transSize) * (N_q + 1) 
   * Multiplies: same as number of adds
   * Total: 2*(N_c * nNodes * transSize) * (N_q + 1) 
   */
  

    if (globalCurve().isCurveValid())
    {       /* ---------- ACI logic ----------- */

    	 Array<double> quadWeightsTmp = quadWeights_[0];
    	 Array<Point> quadPointsTmp = quadPts_[0];
    	 bool isCutByCurve;

    	for (int c=0; c<nCells; c++)
    	{
    	    /* sum untransformed basis combinations over quad points */
    	    for (int i=0; i<swSize; i++) sumWorkspace[i]=0.0;
    	    int fc = 0;
    	    if (nFacetCases() > 1) fc = facetIndex[c];

    	    const Array<double>& w = W_[fc];

    	    quadWeightsTmp = quadWeights_[fc];
    	    quadPointsTmp = quadPts_[fc];
    	    /* call the special integration routine */
    	    quad_.getAdaptedWeights(cellType(), dim(), (*cellLID)[c] , fc ,
    	    		mesh() , globalCurve() , quadPointsTmp , quadWeightsTmp , isCutByCurve );
    	    if (isCutByCurve){
    	    	Array<double> wi;
    	    	if (order()==1){ // one form integral
    	    		wi.resize(nQuad() * nRefDerivTest() * nNodesTest()); //recalculate the special weights
    	    		for (int ii = 0 ; ii < wi.size() ; ii++ ) wi[ii] = 0.0;
    	    		for (int n = 0 ; n < nNodes() ; n++)
    	    			for (int t=0; t<nRefDerivTest(); t++)
    	    				for (int q=0 ; q < quadWeightsTmp.size() ; q++)
    	    					//Indexing: testNode + nNodesTest()*(testDerivDir + nRefDerivTest()*q)
    	    					wi[n + nNodes()*(t + nRefDerivTest()*q)] +=
    	    							chop(quadWeightsTmp[q] * W_ACI_F1_[fc][q][t][n]);
    	    	}else{ // two from integrals
    	    		wi.resize(nQuad() * nRefDerivTest() * nNodesTest()
    	    				* nRefDerivUnk() * nNodesUnk() ); //recalculate the special weights
    	    		for (int ii = 0 ; ii < wi.size() ; ii++ ) wi[ii] = 0.0;

    	    		for (int t=0; t<nRefDerivTest(); t++)
    	    			for (int nt = 0 ; nt < nNodesTest() ; nt++){
    	    				for (int u=0; u<nRefDerivUnk(); u++)
    	    					for (int nu=0; nu<nNodesUnk(); nu++)
    	    						for (int q=0 ; q < quadWeightsTmp.size() ; q++)
    	    							//Indexing: unkNode + nNodesUnk()*(testNode + nNodesTest()*
    	    							//(unkDerivDir + nRefDerivUnk()*(testDerivDir + nRefDerivTest()*q)))
    	    							wi[nu + nNodesUnk()*(nt + nNodesTest()*(u + nRefDerivUnk()*(t + nRefDerivTest()*q)))] +=
    	    									chop(quadWeightsTmp[q] * W_ACI_F2_[fc][q][t][nt][u][nu]);
    	    			}
    	    	}// end two form integral
    	    	for (int q=0; q<nQuad(); q++, coeffPtr++){
    	    		double f = (*coeffPtr);
    	    		for (int n=0; n<swSize; n++){
    	    			sumWorkspace[n] += f*wi[n + q*swSize];
    	    		}
    	    	}
    	    }// end cut by curve
    	    else{
    	    	for (int q=0; q<nQuad(); q++, coeffPtr++)
    	    	{
    	    		double f = (*coeffPtr);
    	    		for (int n=0; n<swSize; n++)
    	    		{
    	    			sumWorkspace[n] += f*w[n + q*swSize];
    	    		}
    	    	}
    	    }
    	}// end from for loop over cells
    }
    else
    {      	    /* ---------- NO ACI logic ----------- */
    	for (int c=0; c<nCells; c++)
    	{
    	    /* sum untransformed basis combinations over quad points */
    	    for (int i=0; i<swSize; i++) sumWorkspace[i]=0.0;
    	    int fc = 0;
    	    if (nFacetCases() > 1) fc = facetIndex[c];

    	    const Array<double>& w = W_[fc];

    	    for (int q=0; q<nQuad(); q++, coeffPtr++)
    	    {
    	    	double f = (*coeffPtr);
    	    	for (int n=0; n<swSize; n++)
    	    	{
    	    		sumWorkspace[n] += f*w[n + q*swSize];
    	    	}
    	    }

    	    /* transform the sum */
    	    const double * const gCell = &(GPtr[transSize*c]);
    	    double* aCell = aPtr + nNodes()*c;
    	    for (int i=0; i<nNodes(); i++)
    	    {
    	    	for (int j=0; j<transSize; j++)
    	    	{
    	    		aCell[i] += sumWorkspace[nNodes()*j + i] * gCell[j];
    	    	}
    	    }
    	}// for loop over cells
    }
  int flops = 2*(nCells * nNodes() * transSize) * (nQuad() + 1) ;
  addFlops(flops);
}

void QuadratureIntegral
::transformSummingLast(int nCells,
  const Array<int>& facetIndex,
  const RCP<Array<int> >& cellLIDs,
  const double* const GPtr,
  const double* const coeff,
  RCP<Array<double> >& A) const
{
  double* aPtr = &((*A)[0]);
  int transSize = 0; 
  const Array<int>* cellLID = cellLIDs.get();
  
  if (order()==2)
  {
    transSize = nRefDerivTest() * nRefDerivUnk();
  }
  else
  {
    transSize = nRefDerivTest();
  }

  /* This workspace is used to store the jacobian values scaled by the coeff
   * at that quad point */
  static Array<double> jWorkspace;
  jWorkspace.resize(transSize);


  /*
   * The number of operations for the sum-last method is 
   * 
   * Adds: N_c * N_q * transSize * nNodes
   * Multiplies: numAdds + N_c * N_q * transSize
   * Total: N_c * N_q * transSize * (1 +  2*nNodes)
   */

    if (globalCurve().isCurveValid()){
        /* ---------- ACI logic ----------- */
    	Array<double> quadWeightsTmp = quadWeights_[0];
    	Array<Point> quadPointsTmp = quadPts_[0];
    	bool isCutByCurve;

    	for (int c=0; c<nCells; c++)
    	{
    	    const double* const gCell = &(GPtr[transSize*c]);
    	    double* aCell = aPtr + nNodes()*c;
    	    int fc = 0;
    	    if (nFacetCases() > 1) fc = facetIndex[c];
    	    const Array<double>& w = W_[fc];

    	    quadWeightsTmp = quadWeights_[fc];
    	    quadPointsTmp = quadPts_[fc];
    	    /* call the special integration routine */
    	    quad_.getAdaptedWeights(cellType(), dim(), (*cellLID)[c] , fc ,
    	    		mesh() , globalCurve() , quadPointsTmp , quadWeightsTmp , isCutByCurve );
    	    if (isCutByCurve){
    	    	Array<double> wi;
    	    	if (order()==1){ // one form integral
    	    		wi.resize(nQuad() * nRefDerivTest() * nNodesTest()); //recalculate the special weights
    	    		for (int ii = 0 ; ii < wi.size() ; ii++ ) wi[ii] = 0.0;
    	    		for (int n = 0 ; n < nNodes() ; n++)
    	    			for (int t=0; t<nRefDerivTest(); t++)
    	    				for (int q=0 ; q < quadWeightsTmp.size() ; q++)
    	    					//Indexing: testNode + nNodesTest()*(testDerivDir + nRefDerivTest()*q)
    	    					wi[n + nNodes()*(t + nRefDerivTest()*q)] +=
    	    							chop(quadWeightsTmp[q] * W_ACI_F1_[fc][q][t][n]);
    	    	}else{ // two from integrals
    	    		wi.resize(nQuad() * nRefDerivTest() * nNodesTest()
    	    				* nRefDerivUnk() * nNodesUnk() ); //recalculate the special weights
    	    		for (int ii = 0 ; ii < wi.size() ; ii++ ) wi[ii] = 0.0;

    	    		for (int t=0; t<nRefDerivTest(); t++)
    	    			for (int nt = 0 ; nt < nNodesTest() ; nt++){
    	    				for (int u=0; u<nRefDerivUnk(); u++)
    	    					for (int nu=0; nu<nNodesUnk(); nu++)
    	    						for (int q=0 ; q < quadWeightsTmp.size() ; q++)
    	    							//Indexing: unkNode + nNodesUnk()*(testNode + nNodesTest()*
    	    							//(unkDerivDir + nRefDerivUnk()*(testDerivDir + nRefDerivTest()*q)))
    	    							wi[nu + nNodesUnk()*(nt + nNodesTest()*(u + nRefDerivUnk()*(t + nRefDerivTest()*q)))] +=
    	    									chop(quadWeightsTmp[q] * W_ACI_F2_[fc][q][t][nt][u][nu]);
    	    			}
    	    	}// end two form integral
    	    	for (int q=0; q<nQuad(); q++){
    	    		double f = coeff[c*nQuad() + q];
    	    		for (int t=0; t<transSize; t++) jWorkspace[t]=f*gCell[t];

    	    		for (int n=0; n<nNodes(); n++){
    	    			for (int t=0; t<transSize; t++){
    	    				aCell[n] += jWorkspace[t]*wi[n + nNodes()*(t + transSize*q)];
    	    			}
    	    		}
    	    	}
    	    }// end cut by curve
    	    else{
    	    	for (int q=0; q<nQuad(); q++){
    	    		double f = coeff[c*nQuad() + q];
    	    		for (int t=0; t<transSize; t++) jWorkspace[t]=f*gCell[t];

    	    		for (int n=0; n<nNodes(); n++){
    	    			for (int t=0; t<transSize; t++){
    	    				aCell[n] += jWorkspace[t]*w[n + nNodes()*(t + transSize*q)];
    	    			}
    	    		}
    	    	}
    	    }
		}// loop over cells
    }
    else       /* ---------- NO ACI logic ----------- */
    {
       	for (int c=0; c<nCells; c++)
        {
        	    const double* const gCell = &(GPtr[transSize*c]);
        	    double* aCell = aPtr + nNodes()*c;
        	    int fc = 0;
        	    if (nFacetCases() > 1) fc = facetIndex[c];
        	    const Array<double>& w = W_[fc];
        	    for (int q=0; q<nQuad(); q++)
        	    {
        	    	double f = coeff[c*nQuad() + q];
        	    	for (int t=0; t<transSize; t++) jWorkspace[t]=f*gCell[t];

        	    	for (int n=0; n<nNodes(); n++)
        	    	{
        	    		for (int t=0; t<transSize; t++)
        	    		{
        	    			aCell[n] += jWorkspace[t]*w[n + nNodes()*(t + transSize*q)];
        	    		}
        	    	}
        	    }
        }
  }
  int flops = nCells * nQuad() * transSize * (1 + 2*nNodes()) ;
  addFlops(flops);
}
