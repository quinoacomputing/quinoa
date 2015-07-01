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

#include "SundanceCurveQuadratureIntegral.hpp"

#include "SundanceCurveIntegralCalc.hpp"
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

static Time& maxCellQuadratureTimer() 
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("max cell quadrature"); 
  return *rtn;
}


CurveQuadratureIntegral::CurveQuadratureIntegral(
  const CellType& cellType,
  const bool isConstantIntegral,
  const QuadratureFamily& quad,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh,
  int verb)
  : ElementIntegral(dimension(cellType), cellType, dimension(cellType),
    cellType, true, globalCurve, mesh, verb),
    quad_(quad),
    quadPts_(),
    quadWeights_(),
    W_(),
    useSumFirstMethod_(true),
    useConstCoeff_(isConstantIntegral)
{
  Tabs tab0(0);
  
  TEUCHOS_TEST_FOR_EXCEPTION( cellType != mesh.cellType(mesh.spatialDim()) , std::runtime_error, " CurveQuadratureIntegral::CurveQuadratureIntegral , only on MAXcell!!");

  SUNDANCE_MSG1(setupVerb(), tab0 << "CurveQuadratureIntegral ctor for 0-form");
  if (setupVerb()) describe(Out::os());
  
  SUNDANCE_MSG1(setupVerb(), tab0 << "quadrature family=" << quad);  

  /* create the quad points and weights */
  quad.getPoints(mesh.cellType(mesh.spatialDim()-1), quadPts_, quadWeights_);
  int nQuad = quadPts_.size();
      
  W_.resize(nQuad);
  quadCurveDerivs_.resize(nQuad);
  quadCurveNormals_.resize(nQuad);
      
  for (int q=0; q<nQuad; q++)
  {
    W_[q] = quadWeights_[q];
  }
}


CurveQuadratureIntegral::CurveQuadratureIntegral(
  const CellType& cellType,
  const bool isConstantIntegral,
  const BasisFamily& testBasis,
  int alpha,
  int testDerivOrder,
  const QuadratureFamily& quad,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh,
  int verb)
  : ElementIntegral(dimension(cellType), cellType, dimension(cellType),
    cellType, testBasis,
    alpha, testDerivOrder, true, globalCurve, mesh, verb),
    quad_(quad),
    quadPts_(),
    quadWeights_(),
    W_(),
    useSumFirstMethod_(true),
    useConstCoeff_(isConstantIntegral)
{
  Tabs tab0(0);
  
  TEUCHOS_TEST_FOR_EXCEPTION( cellType != mesh.cellType(mesh.spatialDim()) , std::runtime_error, " CurveQuadratureIntegral::CurveQuadratureIntegral , only on MAXcell!!");

  SUNDANCE_MSG1(setupVerb(), tab0 << "CurveQuadratureIntegral ctor for 1-form");
  if (setupVerb()) describe(Out::os());
  assertLinearForm();

  
  SUNDANCE_MSG1(setupVerb(), tab0 << "quadrature family=" << quad);  
  
  quad.getPoints(mesh.cellType(mesh.spatialDim()-1), quadPts_, quadWeights_);
  int nQuad = quadPts_.size();
      
  W_.resize(nQuad * nRefDerivTest() * nNodesTest());
  quadCurveDerivs_.resize(nQuad);
  quadCurveNormals_.resize(nQuad);

  SUNDANCE_MSG1(setupVerb(), tab0 << "num nodes for test function " << nNodesTest());
}




CurveQuadratureIntegral::CurveQuadratureIntegral(
  const CellType& cellType,
  const bool isConstantIntegral,
  const BasisFamily& testBasis,
  int alpha,
  int testDerivOrder,
  const BasisFamily& unkBasis,
  int beta,
  int unkDerivOrder,
  const QuadratureFamily& quad,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh,
  int verb)
  : ElementIntegral(dimension(cellType), cellType, dimension(cellType),
    cellType, testBasis,
    alpha, testDerivOrder, unkBasis, beta, unkDerivOrder, true,
    globalCurve, mesh, 
    verb),
    quad_(quad),
    quadPts_(),
    quadWeights_(),
    W_(),
    useSumFirstMethod_(true),
    useConstCoeff_(isConstantIntegral)
{
  Tabs tab0(0);

  TEUCHOS_TEST_FOR_EXCEPTION( cellType != mesh.cellType(mesh.spatialDim()) , std::runtime_error, " CurveQuadratureIntegral::CurveQuadratureIntegral , only on MAXcell!!");
  
  SUNDANCE_MSG1(setupVerb(), tab0 << "CurveQuadratureIntegral ctor for 2-form");
  if (setupVerb()) describe(Out::os());
  assertBilinearForm();

  // store the quadrature points and weights  
  quad.getPoints(mesh.cellType(mesh.spatialDim()-1), quadPts_, quadWeights_);
  int nQuad = quadPts_.size();

  W_.resize(nQuad * nRefDerivTest() * nNodesTest()  
    * nRefDerivUnk() * nNodesUnk());

  quadCurveDerivs_.resize(nQuad);
  quadCurveNormals_.resize(nQuad);
}



void CurveQuadratureIntegral::updateRefCellInformation(int maxCellLID , const ParametrizedCurve& curve) const{

	//get the points on the reference cell, instead of "quadPts_" ,  and update "quadCurveDerivs_"
	Tabs tabs;
    SUNDANCE_MSG3(integrationVerb(), tabs << "CurveQuadratureIntegral::updateRefCellInformation, "
  		  "Update Cell LID: " << maxCellLID << "  curve.myID():" << curve.myID());

	if ( mesh().hasCurvePoints( maxCellLID , curve.myID() ))
	{
	    SUNDANCE_MSG3(integrationVerb(), tabs << "CurveQuadratureIntegral::updateRefCellInformation, has Points, just returning them");
		mesh().getCurvePoints( maxCellLID , curve.myID() , quadPts_ , quadCurveDerivs_ , quadCurveNormals_ );
	}
	else // we have to calculate now the points
	{
		// calculate the intersection points
	    SUNDANCE_MSG3(integrationVerb(), tabs << "CurveQuadratureIntegral::updateRefCellInformation, computing the intersection points");
		CurveIntegralCalc::getCurveQuadPoints( mesh().cellType(mesh().spatialDim()) , maxCellLID , mesh() , curve ,
				quad_ , quadPts_ , quadCurveDerivs_ , quadCurveNormals_ );
		// store the intersection point in the mesh
	    SUNDANCE_MSG3(integrationVerb(), tabs << "CurveQuadratureIntegral::updateRefCellInformation, storing the intersection points in the mesh");
		mesh().setCurvePoints( maxCellLID , curve.myID() , quadPts_ , quadCurveDerivs_ , quadCurveNormals_ );
	}
}


void CurveQuadratureIntegral::updateRefCellIntegralOneForm(int maxCellLID , int cellInBatch) const {

	  /* compute the test functions */
	  Array<Array<Array<Array<double> > > > testBasisVals(nRefDerivTest());

	  Tabs tabs;
      SUNDANCE_MSG2(integrationVerb(), tabs << "CurveQuadratureIntegral::updateRefCellIntegralTwoForm, "
    		  "Update Cell LID: " << maxCellLID);

	  // get the points on the reference cell, instead of "quadPts_" ,  and update "quadCurveDerivs_"
	  updateRefCellInformation( maxCellLID , globalCurve() );

	  for (int r=0; r<nRefDerivTest(); r++)
	  {
	    Tabs tab3;
	    SUNDANCE_MSG1(setupVerb(), tab3
	      << "evaluating basis functions for ref deriv direction " << r);
	    MultiIndex mi;
	    testBasisVals[r].resize(testBasis().dim());
	    if (testDerivOrder()==1) mi[r] = 1;
	    SpatialDerivSpecifier deriv(mi);
	    testBasis().refEval(evalCellType(), quadPts_, deriv,
	      testBasisVals[r], setupVerb());
	  }

	  int vecComp = 0;
	  int nQuad = quadPts_.size();
	  for (int q=0; q<nQuad; q++)
	  {
	    for (int t=0; t<nRefDerivTest(); t++)
	    {
	      for (int nt=0; nt<nNodesTest(); nt++)
	      {
	        wValue(q, t, nt)
	          = chop(quadWeights_[q] * testBasisVals[t][vecComp][q][nt]) ;
	      }
	    }
	  }

      // chop the numerical zero values
	  for (int i=0; i<W_.size(); i++) W_[i] = chop(W_[i]);
}


void CurveQuadratureIntegral::updateRefCellIntegralTwoForm(int maxCellLID , int cellInBatch) const {

	  /* compute the basis functions */
	  Array<Array<Array<Array<double> > > > testBasisVals(nRefDerivTest());
	  Array<Array<Array<Array<double> > > > unkBasisVals(nRefDerivUnk());

	  Tabs tabs;
      SUNDANCE_MSG2(integrationVerb(), tabs << "CurveQuadratureIntegral::updateRefCellIntegralTwoForm, "
    		  "Update Cell LID: " << maxCellLID);

	  // get the points on the reference cell, instead of "quadPts_" and update "quadCurveDerivs_"
	  updateRefCellInformation( maxCellLID , globalCurve() );

	  for (int r=0; r<nRefDerivTest(); r++)
	  {
	    testBasisVals[r].resize(testBasis().dim());
	    MultiIndex mi;
	    if (testDerivOrder()==1) mi[r] = 1;
	    SpatialDerivSpecifier deriv(mi);
	    SUNDANCE_MSG1(setupVerb(), tabs
	      << "evaluating basis functions for ref deriv direction " << r);
	    testBasis().refEval(evalCellType(), quadPts_, deriv,
	      testBasisVals[r], setupVerb());
	  }

	  for (int r=0; r<nRefDerivUnk(); r++)
	  {
	    unkBasisVals[r].resize(unkBasis().dim());
	    MultiIndex mi;
	    if (unkDerivOrder()==1) mi[r] = 1;
	    SpatialDerivSpecifier deriv(mi);
	    SUNDANCE_MSG1(setupVerb(), tabs
	      << "evaluating basis functions for ref deriv direction " << r);
	    unkBasis().refEval(evalCellType(),
	      quadPts_, deriv, unkBasisVals[r], setupVerb());
	  }


	  int vecComp = 0;
	  int nQuad = quadPts_.size();
	  /* form the products of basis functions at each quad pt */
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
	            wValue(q, t, nt, u, nu)
	              = chop(quadWeights_[q] * testBasisVals[t][vecComp][q][nt]
	                * unkBasisVals[u][vecComp][q][nu]);
	          }
	        }
	      }
	    }
	  }

      // chop the numerical zero values
	  for (int i=0; i<W_.size(); i++) W_[i] = chop(W_[i]);
}

void CurveQuadratureIntegral
::transformZeroForm(const CellJacobianBatch& JTrans,  
  const CellJacobianBatch& JVol,
  const Array<int>& isLocalFlag,
  const Array<int>& facetIndex,
  const RCP<Array<int> >& cellLIDs,
  const double constCoeff,
  const double* const coeff,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(maxCellQuadratureTimer());
  Tabs tabs;
  SUNDANCE_MSG1(integrationVerb(), tabs << "doing zero form by quadrature , isLocalFlag.size():" << isLocalFlag.size() );

  TEUCHOS_TEST_FOR_EXCEPTION(order() != 0, std::logic_error,
    "CurveQuadratureIntegral::transformZeroForm() called "
    "for form of order " << order());

  int nQuad = quadWeights_.size();

  // if we have constant coefficients then copy the constant value in the array
  double* coeffPtr = (double*) coeff;
  Array<double> constCoeff_arr( JVol.numCells() * nQuad);
  if (useConstCoeff_){
	  for (int jk=0; jk < JVol.numCells() * nQuad ; jk++ ) { constCoeff_arr[jk] = constCoeff; };
	  coeffPtr = &(constCoeff_arr[0]);
  }

  double& a = (*A)[0];
  SUNDANCE_MSG5(integrationVerb(), tabs << "input A=");
  if (integrationVerb() >= 5) writeTable(Out::os(), tabs, *A, 6);
  const Array<double>& w = W_;
  double curveDerivNorm = 0.0;

  if ( (int)isLocalFlag.size() == 0 )
  {
	 for (int c=0; c<JVol.numCells(); c++)
	 {
  	   // get the points on the reference cell, instead of "quadPts_" ,  and update "quadCurveDerivs_"
	   updateRefCellInformation( (*cellLIDs)[c] , globalCurve() );

       // here we do not have to update W_ since it is only the weight of the quadrature points
	   for (int q=0; q<nQuad; q++, coeffPtr++)
	   {
	   	  // no multiplication with the Jacobian!!!
	   	  // multiply with the norm of the derivative of the curve (this should be stored in the mesh)
	   	  curveDerivNorm = sqrt(quadCurveDerivs_[q]*quadCurveDerivs_[q]);
	      a += w[q]*(*coeffPtr)*curveDerivNorm;
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

    for (int c=0; c<JVol.numCells(); c++)
    {
      // only if this cell is local on the processor
  	  if ( isLocalFlag[c] )
  	  {
  	    // get the points on the reference cell, instead of "quadPts_" ,  and update "quadCurveDerivs_"
  	    updateRefCellInformation( (*cellLIDs)[c] , globalCurve() );

  	    // here we pass the evaluated points(coeffPtr) for the interface (polygon or 3D surface), in case to set the values on the interface
  	    // in this place it should be only evaluated for functional, since isLocalFlag.size() > 0
  	    SUNDANCE_MSG1(integrationVerb(), " call .addEvaluationPointValues for cellLID: " << (*cellLIDs)[c] << " curveID:" << globalCurve().myID() );
  	    globalCurve().addEvaluationPointValues( mesh() , (*cellLIDs)[c] , nQuad , coeffPtr , quadPts_ );
      
        // here we do not have to update W_ since it is only the weight of the quadrature points
        for (int q=0; q<nQuad; q++, coeffPtr++)
        {
    	  // no multiplication with the Jacobian!!!
    	  // multiply with the norm of the derivative of the curve (this should be stored in the mesh)
    	  curveDerivNorm = sqrt(quadCurveDerivs_[q]*quadCurveDerivs_[q]);
          a += w[q]*(*coeffPtr)*curveDerivNorm;
        }
  	  }
  	  else {
  		 // increment the coeff vector, if this cell does not count
  		 coeffPtr += nQuad;
  	  }
    }
  }

  SUNDANCE_MSG5(integrationVerb(), tabs << "output A = ");
  if (integrationVerb() >= 5) writeTable(Out::os(), tabs, *A, 6);

  SUNDANCE_MSG1(integrationVerb(), tabs << "done zero form");
}


void CurveQuadratureIntegral::transformOneForm(const CellJacobianBatch& JTrans,
  const CellJacobianBatch& JVol,
  const Array<int>& facetIndex,
  const RCP<Array<int> >& cellLIDs,
  const double constCoeff,
  const double* const coeff,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(maxCellQuadratureTimer());
  Tabs tabs;
  TEUCHOS_TEST_FOR_EXCEPTION(order() != 1, std::logic_error,
    "CurveQuadratureIntegral::transformOneForm() called for form "
    "of order " << order());
  SUNDANCE_MSG2(integrationVerb(), tabs << "doing one form by quadrature");
  int flops = 0;

  int nQuad = quadWeights_.size();

  /* If the derivative order is zero, the only thing to be done 
   * is to multiply by the cell's Jacobian determinant and sum over the
   * quad points */
  if (testDerivOrder() == 0)
  {
    double* aPtr = &((*A)[0]);
    SUNDANCE_MSG5(integrationVerb(), tabs << "input A = ");
    if (integrationVerb() >= 5) writeTable(Out::os(), tabs, *A, 6);
  
    int offset = 0 ;
    const Array<double>& w = W_;
    double* coeffPtr = (double*) coeff;
    double curveDerivNorm = 0.0 , f = 0.0;
    // if we have constant coefficient then copy the constant value into one array
    Array<double> constCoeff_arr( JVol.numCells() * nQuad);
    if (useConstCoeff_){
  	  for (int jk=0; jk < JVol.numCells() * nQuad ; jk++ ) { constCoeff_arr[jk] = constCoeff; };
  	  coeffPtr = &(constCoeff_arr[0]);
    }

    for (int c=0; c<JVol.numCells(); c++, offset+=nNodes())
    {
    	Tabs tab2;
        SUNDANCE_MSG4(integrationVerb(), tab2 << "c=" << c << " detJ=" << JVol.detJ()[c] );

        // update W_
        updateRefCellIntegralOneForm( (*cellLIDs)[c] , c );

        for (int q=0; q<nQuad; q++, coeffPtr++)
        {
          Tabs tab3;
          // multiply with the norm of the derivative , (not with the determinant)
          curveDerivNorm = sqrt(quadCurveDerivs_[q]*quadCurveDerivs_[q]);
          f = (*coeffPtr) * curveDerivNorm;
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

        if (integrationVerb() >= 4)
        {
          Out::os() << tab2 << "integration results on cell:" << std::endl;
          Out::os() << tab2 << setw(10) << "n" << setw(12) << "I_n" << std::endl;
          for (int n=0; n<nNodes(); n++) 
          {
            Out::os() << tab2 << setw(10) << n 
                      << setw(12) << setprecision(6) << aPtr[offset+n] << std::endl;
          }
        }
        
    } // -- end cell iteration

    SUNDANCE_MSG5(integrationVerb(), tabs << "output A = ");
    if (integrationVerb() >= 5) writeTable(Out::os(), tabs, *A, 6);
  }
  else
  {
    /* If the derivative order is nonzero, then we have to do a transformation. */
    
    createOneFormTransformationMatrix(JTrans, JVol);
    
    SUNDANCE_MSG4(transformVerb(), 
      Tabs() << "transformation matrix=" << G(alpha()));
    
    double* GPtr = &(G(alpha())[0]);      
    
    transformSummingFirst(JVol.numCells(), JVol ,facetIndex, cellLIDs, constCoeff , GPtr, coeff, A);

  }
  addFlops(flops);
}


void CurveQuadratureIntegral::transformTwoForm(const CellJacobianBatch& JTrans,
  const CellJacobianBatch& JVol,
  const Array<int>& facetIndex,
  const RCP<Array<int> >& cellLIDs,
  const double constCoeff,
  const double* const coeff,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(maxCellQuadratureTimer());
  Tabs tabs;
  TEUCHOS_TEST_FOR_EXCEPTION(order() != 2, std::logic_error,
    "CurveQuadratureIntegral::transformTwoForm() called for form "
    "of order " << order());
  SUNDANCE_MSG2(integrationVerb(), tabs << "doing one form by quadrature");

  int nQuad = quadWeights_.size();

  /* If the derivative orders are zero, the only thing to be done 
   * is to multiply by the cell's Jacobian determinant and sum over the
   * quad points */
  if (testDerivOrder() == 0 && unkDerivOrder() == 0)
  {
    double* aPtr = &((*A)[0]);
    int offset = 0 ;
    const Array<double>& w = W_;
    double* coeffPtr = (double*) coeff;
    double curveDerivNorm = 0.0 , f = 0.0;

    // if we have constant coefficient then copy the constant value into one array
    Array<double> constCoeff_arr( JVol.numCells() * nQuad);
    SUNDANCE_MSG2(integrationVerb(), tabs << "CurveQuadratureIntegral::transformTwoForm , use const variabl: " << useConstCoeff_ );
    if (useConstCoeff_){
      SUNDANCE_MSG2(integrationVerb(), tabs << "CurveQuadratureIntegral::transformTwoForm , "
    		  "Using constant coeff: " << constCoeff);
  	  for (int jk=0; jk < JVol.numCells() * nQuad ; jk++ ) { constCoeff_arr[jk] = constCoeff; };
  	  coeffPtr = &(constCoeff_arr[0]);
    }

    for (int c=0; c<JVol.numCells(); c++, offset+=nNodes())
    {
        // ---- update W_ ---
        updateRefCellIntegralTwoForm( (*cellLIDs)[c] , c );

        for (int q=0; q<nQuad; q++, coeffPtr++)
        {
          // multiply with the norm of the derivative , (not with the determinant)
          curveDerivNorm = sqrt(quadCurveDerivs_[q]*quadCurveDerivs_[q]);
          f = (*coeffPtr)*curveDerivNorm;
          for (int n=0; n<nNodes(); n++)
          {
            aPtr[offset+n] += f*w[n + nNodes()*q];
          }
        }
    } // -- end cell iteration

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
        
    transformSummingFirst(JTrans.numCells(), JVol , facetIndex, cellLIDs, constCoeff , GPtr, coeff, A);

  }
}

void CurveQuadratureIntegral
::transformSummingFirst(int nCells,
  const CellJacobianBatch& JVol,
  const Array<int>& facetIndex,
  const RCP<Array<int> >& cellLIDs,
  const double constCoeff,
  const double* const GPtr,
  const double* const coeff,
  RCP<Array<double> >& A) const
{
  double* aPtr = &((*A)[0]);
  int nQuad = quadWeights_.size();
  double* coeffPtr = (double*) coeff;

  // if we have constant coefficient then copy the constant value into one array
  Array<double> constCoeff_arr( JVol.numCells() * nQuad);
  if (useConstCoeff_){
	  for (int jk=0; jk < JVol.numCells() * nQuad ; jk++ ) { constCoeff_arr[jk] = constCoeff; };
	  coeffPtr = &(constCoeff_arr[0]);
  }
  
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
  const Array<double>& w = W_;
  double curveDerivNorm = 0.0 , f = 0.0;
  
  /*
   * The number of operations for the sum-first method is 
   * 
   * Adds: (N_c * nNodes * transSize) * (N_q + 1) 
   * Multiplies: same as number of adds
   * Total: 2*(N_c * nNodes * transSize) * (N_q + 1) 
   */

    /* Sum */
    for (int c=0; c<nCells; c++)
    {
      /* sum untransformed basis combinations over quad points */
      for (int i=0; i<swSize; i++) sumWorkspace[i]=0.0;
      
      // ---- update W_ ----
      if (order() == 2) {
    	  updateRefCellIntegralTwoForm( (*cellLIDs)[c] , c );
      }
      if (order() == 1) {
    	  updateRefCellIntegralOneForm( (*cellLIDs)[c] , c );
      }

      for (int q=0; q<nQuad; q++, coeffPtr++)
      {
        // multiply with the norm of the derivative ,(not with the determinant)
      	curveDerivNorm = sqrt(quadCurveDerivs_[q]*quadCurveDerivs_[q]);
        f = (*coeffPtr);
        for (int n=0; n<swSize; n++)
        {
          sumWorkspace[n] += f*w[n + q*swSize] * curveDerivNorm;
        }
      }
    
      /* transform the sum */
      const double * const gCell = &(GPtr[transSize*c]);
      double* aCell = aPtr + nNodes()*c;
  	  // since the trafo matrix is multiplied by the determinant, we have to multiply this by 1/det(J)
      double detJInv = 1/fabs(JVol.detJ()[c]);
      for (int i=0; i<nNodes(); i++)
      {
        for (int j=0; j<transSize; j++)
        {
          aCell[i] += sumWorkspace[nNodes()*j + i] * gCell[j] * detJInv;
        }
      }
    }
}

