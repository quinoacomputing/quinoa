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

#include "SundanceCurveEvalMediator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceTempStack.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceCurveNormExpr.hpp"
#include "SundanceCellVectorExpr.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceCurveIntegralCalc.hpp"

#include "Teuchos_BLAS.hpp"


using namespace Sundance;
using namespace Teuchos;
using namespace Playa;

TEUCHOS_TIMER(coordEvalTimer, "Quad mediator: coord eval")

CurveEvalMediator::CurveEvalMediator(
		const Mesh& mesh, 	const ParametrizedCurve& paramcurve ,
		int cellDim, const QuadratureFamily& quad) : StdFwkEvalMediator(mesh, cellDim) ,
    numEvaluationCases_(-1),
    quad_(quad),
    numQuadPtsForMaxCell_(0) ,
    paramcurve_(paramcurve) {

	// set the cell types
	maxCellType_ = mesh.cellType(mesh.spatialDim() );
	curveCellType_ = mesh.cellType(mesh.spatialDim()-1);

	// detect the (constant) quadrature points (we must have the same nr quad points per cell!!!)
	Array<Point> quadPoints;
	Array<double> quadWeights;
	quad.getPoints(curveCellType_ , quadPoints , quadWeights );
	numQuadPtsForMaxCell_ = quadWeights.size();
}

Time& CurveEvalMediator::coordEvaluationTimer()
{
  return coordEvalTimer();
}

void CurveEvalMediator::setCellType(const CellType& cellType,
  const CellType& maxCellType,
  bool isInternBdry) 
{
  StdFwkEvalMediator::setCellType(cellType, maxCellType, isInternBdry);

  maxCellType_ = maxCellType;

  Tabs tab;

  SUNDANCE_MSG2(verb(), tab <<  "CurveEvalMediator::setCellType: cellType="
    << cellType << " cellDim=" << cellDim() << " maxCellType=" << maxCellType);
  SUNDANCE_MSG2(verb(), tab << "integration spec: =" 
    << integrationCellSpec());
  SUNDANCE_MSG2(verb(), tab << "forbid cofacet integrations =" 
    << forbidCofacetIntegrations());
  if (isInternalBdry()) 
  {
    SUNDANCE_MSG2(verb(), tab << "working on internal boundary");
  }


  if (cellType != maxCellType)
  {
    // this should not happen
  }
  else 
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb(), tab1 << "no need for facet cases; work on original cell"); 
    numEvaluationCases_ = 1;
  }

  SUNDANCE_MSG2(verb(), tab << "num eval cases =" << numEvaluationCases_); 

}

int CurveEvalMediator::numQuadPts(const CellType& ct) const
{
	if (ct != maxCellType_){
		// throw error
		TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "CurveEvalMediator::numQuadPts , cellType must be maxCellType_:" << ct);
	}
	return numQuadPtsForMaxCell_;
}

void CurveEvalMediator::evalCellDiameterExpr(const CellDiameterExpr* expr,
  RCP<EvalVector>& vec) const 
{
  // this is the same for curve integral as for QuadratureEvalMediator
  Tabs tabs;
  SUNDANCE_MSG2(verb(),tabs 
    << "CurveEvalMediator evaluating cell diameter expr "
    << expr->toString());

  int nQuad = numQuadPts(cellType());
  int nCells = cellLID()->size();

  SUNDANCE_MSG3(verb(),tabs << "number of quad pts=" << nQuad);
  Array<double> diameters;
  mesh().getCellDiameters(cellDim(), *cellLID(), diameters);

  vec->resize(nQuad*nCells);
  double * const xx = vec->start();
  int k=0;
  for (int c=0; c<nCells; c++)
  {
    double h = diameters[c];
    for (int q=0; q<nQuad; q++, k++) 
    {
      xx[k] = h;
    }
  }
}

void CurveEvalMediator::evalCellVectorExpr(const CellVectorExpr* expr,
  RCP<EvalVector>& vec) const 
{
  Tabs tabs;
  SUNDANCE_MSG2(verb(),tabs 
    << "CurveEvalMediator evaluating cell vector expr "
    << expr->toString());

  int nQuad = numQuadPts(cellType());
  int nCells = cellLID()->size();

  vec->resize(nQuad*nCells);

  SUNDANCE_MSG3(verb(),tabs << "number of quad pts=" << nQuad);
  int dir = expr->componentIndex();

  Array<Point> vectors;
  if (expr->isNormal())
  { 
    mesh().outwardNormals(*cellLID(), vectors);
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(cellDim() != 1, std::runtime_error,
      "unable to compute tangent vectors for cell dim = " << cellDim());
    mesh().tangentsToEdges(*cellLID(), vectors);
  }
    
  double * const xx = vec->start();
  int k=0;
  for (int c=0; c<nCells; c++)
  {
    double n = vectors[c][dir];
    for (int q=0; q<nQuad; q++, k++) 
    {
      xx[k] = n;
    }
  }
}

void CurveEvalMediator::evalCoordExpr(const CoordExpr* expr,
  RCP<EvalVector>& vec) const 
{
  Tabs tabs;
  SUNDANCE_MSG2(verb(),tabs
    << "CurveEvalMediator evaluating coord expr "
    << expr->toString());

  TimeMonitor timer(coordEvalTimer());
  
  int nQuad = numQuadPts(cellType());
  int nCells = cellLID()->size();
  int d = expr->dir();
  
  SUNDANCE_MSG3(verb(),tabs << "number of quad pts=" << nQuad << ", nCells:" << nCells);

  // the size must be correct!!
  vec->resize(nQuad*nCells);
  double * const xx = vec->start();
  int k=0;
  int maxDim = mesh().spatialDim() ;
  Array<int> cellLIDs_tmp(1);
  Array<Point> phyPoints;
  Array<Point> refPoints;
  Array<Point> refDevs;
  Array<Point> refNormal;

  // iterate over each cell
  for (int c=0; c<nCells; c++)
  {
	int maxCellLID = (*cellLID())[c];
	cellLIDs_tmp[0] = (*cellLID())[c];
	// see if the mesh already contains the intersection points
	if ( mesh().hasCurvePoints( maxCellLID , paramcurve_.myID() ))
	{
		mesh().getCurvePoints( maxCellLID , paramcurve_.myID() , refPoints , refDevs , refNormal );
	}
	else // we have to calculate now the points
	{
		// calculate the intersection points
		CurveIntegralCalc::getCurveQuadPoints( maxCellType_ , maxCellLID , mesh() , paramcurve_ , quad_ ,
				refPoints, refDevs , refNormal);
		// store the intersection point in the mesh
		mesh().setCurvePoints( maxCellLID , paramcurve_.myID() , refPoints, refDevs , refNormal );
	}

	// calculate the physical coordinates
	mesh().pushForward( maxDim , cellLIDs_tmp , refPoints , phyPoints );

	// ---- return the desired component of the intersection point ----
    for (int q=0; q<nQuad; q++, k++)
    {
      xx[k] = phyPoints[q][d];
    }
  }
}


void CurveEvalMediator::evalCurveNormExpr(const CurveNormExpr* expr,
  RCP<EvalVector>& vec) const {

	  Tabs tabs;
	  SUNDANCE_MSG2(verb(),tabs
	    << "CurveEvalMediator evaluating curve normal expr "
	    << expr->toString());

	  int nQuad = numQuadPts(cellType());
	  int nCells = cellLID()->size();
	  int d = expr->dir();

	  SUNDANCE_MSG3(verb(),tabs << "number of quad pts=" << nQuad);

	  // the size of this vector must be correct !!!
	  vec->resize(nCells*nQuad);

	  double * const xx = vec->start();
	  int k=0;
	  Array<int> cellLIDs_tmp(1);
	  Array<Point> phyPoints;
	  Array<Point> refPoints;
	  Array<Point> refDevs;
	  Array<Point> refNormal;

	  // iterate over each cell
	  for (int c=0; c<nCells; c++)
	  {
		int maxCellLID = (*cellLID())[c];
		cellLIDs_tmp[0] = (*cellLID())[c];
		// see if the mesh already contains the intersection points
		if ( mesh().hasCurvePoints( maxCellLID , paramcurve_.myID() ))
		{
			mesh().getCurvePoints( maxCellLID , paramcurve_.myID() , refPoints , refDevs , refNormal );
		}
		else // we have to calculate now the points
		{
			// calculate the intersection points
			CurveIntegralCalc::getCurveQuadPoints( maxCellType_ , maxCellLID , mesh() , paramcurve_ , quad_ ,
					refPoints, refDevs , refNormal);
			// store the intersection point in the mesh
			mesh().setCurvePoints( maxCellLID , paramcurve_.myID() , refPoints, refDevs , refNormal );
		}

		// ---- return the desired component of the unit normal vector ----
	    for (int q=0; q<nQuad; q++, k++)
	    {
	      xx[k] = refNormal[q][d];
	    }
	  }
}


void CurveEvalMediator
::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
  const Array<MultiIndex>& multiIndices,
  Array<RCP<EvalVector> >& vec) const
{

	  int verbo = dfVerb();
	  Tabs tab1;
	  SUNDANCE_MSG2(verbo , tab1
	    << "CurveEvalMediator evaluating Discrete Function expr " << expr->toString());

	  const DiscreteFunctionData* f = DiscreteFunctionData::getData(expr);

	  TEUCHOS_TEST_FOR_EXCEPTION(f==0, std::logic_error,
	    "QuadratureEvalMediator::evalDiscreteFuncElement() called "
	    "with expr that is not a discrete function");

	  SUNDANCE_MSG2(verbo , tab1 << "After casting DiscreteFunctionData" << expr->toString());

	  RCP<Array<Array<double> > > localValues;
	  Array<int> cellLIDs_tmp(1);
	  Array<Point> phyPoints;
	  Array<Point> refPoints;
	  Array<Point> refDevs;
	  Array<Point> refNormal;
	  int nCells = cellLID()->size();

      RCP<const MapStructure> mapStruct;
	  int myIndex = expr->myIndex();
	  int nQuad = numQuadPtsForMaxCell_;
	  Array<int> k(multiIndices.size(),0);

	  Teuchos::BLAS<int,double> blas;

	  SUNDANCE_MSG2(verbo , tab1 << "After declaring BLAS: " << expr->toString());

	  // resize correctly the result vector
	  for (int i=0; i<multiIndices.size(); i++)
	  {
		  vec[i]->resize(nCells*nQuad);
	  }

	  // loop over each cell
	  for (int c=0; c<nCells; c++)
	  {
			int maxCellLID = (*cellLID())[c];
		    localValues = rcp(new Array<Array<double> >());
			SUNDANCE_MSG2(verbo , tab1 <<  "Cell:" << c <<  " of " << nCells << " , maxCellLID:" << maxCellLID );

			cellLIDs_tmp[0] = (*cellLID())[c];
			SUNDANCE_MSG2(verbo , tab1 <<  " Before calling f->getLocalValues:" << cellLIDs_tmp.size() <<
					" f==0 : " << (f==0) << " tmp:" << f->mesh().spatialDim());
			// - get local values from the DiscreteFunctionElementData
			mapStruct = f->getLocalValues(maxCellDim(), cellLIDs_tmp , *localValues);
			SUNDANCE_MSG2(verbo , tab1 <<  " After getting mapStruct:" << maxCellLID );
			SUNDANCE_MSG2(verbo , tab1 <<  " mapStruct->numBasisChunks():" << mapStruct->numBasisChunks() );
		    int chunk = mapStruct->chunkForFuncID(myIndex);
		    int funcIndex = mapStruct->indexForFuncID(myIndex);
		    int nFuncs = mapStruct->numFuncs(chunk);
			SUNDANCE_MSG2(verbo , tab1 <<  " chunk:" << chunk );
			SUNDANCE_MSG2(verbo , tab1 <<  " funcIndex:" << funcIndex );
			SUNDANCE_MSG2(verbo , tab1 <<  " nFuncs:" << nFuncs );

			// the chunk of the function
		    BasisFamily basis = rcp_dynamic_cast<BasisFamilyBase>(mapStruct->basis(chunk));
		    int nNodesTotal = basis.nReferenceDOFsWithFacets(maxCellType(), maxCellType());

			// - get intersection (reference)points from the mesh (if not existent than compute them)
			if ( mesh().hasCurvePoints( maxCellLID , paramcurve_.myID() ))
			{
				mesh().getCurvePoints( maxCellLID , paramcurve_.myID() , refPoints , refDevs , refNormal );
			}
			else // we have to calculate now the points
			{
				// calculate the intersection points
				CurveIntegralCalc::getCurveQuadPoints( maxCellType_ , maxCellLID , mesh() , paramcurve_ , quad_ ,
						refPoints, refDevs , refNormal);
				// store the intersection point in the mesh
				mesh().setCurvePoints( maxCellLID , paramcurve_.myID() , refPoints, refDevs , refNormal );
			}

			// loop over each multi-index
			SUNDANCE_MSG2(verbo , tab1 << " multiIndices.size()" << multiIndices.size() );
			for (int i=0; i<multiIndices.size(); i++)
			{

				int nDerivResults = 1;
				if ( multiIndices[i].order() == 1 ) nDerivResults = maxCellDim();

			    int pDir = 0;
			    int derivNum = 1;
				MultiIndex mi;
				SUNDANCE_MSG2(verbo , tab1 << " before asking anything i = " << i);
				SUNDANCE_MSG2(verbo , tab1 << " multiindex order : " << multiIndices[i].order());
				SUNDANCE_MSG2(verbo , tab1 << " multiindex : " << multiIndices[i] );
				if (multiIndices[i].order() > 0){
					pDir = multiIndices[i].firstOrderDirection();
					mi[pDir] = 1;
					derivNum = mesh().spatialDim();
				}

				Array<Array<double> > result(nQuad*derivNum);
				Array<Array<Array<double> > > tmp;

				int offs = nNodesTotal;

				// resize the result vector
				for (int deriv = 0 ; deriv < derivNum ; deriv++)
				{
                    // test weather we have to compute derivative
					if (multiIndices[i].order() > 0){
						// in case of derivatives we set one dimension
						MultiIndex mi_tmp;
						mi_tmp[deriv] = 1;
						SpatialDerivSpecifier deriv(mi_tmp);
						SUNDANCE_MSG2(verbo , tab1 << "computing derivatives : " << deriv << " on reference cell ");
						basis.refEval( maxCellType_ , refPoints , deriv, tmp , verbo );
					} else {
						SpatialDerivSpecifier deriv(mi);
						 // --- end eval basis functions
						SUNDANCE_MSG2(verbo , tab1 << "computing values reference cell ");
						basis.refEval( maxCellType_ , refPoints , deriv, tmp , verbo );
					}

					SUNDANCE_MSG2(verbo , tab1 << "resize result vector , offs:" << offs);
					for (int q=0; q<nQuad; q++){
						result[nQuad*deriv + q].resize(offs);
					}

				   	// copy the result in an other format
					SUNDANCE_MSG2(verbo , tab1 << "copy results ");
					int  offs1 = 0;
				    for (int q=0; q<nQuad; q++)
				    {
				    	offs1 = 0;
						for (int d=0; d<basis.dim(); d++)
						{
						   int nNodes = tmp[d][q].size();
					       for (int n=0; n<nNodes; n++ , offs1++ )
					       {
					          	result[nQuad*deriv + q][offs1] = tmp[d][q][n];
					       }
					    }
					}
				}// loop over all dimensional derivative

                // multiply the local results with the coefficients, (matrix vector OP)
				SUNDANCE_MSG2(verbo , tab1 << "summing up values , funcIndex:" << funcIndex << " offs:" << offs);
				for (int deriv = 0 ; deriv < derivNum ; deriv++)
				{
					for (int q=0; q<nQuad; q++)
					{
						double sum = 0.0;
						// sum over nodes
						for (int n = 0 ; n < offs ; n++){
							sum = sum + result[nQuad*deriv + q][n] * (*localValues)[chunk][funcIndex*offs + n];
						}
						// sum up the result in the 0th element
						result[nQuad*deriv + q][0] = sum;
					}
				}

			    // multiply the result if necesary with the inverse of the Jacobian
				const CellJacobianBatch& J = JTrans();
				if (mi.order()==1)
				{
					Tabs tab1;
					Tabs tab2;
				    SUNDANCE_MSG2(verbo, tab2 << "Jacobian batch nCells=" << J.numCells());
			        SUNDANCE_MSG2(verbo, tab2 << "Jacobian batch cell dim=" << J.cellDim());
			        SUNDANCE_MSG2(verbo, tab2 << "Jacobian batch spatial dim=" << J.spatialDim());
                    // we just multiply the derivative direction component
			        Array<double> invJ;
			        J.getInvJ(c, invJ);
					for (int q=0; q<nQuad; q++)
					{
						double sum = 0.0;
						for (int deriv = 0 ; deriv < derivNum ; deriv++)
						{
							// multiply one row from the J^{-T} matrix with the gradient vector
							sum = sum + result[nQuad*deriv + q][0] * invJ[derivNum*pDir + deriv];
						}
						// the resulting derivative on the physical cell in the "pDir" direction
						result[q][0] = sum;
					}
			    }

				// --- just copy the result to the "vec" back, the result should be in the  "result[q][0]" place----
				//SUNDANCE_MSG2(verbo , tab1 << "copy results back ");
				double* vecPtr = vec[i]->start();
				for (int q=0; q<nQuad; q++, k[i]++)
				{
					vecPtr[k[i]] = result[q][0];
				}
				SUNDANCE_MSG2(verbo , tab1 << " END copy results back ");
			} // --- end loop multiindex
			SUNDANCE_MSG2(verbo , tab1 << " END loop over multiindex ");
	}// --- end loop over cells
	SUNDANCE_MSG2(verbo , tab1 << " END loop over cells ");
}


void CurveEvalMediator::print(std::ostream& os) const
{
	// todo: implement this
}


