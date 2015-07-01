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


/*
 * SundanceTrapesoidQuadrature.cpp
 *
 *  Created on: Jan 20, 2011
 *      Author: benk
 */

#include "SundanceTrapesoidQuadrature.hpp"

using namespace Sundance;

TrapesoidQuadrature::TrapesoidQuadrature(int resolution) :
   QuadratureFamilyBase(1) , resolution_(resolution) {
	// this quadrature has always order 1 , the minimum number of points
	resolution_ = (resolution_<2) ? 2 : resolution_;
}

XMLObject TrapesoidQuadrature::toXML() const
{
	XMLObject rtn("TrapesoidQuadrature");
	rtn.addAttribute("order", Teuchos::toString(order()));
	return rtn;
}

void TrapesoidQuadrature::getLineRule(
		Array<Point>& quadPoints,
		Array<double>& quadWeights) const {

	int nrVirInnerPoints = resolution_ - 1;
	double innerWeights = 1/((double)nrVirInnerPoints);

	int nrPoints = nrVirInnerPoints + 1;
	quadPoints.resize(nrPoints);
	quadWeights.resize(nrPoints);

	// create the first point the beginning
	quadPoints[0] = Point(0.0);
	quadWeights[0] = innerWeights/2.0;

	for (int tmp = 1 ; tmp < nrPoints-1 ; tmp++){
		quadPoints[tmp] = Point( (double)tmp/(double)(nrPoints-1) );
		quadWeights[tmp] = innerWeights;
	}

	// create the last point the end of the line
	quadPoints[nrPoints-1] = Point(1.0);
	quadWeights[nrPoints-1] = innerWeights/2.0;

}


void TrapesoidQuadrature::getTriangleRule(
		Array<Point>& quadPoints,
		Array<double>& quadWeights) const {
   // todo: implement this
	SUNDANCE_ERROR("Triangle rule not available for " << toXML());
}


void TrapesoidQuadrature::getQuadRule(
		Array<Point>& quadPoints,
		Array<double>& quadWeights) const {

	Array<Point> quadPointsLine;
	Array<double> quadWeightsLine;
	// get the line rule
    this->getLineRule( quadPointsLine, quadWeightsLine );

    int nrPointPerAxis = quadPointsLine.size();
    // we simply take the tensor product
    quadPoints.resize(nrPointPerAxis*nrPointPerAxis);
    quadWeights.resize(nrPointPerAxis*nrPointPerAxis);

    int pcount = 0;
    for (int ix = 0 ; ix < nrPointPerAxis ; ix ++){
    	for (int iy = 0 ; iy < nrPointPerAxis ; iy ++){
    		// here we take the tensor product of the
    		quadPoints[pcount] = Point( quadPointsLine[ix][0] , quadPointsLine[iy][0] );
    		quadWeights[pcount] = quadWeightsLine[ix] * quadWeightsLine[iy];
    		pcount++;
    	}
    }

}


void TrapesoidQuadrature::getTetRule(
		Array<Point>& quadPoints,
		Array<double>& quadWeights) const {
   // todo: implement this
   SUNDANCE_ERROR("Tet rule not available for " << toXML());
}


void TrapesoidQuadrature::getBrickRule(
		Array<Point>& quadPoints,
		Array<double>& quadWeights) const {

	Array<Point> quadPointsLine;
	Array<double> quadWeightsLine;
	// get the line rule
    getLineRule( quadPointsLine, quadWeightsLine );

    int nrPointPerAxis = quadPointsLine.size();
    // we simply take the tensor product
    quadPoints.resize(nrPointPerAxis*nrPointPerAxis*nrPointPerAxis);
    quadWeights.resize(nrPointPerAxis*nrPointPerAxis*nrPointPerAxis);

    int pcount = 0;
    for (int ix = 0 ; ix < nrPointPerAxis ; ix ++){
    	for (int iy = 0 ; iy < nrPointPerAxis ; iy ++){
    		for (int iz = 0 ; iz < nrPointPerAxis ; iz ++){
    			// here we take the tensor product of the
    			quadPoints[pcount] = Point( quadPointsLine[ix][0] , quadPointsLine[iy][0] , quadPointsLine[iz][0]);
    			quadWeights[pcount] = quadWeightsLine[ix] * quadWeightsLine[iy] * quadWeightsLine[iz];
    			pcount++;
    		}
    	}
    }
}


void TrapesoidQuadrature::getAdaptedWeights(
		const CellType& cellType, int cellDim,
		int cellLID, int facetIndex, const Mesh& mesh,
		const ParametrizedCurve& globalCurve,
		Array<Point>& quadPoints, Array<double>& quadWeights,
		bool &weightsChanged) const {

	weightsChanged = true; // todo: this not always might be true to save computations we might test if this is the case
	// get the weights from the mesh
	if (mesh.IsSpecialWeightValid() && mesh.hasSpecialWeight(cellDim, cellLID))
	{
			mesh.getSpecialWeight(cellDim, cellLID, quadWeights);
			//SUNDANCE_MSG3(verb, tabs << "GaussianQuadrature::getAdaptedWeights Found cached weights for cell LID " << cellLID)
			return;
	}

	// get the points
	if (quadPoints.size() <= 1) getPoints(cellType,  quadPoints, quadWeights);

	// we say that the cell is always cut by the curve so these weights will be used
	weightsChanged = true;
	Array<Point> tmpPoints = quadPoints;

    Array<int> celLIDs(1);
    celLIDs[0] = cellLID;
    // transform the ref points to real coordinates
	mesh.pushForward( cellDim, celLIDs, quadPoints, tmpPoints );

	quadWeights.resize(tmpPoints.size());
    // simple weight calculation, by just multiplying the weight by the integration factor(from the curve)
	for (int i=0 ; i < tmpPoints.size() ; i++){
		quadWeights[i] = quadWeights[i] * globalCurve.integrationParameter(tmpPoints[i]);
	}

	// store the weights in the mesh
	mesh.setSpecialWeight(cellDim, cellLID, quadWeights);
}
