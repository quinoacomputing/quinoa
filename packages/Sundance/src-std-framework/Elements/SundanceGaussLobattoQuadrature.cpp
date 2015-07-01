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
 * SundanceGaussLobattoQuadrature.cpp
 *
 *  Created on: Jan 20, 2011
 *      Author: benk
 */

#include "SundanceGaussLobattoQuadrature.hpp"
#include "SundanceCurveIntegralCalc.hpp"
#include "SundancePolygon2D.hpp"
#include "SundanceSurf3DCalc.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;


int const GaussLobattoQuadrature::quadsEdgesPoints[4][2] = { {0,1} , {0,2} , {1,3} , {2,3} };
                                                       //  0   1   2   3   4   5   6   7   8   9  10  11
int const GaussLobattoQuadrature::edge3DProjection[12] = { 0 , 1 , 2 , 1 , 2 , 0 , 2 , 2 , 0 , 1 , 1 , 0};


GaussLobattoQuadrature::GaussLobattoQuadrature(int order) :
   QuadratureFamilyBase(order) {
	nrPointin1D_ = order+1;
	verb_ = 0;
}

XMLObject GaussLobattoQuadrature::toXML() const
{
	XMLObject rtn("GaussLobattoQuadrature");
	rtn.addAttribute("order", Teuchos::toString(order()));
	return rtn;
}


void GaussLobattoQuadrature::getTriangleRule(
		Array<Point>& quadPoints,
		Array<double>& quadWeights) const {
    // todo: implement this
	SUNDANCE_ERROR("Triangle rule not available for " << toXML());
}

void GaussLobattoQuadrature::getQuadRule(
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


void GaussLobattoQuadrature::getTetRule(
		Array<Point>& quadPoints,
		Array<double>& quadWeights) const {
   // todo: implement this
   SUNDANCE_ERROR("Tet rule not available for " << toXML());
}


void GaussLobattoQuadrature::getBrickRule(
		Array<Point>& quadPoints,
		Array<double>& quadWeights) const {
	Array<Point> quadPointsLine;
	Array<double> quadWeightsLine;
	// get the line rule
    this->getLineRule( quadPointsLine, quadWeightsLine );

    int nrPointPerAxis = quadPointsLine.size();
    // we simply take the tensor product
    quadPoints.resize(nrPointPerAxis*nrPointPerAxis*nrPointPerAxis);
    quadWeights.resize(nrPointPerAxis*nrPointPerAxis*nrPointPerAxis);

    int pcount = 0;
    for (int iz = 0 ; iz < nrPointPerAxis ; iz ++){
	  for (int iy = 0 ; iy < nrPointPerAxis ; iy ++){
        for (int ix = 0 ; ix < nrPointPerAxis ; ix ++){
    		// here we take the tensor product of the
    		quadPoints[pcount] = Point( quadPointsLine[ix][0] , quadPointsLine[iy][0] , quadPointsLine[iz][0]);
    		quadWeights[pcount] = quadWeightsLine[ix] * quadWeightsLine[iy] * quadWeightsLine[iz];
    		pcount++;
    	}
      }
    }
}


void GaussLobattoQuadrature::getAdaptedWeights(
		const CellType& cellType, int cellDim,
		int cellLID, int facetIndex, const Mesh& mesh,
		const ParametrizedCurve& globalCurve,
		Array<Point>& quadPoints, Array<double>& quadWeights,
		bool &weightsChanged) const {

	Tabs tabs;

	weightsChanged = true; // todo: this not always might be true to save computations we might test if this is the case
	if (mesh.IsSpecialWeightValid() && mesh.hasSpecialWeight(cellDim, cellLID))
	{
		mesh.getSpecialWeight(cellDim, cellLID, quadWeights);
		//SUNDANCE_MSG3(verb, tabs << "GaussianQuadrature::getAdaptedWeights Found cached weights for cell LID " << cellLID)
		return;
	}
	// if we have no quad points then get them
	if (quadPoints.size() <= 1) getPoints(cellType,  quadPoints, quadWeights);

	// Maximal cell type
	CellType maxCellType = mesh.cellType(mesh.spatialDim());

	switch (maxCellType)
	{
	// We have a triangle based mesh
	case TriangleCell:
	{
		switch (cellType)
		{
		case TriangleCell:
		{
			// todo: implement this
			break;
		}
		default:
#ifndef TRILINOS_7
			SUNDANCE_ERROR ("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for submaximal cell " << maxCellType << " in a triangle mesh")
			;
#else
			SUNDANCE_ERROR7("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for submaximal cell " << maxCellType << " in a triangle mesh");
#endif
		}
		break;
	}
	case QuadCell:
	{
		switch(cellType)
		{
		case QuadCell:
		{
			// call the method to integrate the Quad

			if (CurveIntegralCalc::usePolygoneCurveQuadrature( cellType , mesh , globalCurve)){
				// if the curve is a polygon then use a different
				getAdaptedQuadWeights_polygon(cellLID, mesh, globalCurve, quadPoints, quadWeights, weightsChanged);
			}
			else
			{
				getAdaptedQuadWeights(cellLID, mesh, globalCurve, quadPoints, quadWeights, weightsChanged);
			}
			break;
		}
		default:
#ifndef TRILINOS_7
		SUNDANCE_ERROR ("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for submaximal cell " << maxCellType << " in a triangle mesh")
		;
#else
		SUNDANCE_ERROR7("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for submaximal cell " << maxCellType << " in a triangle mesh");
#endif
		}
		break;
	}
	case BrickCell:
	{
		switch(cellType)
		{
		case BrickCell:
		{
			// call the method which does the job for 3D
			getAdaptedQuadWeights_surf(cellLID, mesh, globalCurve, quadPoints, quadWeights, weightsChanged);
			break;
		}
		default: {}
	    }
		break;
	}
	default:
#ifndef TRILINOS_7
		SUNDANCE_ERROR("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for cell type " << cellType)
		;
#else
		SUNDANCE_ERROR7("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for cell type " << cellType);
#endif
	}

	//store the weights in the mesh
	mesh.setSpecialWeight(cellDim, cellLID, quadWeights);
}



void GaussLobattoQuadrature::getAdaptedQuadWeights_polygon(int cellLID, const Mesh& mesh,
		const ParametrizedCurve& globalCurve, Array<Point>& quadPoints,
		Array<double>& quadWeights, bool& weightsChanged) const {

	// this works only in 2D for polygons, it is important that we ask the underlying object and not the object itself
	const Polygon2D* poly = dynamic_cast<const Polygon2D* > ( globalCurve.ptr().get()->getUnderlyingCurveObj() );
	if ( poly == 0 ) {
		TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error , " getAdaptedQuadWeights_polygon , paramCurve must be a Polygon and the cell must be a quad" );
		return;
	}

	// the main idea of the algorithm is to get a raw of points
	// (P1,P2,P3,P4,P5) and (alpha1,aplpha2), the coefficients
	// this raw of points must be ordered correctly and contains intersection and polygon points
	// by integrating the left side of this curve we can calculate the exact integral

	// - get the intersection points
	// - get the polygon points
	// form out of these points a row of points with the shortest distance
	// choose as starting point the one with the smallest x and y coordinate

	// use this raw of points for quadrature
	//  - see in which direction is always increasing the coordinate
	//  - choose one direction (x or y) to parse the raw of points (if not always increasing then take the sign in consideration)

	Array<Point> quadOfPoints(4);
	Array<Point> intesectPoints(8);
	Array<Point> pointStack(0);
	Array<bool>  isInnerPoint(0);
	Array<Point> pointStackTmp(0);
	Array<Point> allPolygonPoints(0);
	int nrIPoints , pointsNr , polyNrPoint ;
	double eps = 1e-14;

	SUNDANCE_MSG3(verb_, " ---------- START METHOD -----------  ");

	// get the points of the quad
	for (int i = 0 ; i < 4 ; i++){
		quadOfPoints[i] = mesh.nodePosition( mesh.facetLID(2,cellLID , 0 , i , nrIPoints ) );
	}

	// get all the intersection points
	pointsNr = 0;
	for (int i = 0 ; i < 4 ; i++)
	{
		// for each edge test intersection point
		globalCurve.returnIntersectPoints( quadOfPoints[quadsEdgesPoints[i][0]] , quadOfPoints[quadsEdgesPoints[i][1]],
				 nrIPoints , intesectPoints );
		//SUNDANCE_MSG3(verb_, " Test edge: "<< quadOfPoints[quadsEdgesPoints[i][0]] << " -- " << quadOfPoints[quadsEdgesPoints[i][1]] );
		for (int j = 0 ; j < nrIPoints ; j++){
			bool addThisPoint = true;
			// add these points only when they are not duplicate
			for (int tmpP = 0 ; tmpP < pointsNr ; tmpP++){
				double dist = ::sqrt( (pointStackTmp[tmpP] - intesectPoints[j])*(pointStackTmp[tmpP] - intesectPoints[j]) );
				if ( dist < eps) addThisPoint = false;
			}
			// if the point is not a duplicate then add
			if (addThisPoint) {
				pointStack.resize( pointsNr + 1  );
				pointStackTmp.resize( pointsNr + 1 );
				isInnerPoint.resize( pointsNr + 1 );
				// transform to reference coordinates back
				pointStack[pointsNr] = Point( ( intesectPoints[j][0] - quadOfPoints[0][0])/( quadOfPoints[3][0]-quadOfPoints[0][0] ) ,
					                      ( intesectPoints[j][1] - quadOfPoints[0][1])/( quadOfPoints[3][1]-quadOfPoints[0][1] )  );
				// we store the intersection point also in real coordinates
				pointStackTmp[pointsNr] = intesectPoints[j];
				isInnerPoint[pointsNr] = false;
				//SUNDANCE_MSG3(verb_, " adding Ints Points:"<< pointsNr << " P:"
				//		<< pointStack[pointsNr] << " Real P:" << intesectPoints[j] );
				pointsNr = pointsNr + 1;
			}
		}
	}

	// get the point
	poly->getCellsPolygonesPoint( cellLID ,allPolygonPoints);
	polyNrPoint = allPolygonPoints.size();
	//SUNDANCE_MSG3(verb_, " Polygon has :"<< polyNrPoint << " points inside the cell");
	// add also the points which come from the polygon
	for (int j = 0 ; j < polyNrPoint ; j++){
		bool addThisPoint = true;
		// add these points only when they are not duplicate
		for (int tmpP = 0 ; tmpP < pointsNr ; tmpP++){
			double dist = ::sqrt( (pointStackTmp[tmpP] - allPolygonPoints[j])*(pointStackTmp[tmpP] - allPolygonPoints[j]) );
			if ( dist < eps) addThisPoint = false;
		}
		// if the point is not a duplicate then add
		if (addThisPoint) {
			pointStack.resize( pointsNr + 1  );
			pointStackTmp.resize( pointsNr + 1 );
			isInnerPoint.resize( pointsNr + 1 );
			// transform to reference coordinates back
			pointStack[pointsNr] = Point( ( allPolygonPoints[j][0] - quadOfPoints[0][0])/( quadOfPoints[3][0]-quadOfPoints[0][0] ) ,
				                      ( allPolygonPoints[j][1] - quadOfPoints[0][1])/( quadOfPoints[3][1]-quadOfPoints[0][1] )  );
			// we store the intersection point also in real coordinates
			pointStackTmp[pointsNr] = allPolygonPoints[j];
			isInnerPoint[pointsNr] = true;
			//SUNDANCE_MSG3(verb_, " adding Polygon Points:"<< pointsNr << " P:"
			//		<< pointStack[pointsNr] << " Real P:" << allPolygonPoints[j] );
			pointsNr = pointsNr + 1;
		}
	}

	// now we have all the points we just need to get them in a sorted way, after distance (Manhattan distance)
	// first we select the point which is closes to the (0,0) Point and is an intersection point of the edges
	double dist = 1e+100 , firstind = -1;
	for (int tmpP = 0 ; tmpP < pointsNr ; tmpP++){
		// we look only for intersection points to have as starting point
		if ( (pointStack[tmpP][0] + pointStack[tmpP][1] < dist) && (!isInnerPoint[tmpP])){
			dist = pointStack[tmpP][0] + pointStack[tmpP][1];
			firstind = tmpP;
		}
	}

	Array<Point> orderedPoints(1);
	Array<Point> orderedPointsRealCoords(1);
	Array<bool> orderedPointsIsIntersection(1);
	orderedPoints[0] = pointStack[firstind];
	orderedPointsRealCoords[0] = pointStackTmp[firstind];
	orderedPointsIsIntersection[0] = (!isInnerPoint[firstind]);
	//SUNDANCE_MSG3(verb_, "Add first point to position 0  from firstind:"<< firstind << " , P:" << orderedPoints[0]);

	int orderedP = 1;

	// with this for loop we create an X or Y ordered list of points, we take the
	// first point and insert the next one there, where the distance to the neughbors is minimal
	for (int tmpP = 0 ; tmpP < pointsNr ; tmpP++){
		// the first point which we added we do not want to add again
		if (tmpP == firstind) continue;

		// add all the 1..pointsNr points to the ordered list
		double dist = 1e+100 , dist_tmp;
		int ind = -1;
		//SUNDANCE_MSG3(verb_, " Trying to insert: "<< tmpP << " out of " << pointsNr);
		for (int p = 0 ; p <= orderedP ; p++){
			dist_tmp = dist;
			//SUNDANCE_MSG3(verb_, " probe position " << p << " , orderedP:"<<orderedP);
			if (( p == orderedP ) || (p == 0) ){
				if (!isInnerPoint[tmpP]){
					Point poin = orderedPoints[(p == orderedP)?p-1:p] - pointStack[tmpP];
					// check this, what kind of distance * 1.0 should be here added
					dist_tmp = 1.0 * ::sqrt(poin*poin);
				}
			}
			else {
				Point poin1 = orderedPoints[p-1] - pointStack[tmpP];
				Point poin2 = orderedPoints[p] - pointStack[tmpP];
				Point prev_dist = orderedPoints[p] - orderedPoints[p-1];
				dist_tmp = ::sqrt(poin1*poin1) + ::sqrt(poin2*poin2) - ::sqrt(prev_dist*prev_dist);
			}
			//SUNDANCE_MSG3(verb_, " dist_tmp: " << dist_tmp << " , dist:" << dist );
			// if this distance is shorter than mark this position
			if (dist_tmp < dist){
				dist = dist_tmp;
				ind = p;
			}
		}

		// the first point should be inserted always to the back
		if ( (ind == 0) && ( orderedP <= 1) ){ ind = 1;}

		// if "ind"
		TEUCHOS_TEST_FOR_EXCEPTION( ind < 0 , std::runtime_error , " GaussLobattoQuadrature::getAdaptedQuadWeights_polygon ind < 0 ind:"<<ind );
		// insert the point to the "ind" position
		orderedPoints.resize(orderedP+1);
		orderedPointsRealCoords.resize(orderedP+1);
		orderedPointsIsIntersection.resize(orderedP+1);
		for (int p = orderedP-1 ; p >= ind ; p--){
			orderedPoints[p+1] = orderedPoints[p];
			orderedPointsRealCoords[p+1] = orderedPointsRealCoords[p];
			orderedPointsIsIntersection[p+1] = orderedPointsIsIntersection[p];
		}
		//SUNDANCE_MSG3(verb_, "Add point to position ind:"<< ind <<  " , tmpP:" << tmpP << " , orderedP:" << orderedP );
		//SUNDANCE_MSG3(verb_, "pointStack[tmpP]:"<< pointStack[tmpP] <<  " , pointStackTmp[tmpP]:" << pointStackTmp[tmpP] );
		orderedPoints[ind] = pointStack[tmpP];
		orderedPointsRealCoords[ind] = pointStackTmp[tmpP];
		orderedPointsIsIntersection[ind] = (!isInnerPoint[tmpP]);
		orderedP = orderedP + 1;
	}

	// at this stage we have the ordered raw of point we just have to use it
	// the ordered
	bool increaseX = true;
	bool increaseY = true;
	double minX = orderedPoints[0][0], maxX = orderedPoints[0][0], minY = orderedPoints[0][1], maxY = orderedPoints[0][1];
	for (int p = 1 ; p < orderedP ; p++){
		if ( orderedPoints[p-1][0] - eps > orderedPoints[p][0]) increaseX = false;
		if ( orderedPoints[p-1][1] - eps > orderedPoints[p][1]) increaseY = false;
		minX = (orderedPoints[p][0] < minX) ? orderedPoints[p][0] : minX;
		maxX = (orderedPoints[p][0] > maxX) ? orderedPoints[p][0] : maxX;
		minY = (orderedPoints[p][1] < minY) ? orderedPoints[p][1] : minY;
		maxY = (orderedPoints[p][1] > maxY) ? orderedPoints[p][1] : maxY;
	}

	// at this point we know in which dimension we are going to parse the points (X or Y)
	// we have to make sure that the points are incrementally ordered with respect to this dimension
	//(or at least the with respect to the first and last element)
	if ( (increaseX) || (increaseY == false)) { // ordered in X direction
		if ( orderedPoints[0][0] > orderedPoints[orderedP-1][0]+eps) {
			SUNDANCE_MSG3(verb_, " Order in X direction ");
			for (int p = 0 ; p < (orderedP/2) ; p++){ // swap the three elements
				Point tmp_pp = orderedPoints[p];    orderedPoints[p] = orderedPoints[orderedP-1-p];
				orderedPoints[orderedP-1-p] = tmp_pp;
				tmp_pp = orderedPointsRealCoords[p];   orderedPointsRealCoords[p] = orderedPointsRealCoords[orderedP-1-p];
				orderedPointsRealCoords[orderedP-1-p] = tmp_pp;
				bool tmp_bb = orderedPointsIsIntersection[p];  orderedPointsIsIntersection[p] = orderedPointsIsIntersection[orderedP-1-p];
				orderedPointsIsIntersection[orderedP-1-p] = tmp_bb;
			}
		}
	} else {
		// ordered in Y direction
		if ( orderedPoints[0][1] > orderedPoints[orderedP-1][1]+eps) {
			SUNDANCE_MSG3(verb_, " Order in Y direction ");
			for (int p = 0 ; p < (orderedP/2) ; p++){  // swap the three elements
				Point tmp_pp = orderedPoints[p];    orderedPoints[p] = orderedPoints[orderedP-1-p];
				orderedPoints[orderedP-1-p] = tmp_pp;
				tmp_pp = orderedPointsRealCoords[p];   orderedPointsRealCoords[p] = orderedPointsRealCoords[orderedP-1-p];
				orderedPointsRealCoords[orderedP-1-p] = tmp_pp;
				bool tmp_bb = orderedPointsIsIntersection[p];   orderedPointsIsIntersection[p] = orderedPointsIsIntersection[orderedP-1-p];
				orderedPointsIsIntersection[orderedP-1-p] = tmp_bb;
			}
		}
	}

	// print the points which were ordered
	for (int pi = 0 ; pi < orderedP ; pi++){
		//SUNDANCE_MSG3(verb_, "Ordered Points pi:"<< pi << " P:" << orderedPoints[pi] );
	}


	// get all the possible weights and quadrature points
	Array<Point> linePoints;
	Array<double> lineWeights;
	getLineRule( linePoints , lineWeights);
	int nr1DPoints = linePoints.size();
	int nr2DPoints = nr1DPoints * nr1DPoints;
	Array<Point> quadQuadPoints;
	Array<double> quadQuadWeights;
	getQuadRule( quadQuadPoints , quadQuadWeights);
	Array<Point> trianglePoints;
	Array<double> triangleWeights;
	getTriangleQuadPoints( trianglePoints , triangleWeights);

	// the integration coefficients , we get the coefficient left from the curve and right from the curve
	double alphaLeft , alphaRight;
	Point curvderi = orderedPointsRealCoords[0] - orderedPointsRealCoords[1];
	Point norm_vec = Point( -curvderi[1] , curvderi[0] );
	norm_vec = norm_vec / ::sqrt(norm_vec*norm_vec);
	double relative_length = 1e-2 * ::sqrt(curvderi*curvderi);
	//SUNDANCE_MSG3(verb_, " First segment P0: " << orderedPointsRealCoords[0] << " , P1:" <<orderedPointsRealCoords[1]);
	//SUNDANCE_MSG3(verb_, " Alpha left P: " << (0.5*orderedPointsRealCoords[0] + 0.5*orderedPointsRealCoords[1]) + relative_length*norm_vec );
	//SUNDANCE_MSG3(verb_, " Alpha right P: " << (0.5*orderedPointsRealCoords[0] + 0.5*orderedPointsRealCoords[1]) - relative_length*norm_vec );
	alphaLeft = globalCurve.integrationParameter( (0.5*orderedPointsRealCoords[0] + 0.5*orderedPointsRealCoords[1]) + relative_length*norm_vec );
	alphaRight = globalCurve.integrationParameter( (0.5*orderedPointsRealCoords[0] + 0.5*orderedPointsRealCoords[1]) - relative_length*norm_vec );

	//SUNDANCE_MSG3(verb_, "minX:"<< minX << " , maxX:" << maxX << " , minY:" << minY << " , maxY:" << maxY <<
	//		" , alphaLeft:" << alphaLeft << " , alphaRight:" << alphaRight);
	//SUNDANCE_MSG3(verb_, "increaseX:"<< increaseX << " , increaseY:" << increaseY );

	// first quadrate the whole quad
	Array<double> wholeWeightsQuad( quadWeights.size() , 0.0 );
    makeInterpolantQuad( 0.0, 0.0 , 1.0 , 1.0, nr1DPoints , nr2DPoints ,linePoints ,
    		             quadQuadPoints , quadQuadWeights ,wholeWeightsQuad , 1.0 );

	Array<double> leftWeightsQuad( quadWeights.size() , 0.0 );
	// here we have a big branch, eighter to parse in the X direction or in the Y direction
	if ( (increaseX) || (increaseY == false)){
		// quadrate in the X direction
		// first make a quadrature from 0 to MinX
		SUNDANCE_MSG3(verb_, " Increase X , integrate initial quad");
	    makeInterpolantQuad( 0.0, orderedPoints[0][1] ,
	    		             minX , 1.0-orderedPoints[0][1],
	    		             nr1DPoints , nr2DPoints ,linePoints ,
	    		             quadQuadPoints , quadQuadWeights ,leftWeightsQuad , 1.0 );

	    // for each segment make a triangle and quad quadrature
	    for (int seg = 1 ; seg < orderedP ; seg++){
	    	//SUNDANCE_MSG3(verb_, " Integrate segment : " << seg );
	    	Point p0 = orderedPoints[seg-1];
	    	Point p1 = orderedPoints[seg];
	    	if (p0[0] < p1[0] ){
	    		// X is increasing
	    		if (p0[1] < p1[1]){
	    			// Y increasing
	    		    makeInterpolantQuad( p0[0] , p1[1] ,
	    		    		             p1[0] - p0[0] , 1.0-p1[1] ,
	    		    		             nr1DPoints , nr2DPoints ,linePoints ,
	    		    		             quadQuadPoints , quadQuadWeights ,leftWeightsQuad , 1.0 );
	    		    makeInterpolantQuad( p0[0] , p1[1] ,
	    		                         p1[0] - p0[0] , p0[1] - p1[1] ,
	    		    		             nr1DPoints , nr2DPoints ,linePoints ,
	    		    		             trianglePoints , triangleWeights ,leftWeightsQuad , 2.0  );
	    		}else{
	    			// Y decreasing
	    		    makeInterpolantQuad( p0[0] , p0[1] ,
	    		    		             p1[0] - p0[0] , 1.0-p0[1] ,
	    		    		             nr1DPoints , nr2DPoints ,linePoints ,
	    		    		             quadQuadPoints , quadQuadWeights ,leftWeightsQuad , 1.0 );
	    		    makeInterpolantQuad( p1[0] , p0[1] ,
	    		                         p0[0] - p1[0] , p1[1] - p0[1] ,
	    		    		             nr1DPoints , nr2DPoints ,linePoints ,
	    		    		             trianglePoints , triangleWeights ,leftWeightsQuad , 2.0  );
	    		}
	    	}
	    	else{
	    		// X is decreasing (this should happen rather seldom)
	    		SUNDANCE_MSG3(verb_, " Decrease X ");
	    		if (p0[1] < p1[1]){
	    			// Y increasing
	    		    makeInterpolantQuad( p1[0] , 0.0 ,
	    		    		             p1[0] - p0[0] , p0[1] ,
	    		    		             nr1DPoints , nr2DPoints ,linePoints ,
	    		    		             quadQuadPoints , quadQuadWeights ,leftWeightsQuad , 1.0 );
	    		    makeInterpolantQuad( p1[0] , p0[1] ,
	    		                         p0[0] - p1[0] , p1[1] - p0[1] ,
	    		    		             nr1DPoints , nr2DPoints ,linePoints ,
	    		    		             trianglePoints , triangleWeights ,leftWeightsQuad , 2.0  );
	    		}else{
	    			// Y decreasing
	    		    makeInterpolantQuad( p1[0] , 0.0 ,
	    		    		             p0[0] - p1[0] , p1[1] ,
	    		    		             nr1DPoints , nr2DPoints ,linePoints ,
	    		    		             quadQuadPoints , quadQuadWeights ,leftWeightsQuad , 1.0 );
	    		    makeInterpolantQuad( p0[0] , p1[1] ,
	    		                         p1[0] - p0[0] , p0[1] - p1[1] ,
	    		    		             nr1DPoints , nr2DPoints ,linePoints ,
	    		    		             trianglePoints , triangleWeights ,leftWeightsQuad , 2.0  );
	    		}
	    	}
	    }

	    //SUNDANCE_MSG3(verb_, " Integrate rest of the quad ");
	    // quadrate the rest of the quad
	    makeInterpolantQuad( maxX , orderedPoints[orderedP-1][1] ,
	    		             1.0-maxX , 1.0-orderedPoints[orderedP-1][1] ,
	                         nr1DPoints , nr2DPoints ,linePoints ,
	    		             quadQuadPoints , quadQuadWeights ,leftWeightsQuad , 1.0 );
	}
	else{
		SUNDANCE_MSG3(verb_, " Increase Y , integrate initial quad");
		// quadrate in Y direction
		// first make a quadrature from 0 to MinY
	    makeInterpolantQuad( 0.0 , 0.0 ,
	    		             orderedPoints[0][0] , minY,
	    		             nr1DPoints , nr2DPoints ,linePoints ,
	    		             quadQuadPoints , quadQuadWeights ,leftWeightsQuad , 1.0 );

	    // for each segment make a triangle and quad quadrature
	    for (int seg = 1 ; seg < orderedP ; seg++){
	    	//SUNDANCE_MSG3(verb_, " Integrate segment : " << seg );
	    	Point p0 = orderedPoints[seg-1];
	    	Point p1 = orderedPoints[seg];
	    	// Y must be always increasing !!!
	    	if (p0[0] < p1[0]){
	    		// X increasing
	    	    makeInterpolantQuad( 0.0 , p0[1] ,
	    	    		             p0[0] , p1[1]-p0[1] ,
	    	    		             nr1DPoints , nr2DPoints ,linePoints ,
	    	    		             quadQuadPoints , quadQuadWeights ,leftWeightsQuad , 1.0 );
	    	    makeInterpolantQuad( p0[0] , p1[1] ,
	    	                         p1[0] - p0[0] , p0[1] - p1[1] ,
	    	    		             nr1DPoints , nr2DPoints ,linePoints ,
	    	    		             trianglePoints , triangleWeights ,leftWeightsQuad , 2.0  );
	    	}else{
	    		// X decreasing
	    	    makeInterpolantQuad( 0.0 , p0[1] ,
	    	    		             p1[0] , p1[1]-p0[1] ,
	    	    		             nr1DPoints , nr2DPoints ,linePoints ,
	    	    		             quadQuadPoints , quadQuadWeights ,leftWeightsQuad , 1.0 );
	    	    makeInterpolantQuad( p1[0] , p0[1] ,
	    	                         p0[0] - p1[0] , p1[1] - p0[1] ,
	    	    		             nr1DPoints , nr2DPoints ,linePoints ,
	    	    		             trianglePoints , triangleWeights ,leftWeightsQuad , 2.0  );
	    	}
	    }

	    //SUNDANCE_MSG3(verb_, " Integrate rest of the quad ");
	    // quadrate the rest of the quad
	    makeInterpolantQuad( 0.0 , maxY ,
	    		             orderedPoints[orderedP-1][0] , 1.0-maxY ,
	                         nr1DPoints , nr2DPoints ,linePoints ,
	    		             quadQuadPoints , quadQuadWeights ,leftWeightsQuad , 1.0 );
	}


	// here just add up   W(i) = alphaLeft*leftW(i) + alphaRight*(wholeW(i)-leftW(i))
	double sumWeights = 0.0;
	for (int q = 0 ; q < quadPoints.size() ; q++ ){
		quadWeights[q] = alphaLeft * leftWeightsQuad[q] + alphaRight*( wholeWeightsQuad[q]- leftWeightsQuad[q] );
		sumWeights = sumWeights + quadWeights[q];
	}
	SUNDANCE_MSG3(verb_, " ---------- END METHOD -----------  sum weights = " << sumWeights
			<< "   area:" << sumWeights*(quadOfPoints[3][0]-quadOfPoints[0][0])*(quadOfPoints[3][1]-quadOfPoints[0][1]));
}



void GaussLobattoQuadrature::getAdaptedQuadWeights(int cellLID, const Mesh& mesh,
		const ParametrizedCurve& globalCurve, Array<Point>& quadPoints,
		Array<double>& quadWeights, bool& weightsChanged) const{

	int maxStack = 1000;
	//int maxLevel = 5;
	int stackIndex = 0;
	Tabs tabs;

	Array< Array<Point> > pointStack(maxStack);     // -> the intersection points of each quad cell in the stack (in ref coords)
	Array< Array<Point> > quadPointStack(maxStack); // -> the points which form the quad cell (in ref coords)
	Array< Array<int> > intersectionPointStack(maxStack); // -> which edges are intersected for the actual quad cell
	Array<int> levelStack(maxStack);                      // -> the refinement level of the actual quad cell
	Array<int> refinedStack(maxStack,-1);                 // -> if a quad cell has been replaced then this will be > 0 , -1 otherwise
	Array<int> intersectioncaseStack(maxStack,-1);        // -> the intersection case can be from 0,..,6 , if it is something else then refine
	Array<double> alpha1(maxStack);                       // -> the integration coefficient for the first part of the integral
	Array<double> alpha2(maxStack);                       // -> the integration coefficient for the second part of the integral
	// the actual points of the quad cell and the actual intersection points in real coords
	Array<Point> quadOfPoints(4);
	Array<Point> intesectPoints(8);
	int nrIPoints , tmp;

	// first add the initial quad to the quad stack
	for (int i = 0 ; i < 4 ; i++){
		quadOfPoints[i] = mesh.nodePosition( mesh.facetLID(2,cellLID , 0 , i , nrIPoints ) );
	}

	// first divide the parent cells in smaller quads ,
	pointStack[0].resize(0);
	quadPointStack[0].resize(4);
	quadPointStack[0][0] = Point(0.0,0.0); quadPointStack[0][1] = Point(1.0,0.0);
	quadPointStack[0][2] = Point(0.0,1.0); quadPointStack[0][3] = Point(1.0,1.0);
	intersectionPointStack[0].resize(0);
	levelStack[0] = 0;
	tmp = 0;
	// for each edge test for intersection
	for (int i = 0 ; i < 4 ; i++)
	{
		// for each edge test intersection point
		globalCurve.returnIntersectPoints( quadOfPoints[quadsEdgesPoints[i][0]] , quadOfPoints[quadsEdgesPoints[i][1]],
				 nrIPoints , intesectPoints );
		SUNDANCE_MSG3(verb_, " Test edge: "<< quadOfPoints[quadsEdgesPoints[i][0]] << " -- " << quadOfPoints[quadsEdgesPoints[i][1]] );
		pointStack[0].resize( tmp + nrIPoints );
		intersectionPointStack[0].resize( tmp + nrIPoints );
		for (int j = 0 ; j < nrIPoints ; j++){
			// transform to reference coordinates back
			pointStack[0][tmp+j] = Point( ( intesectPoints[j][0] - quadOfPoints[0][0])/( quadOfPoints[3][0]-quadOfPoints[0][0] ) ,
					                      ( intesectPoints[j][1] - quadOfPoints[0][1])/( quadOfPoints[3][1]-quadOfPoints[0][1] )  );
			intersectionPointStack[0][tmp+j] = i;
			SUNDANCE_MSG3(verb_, " adding Ints Points:"<< tmp+j << " edge:" << intersectionPointStack[0][tmp+j] << " P:"
					<< pointStack[0][tmp+j] << " Real P:" << intesectPoints[j] );
		}
		tmp = tmp + nrIPoints;
	}

	// select the correct case
	intersectioncaseStack[0] = -1;
	if (intersectionPointStack[0].size() == 0){
		intersectioncaseStack[0] = 6; // no intersection point
		alpha1[0] = alpha2[0] = globalCurve.integrationParameter(quadOfPoints[0]);
	}
	if (intersectionPointStack[0].size() == 2){
       if ((intersectionPointStack[0][0] == 0 ) && (intersectionPointStack[0][1] == 1 )) {
    	   intersectioncaseStack[0] = 0;
    	   alpha1[0] = globalCurve.integrationParameter(quadOfPoints[0]);
    	   alpha2[0] = globalCurve.integrationParameter(quadOfPoints[3]);
       }
       if ((intersectionPointStack[0][0] == 0 ) && (intersectionPointStack[0][1] == 2 )) {
    	   intersectioncaseStack[0] = 2;
    	   alpha1[0] = globalCurve.integrationParameter(quadOfPoints[0]);
    	   alpha2[0] = globalCurve.integrationParameter(quadOfPoints[1]);
       }
       if ((intersectionPointStack[0][0] == 0 ) && (intersectionPointStack[0][1] == 3 )) {
    	   alpha1[0] = globalCurve.integrationParameter(quadOfPoints[0]);
    	   alpha2[0] = globalCurve.integrationParameter(quadOfPoints[3]);
    	   if (pointStack[0][0][0] > pointStack[0][1][0]) intersectioncaseStack[0] = 41;
    	   else intersectioncaseStack[0] = 42;
       }
       if ((intersectionPointStack[0][0] == 1 ) && (intersectionPointStack[0][1] == 2 )) {
    	   alpha1[0] = globalCurve.integrationParameter(quadOfPoints[0]);
    	   alpha2[0] = globalCurve.integrationParameter(quadOfPoints[3]);
    	   if (pointStack[0][0][1] > pointStack[0][1][1]) intersectioncaseStack[0] = 51;
    	   else intersectioncaseStack[0] = 52;
       }
       if ((intersectionPointStack[0][0] == 1 ) && (intersectionPointStack[0][1] == 3 )) {
    	   alpha1[0] = globalCurve.integrationParameter(quadOfPoints[0]);
    	   alpha2[0] = globalCurve.integrationParameter(quadOfPoints[2]);
    	   intersectioncaseStack[0] = 1;
       }
       if ((intersectionPointStack[0][0] == 2 ) && (intersectionPointStack[0][1] == 3 )) {
    	   alpha1[0] = globalCurve.integrationParameter(quadOfPoints[0]);
    	   alpha2[0] = globalCurve.integrationParameter(quadOfPoints[3]);
    	   intersectioncaseStack[0] = 3;
       }
	}
	SUNDANCE_MSG3(verb_, " intersectioncaseStack[0]:"<< intersectioncaseStack[0] << " , alpha1[0]:" << alpha1[0]
	                       << " , alpha2[0]" << alpha2[0]);
	stackIndex = 1;
	//quadsEdgesPoints


	// the criterion for division is if the are after the refinement is different significantly
	// then before the division
	// put these quads with the intersection points in one list
	// todo: here we make no while loop since then also the curve integrals should be adaptive


	// get all the possible weights and quadrature points
	Array<Point> linePoints;
	Array<double> lineWeights;
	getLineRule( linePoints , lineWeights);
	Array<Point> quadQuadPoints;
	Array<double> quadQuadWeights;
	getQuadRule( quadQuadPoints , quadQuadWeights);
	Array<Point> trianglePoints;
	Array<double> triangleWeights;
	getTriangleQuadPoints( trianglePoints , triangleWeights);
	double summWeights = 0.0;

	int nr1DPoints = linePoints.size();
	int nr2DPoints = quadQuadPoints.size();

	//those must be equal
	TEUCHOS_TEST_FOR_EXCEPTION( quadPoints.size() != quadQuadPoints.size() , std::runtime_error ,
			"quadPoints.size() != quadQuadPoints.size() , size1:" << quadPoints.size() << " , size2:" << quadQuadPoints.size());
	for (int q = 0; q < nr2DPoints; q++) {
		SUNDANCE_MSG3(verb_, " Quad point quadWeights["<<q<<"]="<<quadWeights[q]);
		summWeights = summWeights + quadWeights[q];
		quadWeights[q] = 0.0;
	}
	SUNDANCE_MSG3(verb_, " Summ old weights = " << summWeights );

    // for all elements in the list make the quadrature
	for (int listI = 0 ; listI < stackIndex ; listI++){
		// look at this quad only when it is not refined
		if (refinedStack[listI] < 0 ){
			SUNDANCE_MSG3(verb_, tabs << "getAdaptedQuadWeights Integrate quad listI:" << listI <<
					",intersectioncaseStack[listI]:" << intersectioncaseStack[listI]);
			// ========== calculate first the integral over the whole quad
			Array<double> tmpWeightsQuad( quadWeights.size() , 0.0 );
			double  ofx , ofy , px , py;
			// integrate the whole quad, which result will be later used
            ofx = (quadPointStack[listI][1][0] - quadPointStack[listI][0][0]);
            ofy = (quadPointStack[listI][2][1] - quadPointStack[listI][0][1]);
            px = quadPointStack[listI][0][0]; py = quadPointStack[listI][0][1];
            // call the function which makes the quadrature
            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
            		             quadQuadPoints , quadQuadWeights ,tmpWeightsQuad , 1.0 );
            SUNDANCE_MSG3(verb_, tabs << " end quadring the whole quad, now quadring the remaining parts " );
			switch (intersectioncaseStack[listI]){
			// =================================================================================================================
			case 0:{
				Array<double> tmpWeightsTriangle( quadWeights.size() , 0.0 );
				// integrate the triangle
                ofx = (pointStack[listI][0][0] - quadPointStack[listI][0][0]);
                ofy = (pointStack[listI][1][1] - quadPointStack[listI][0][1]);
                px = quadPointStack[listI][0][0]; py = quadPointStack[listI][0][1];
                // call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		trianglePoints , triangleWeights ,tmpWeightsTriangle , 2.0 );
				// now add up the two diferent arreas
				for (int i = 0 ; i < nr2DPoints ; i++){
					SUNDANCE_MSG3(verb_, tabs << "i=" << i << ",Wquad[i]:" << tmpWeightsQuad[i] << " , Wtriag[i]:"<< tmpWeightsTriangle[i] );
					quadWeights[i] = alpha1[listI]*tmpWeightsTriangle[i] + alpha2[listI]*( tmpWeightsQuad[i] - tmpWeightsTriangle[i] );
				}
			    break;}
			// =================================================================================================================
			case 1:{
				Array<double> tmpWeightsTriangle( quadWeights.size() , 0.0 );
				// integrate the triangle
                ofx = (pointStack[listI][1][0] - quadPointStack[listI][0][0]);
                ofy = (pointStack[listI][0][1] - quadPointStack[listI][2][1]);
                px = quadPointStack[listI][2][0]; py = quadPointStack[listI][2][1];
                // call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		trianglePoints , triangleWeights ,tmpWeightsTriangle , 2.0 );
				// now add up the two different areas
				for (int i = 0 ; i < nr2DPoints ; i++){
					SUNDANCE_MSG3(verb_, tabs << "i=" << i << ",Wquad[i]:" << tmpWeightsQuad[i] << " , Wtriag[i]:"<< tmpWeightsTriangle[i] );
					quadWeights[i] = alpha2[listI]*tmpWeightsTriangle[i] + alpha1[listI]*( tmpWeightsQuad[i] - tmpWeightsTriangle[i] );
				}
			    break;}
			// =================================================================================================================
			case 2:{
				Array<double> tmpWeightsTriangle( quadWeights.size() , 0.0 );
				// integrate the triangle
                ofx = (pointStack[listI][0][0] - quadPointStack[listI][1][0]);
                ofy = (pointStack[listI][1][1] - quadPointStack[listI][1][1]);
                px = quadPointStack[listI][1][0]; py = quadPointStack[listI][1][1];
                // call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		trianglePoints , triangleWeights ,tmpWeightsTriangle , 2.0 );
				// now add up the two diferent arreas
				for (int i = 0 ; i < nr2DPoints ; i++){
					SUNDANCE_MSG3(verb_, tabs << "i=" << i << ",Wquad[i]:" << tmpWeightsQuad[i] << " , Wtriag[i]:"<< tmpWeightsTriangle[i] );
					quadWeights[i] = alpha2[listI]*tmpWeightsTriangle[i] + alpha1[listI]*( tmpWeightsQuad[i] - tmpWeightsTriangle[i] );
				}
			    break;}
			// =================================================================================================================
			case 3:{
				Array<double> tmpWeightsTriangle( quadWeights.size() , 0.0 );
				// integrate the triangle
                ofx = (pointStack[listI][1][0] - quadPointStack[listI][3][0]);
                ofy = (pointStack[listI][0][1] - quadPointStack[listI][3][1]);
                px = quadPointStack[listI][3][0]; py = quadPointStack[listI][3][1];
                // call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		trianglePoints , triangleWeights ,tmpWeightsTriangle , 2.0 );
				// now add up the two diferent arreas
				for (int i = 0 ; i < nr2DPoints ; i++){
					SUNDANCE_MSG3(verb_, tabs << "i=" << i << ",Wquad[i]:" << tmpWeightsQuad[i] << " , Wtriag[i]:"<< tmpWeightsTriangle[i] );
					quadWeights[i] = alpha2[listI]*tmpWeightsTriangle[i] + alpha1[listI]*( tmpWeightsQuad[i] - tmpWeightsTriangle[i] );
				}
			    break;}
			// =================================================================================================================
			case 41: case 42:{
				// integrate the quad
				Array<double> tmpWeightsQuad2( quadWeights.size() , 0.0 );
				if (intersectioncaseStack[listI] == 41){
					ofx = ( pointStack[listI][1][0] - quadPointStack[listI][0][0]);
				} else {
					ofx = ( pointStack[listI][0][0] - quadPointStack[listI][0][0]);
				}
	            ofy = (quadPointStack[listI][2][1] - quadPointStack[listI][0][1]);
	            px = quadPointStack[listI][0][0]; py = quadPointStack[listI][0][1];
	            // call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		quadQuadPoints , quadQuadWeights ,tmpWeightsQuad2 , 1.0 );
				// integrate the triangle
				Array<double> tmpWeightsTriangle( quadWeights.size() , 0.0 );
				if (intersectioncaseStack[listI] == 41){
					ofx = ( pointStack[listI][0][0] - pointStack[listI][1][0]);
	                ofy = ( quadPointStack[listI][3][1] - quadPointStack[listI][0][1] );
					px = pointStack[listI][1][0]; py = pointStack[listI][0][1];
				} else {
					ofx = ( pointStack[listI][1][0] - pointStack[listI][0][0]);
					ofy = -( quadPointStack[listI][3][1] - quadPointStack[listI][0][1] );
					px = pointStack[listI][0][0]; py = pointStack[listI][1][1];
				}
				// call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		trianglePoints , triangleWeights ,tmpWeightsTriangle , 2.0 );
				// now add up the two different areas
				for (int i = 0 ; i < nr2DPoints ; i++){
					SUNDANCE_MSG3(verb_, tabs << "i=" << i << " , Wquad[i]:" << tmpWeightsQuad[i] <<
							" , Wquad2[i]:" << tmpWeightsQuad2[i] << " , Wtriag[i]:"<< tmpWeightsTriangle[i] );
					quadWeights[i] = alpha1[listI]*( tmpWeightsQuad2[i] + tmpWeightsTriangle[i] ) +
							         alpha2[listI]*( tmpWeightsQuad[i] - tmpWeightsQuad2[i] - tmpWeightsTriangle[i]);
				}
			    break;}
			// =================================================================================================================
			case 51: case 52:{
				// integrate the quad
				Array<double> tmpWeightsQuad2( quadWeights.size() , 0.0 );
				if (intersectioncaseStack[listI] == 52){
					ofy = ( pointStack[listI][0][1] - quadPointStack[listI][0][1]);
				} else {
					ofy = ( pointStack[listI][1][1] - quadPointStack[listI][0][1]);
				}
	            ofx = (quadPointStack[listI][1][0] - quadPointStack[listI][0][0]);
	            px = quadPointStack[listI][0][0]; py = quadPointStack[listI][0][1];
	            // call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		quadQuadPoints , quadQuadWeights ,tmpWeightsQuad2 , 1.0 );
				// integrate the triangle
				Array<double> tmpWeightsTriangle( quadWeights.size() , 0.0 );
				if (intersectioncaseStack[listI] == 51){
	                ofx = ( quadPointStack[listI][3][0] - quadPointStack[listI][0][0] );
					ofy = ( pointStack[listI][0][1] - pointStack[listI][1][1]);
					px = pointStack[listI][0][0]; py = pointStack[listI][1][1];
				} else {
					ofx = -( quadPointStack[listI][3][0] - quadPointStack[listI][0][0] );
					ofy = ( pointStack[listI][1][1] - pointStack[listI][0][1]);
					px = pointStack[listI][1][0]; py = pointStack[listI][0][1];
				}
				// call the function which makes the quadrature
	            makeInterpolantQuad( px, py, ofx, ofy, nr1DPoints , nr2DPoints ,linePoints ,
	            		trianglePoints , triangleWeights ,tmpWeightsTriangle , 2.0 );
				// now add up the two different areas
				for (int i = 0 ; i < nr2DPoints ; i++){
					SUNDANCE_MSG3(verb_, tabs << "i=" << i << " , Wquad[i]:" << tmpWeightsQuad[i] <<
							" , Wquad2[i]:" << tmpWeightsQuad2[i] << " , Wtriag[i]:"<< tmpWeightsTriangle[i] );
					quadWeights[i] = alpha1[listI]*(tmpWeightsQuad2[i]+tmpWeightsTriangle[i]) +
							         alpha2[listI]*(tmpWeightsQuad[i] - tmpWeightsQuad2[i] - tmpWeightsTriangle[i]);
				}
			    break;}
			// =================================================================================================================
			case 6:{
				// no intersection point just integrate the quad elem
				for (int i = 0 ; i < nr2DPoints ; i++){
				   SUNDANCE_MSG3(verb_, tabs << "i=" << i << " , Wquad[i]:" << tmpWeightsQuad[i] );
                   quadWeights[i] = quadWeights[i] + alpha1[listI]*tmpWeightsQuad[i];
				}
			    break;}
			default:{
				// throw error
				TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error  , "Quad cell not integrable:" << intersectioncaseStack[listI]);
			    break; }
			}
		}
	}
	// just print the weights
	summWeights = 0.0;
	for (int q = 0; q < nr2DPoints; q++) {
		summWeights = summWeights + quadWeights[q];
		SUNDANCE_MSG3(verb_, " New weights quadWeights["<<q<<"]="<<quadWeights[q]);
	}
	SUNDANCE_MSG3(verb_, " Summ new weights = " << summWeights );
}


// =========================================================================================================================

void GaussLobattoQuadrature::getAdaptedQuadWeights_surf(
		    int cellLID, const Mesh& mesh,
			const ParametrizedCurve& globalCurve, Array<Point>& quadPoints,
			Array<double>& quadWeights, bool& weightsChanged) const {


	int verb = 1;
	int nr_point = mesh.numFacets( mesh.spatialDim() , 0 , 0);
	Array<Point> maxCellsPoints(nr_point);
	int tmp_i , point_LID;
	weightsChanged = true;

	// - get all the brick cells point
	SUNDANCE_MSG3(verb, "CurveIntegralCalc::getCurveQuadPoints nr points per cell: " << nr_point)
	for (int jj = 0 ; jj < nr_point ; jj++){
		point_LID = mesh.facetLID( mesh.spatialDim() , cellLID , 0 , jj , tmp_i );
		maxCellsPoints[jj] = mesh.nodePosition(point_LID);
		SUNDANCE_MSG3(verb, " max cell point p[" << jj << "]:"<< maxCellsPoints[jj]);
	}

	// - use the existing method to get the triangles
	Array<Point> intersectPoints;
	Array<int> edgeIndex;
	Array<int> triangleIndex;
    // the results in this form will be in physical coordinates
	SundanceSurf3DCalc::getCurveQuadPoints( BrickCell , cellLID , mesh , globalCurve, maxCellsPoints,
						 intersectPoints ,  edgeIndex , triangleIndex );

	// initialize the output variables
	getBrickRule(quadPoints,quadWeights);

	// if the cell has no intersection points then just return the standard weights
	if (intersectPoints.size() < 1){
		double intCoeff = globalCurve.integrationParameter( 0.5*(maxCellsPoints[0]+maxCellsPoints[7]) );
		SUNDANCE_MSG1(verb, " WARNING cel has not been intersected , alpha = " << intCoeff );
		for (int ii = 0 ; ii < quadWeights.size() ; ii++ ) { quadWeights[ii] = quadWeights[ii] * intCoeff;}
		return;
	}

	// - determine the projection direction , and the points which overlap with the projected quad cell
	int nrIntPoint = intersectPoints.size() , projDir = -1 , tmpI , coveredNodes = 0 , firstCoveredNode = -1;
	int possibleProject[3] = { -1 , -1 , -1};
	int otherDims[2] = {-1,-1};
	int projectedEdgeCovered[4] = { -1 , -1 , -1 , -1 };
	double dist_tmp = 1e+100;

	// find out the possible projections (there is only three possible projections)
	for (int ii = 0 ; ii < nrIntPoint ; ii++){
		// for each edge index we set the projection directory to one
		possibleProject[ edge3DProjection[ edgeIndex[ii] ] ] = 1;
	}
	// here we find the projection directory
	for (int ii = 0 ; ii < 3 ; ii++){
		if (possibleProject[ii] > 0) {
			SUNDANCE_MSG3(verb, " Projected dir = " << ii);
			projDir = ii; break;
		}
	}
	// here we store the two other directions which are not projected
	tmpI = 0;
	for (int ii = 0 ; ii < 3 ; ii++){
		if ( ii != projDir ) {
			SUNDANCE_MSG3(verb, " otherDims["<< tmpI << "] = " << ii);
			otherDims[tmpI] = ii; tmpI++;
		}
	}

	// - transform all intersection points to reference coordinates
	for (int ii = 0 ; ii < nrIntPoint ; ii++){
		intersectPoints[ii] = Point( (intersectPoints[ii][0] - maxCellsPoints[0][0])/(maxCellsPoints[7][0] - maxCellsPoints[0][0]) ,
				 (intersectPoints[ii][1] - maxCellsPoints[0][1])/(maxCellsPoints[7][1] - maxCellsPoints[0][1]) ,
				 (intersectPoints[ii][2] - maxCellsPoints[0][2])/(maxCellsPoints[7][2] - maxCellsPoints[0][2]));
		SUNDANCE_MSG3(verb, "REF intersectPoints[" << ii << "] = "<< intersectPoints[ii]);
	}

	// - the quad points on the projected
	Array<Point> quadProjectPoints(4);
	quadProjectPoints[0] = Point(0.0,0.0); quadProjectPoints[1] = Point(1.0,0.0);
	quadProjectPoints[2] = Point(0.0,1.0); quadProjectPoints[3] = Point(1.0,1.0);

	// - determine which out of the 4 nodes has been covered
	coveredNodes = 0;
	for (int ii = 0 ; ii < nrIntPoint ; ii++) {
		Point tmpP( intersectPoints[ii][otherDims[0]] , intersectPoints[ii][otherDims[1]] );
		for (int tt = 0 ; tt < 4 ; tt++ ) {
			if ( ((tmpP-quadProjectPoints[tt])*(tmpP-quadProjectPoints[tt]) < 1e-12) && (projectedEdgeCovered[tt] < 0) )
			{
				projectedEdgeCovered[tt] = 1; // this node is covered
				//SUNDANCE_MSG3(verb, "COVERED node[" << tt << "] = "<< quadProjectPoints[tt]);
				coveredNodes++;
				// to store the first node which is covered is important to measure the integration coefficient
				if (firstCoveredNode < 0) { firstCoveredNode = tt; }
			}
		}
	}
	SUNDANCE_MSG3( verb , "covered nodes = "<< coveredNodes );

	// if not all nodes in the projected quad are covered then we have to add them
	// so we loop as long this is done
	while (coveredNodes < 4) {
	    // - look for one edge, which has at least one neighbor which is covered
		int foundEdgeIndex[3] = {-1,-1,-1} , foundI = 0 , pF = -1 , pL = -1;
		int neighbourNodes[4][2] = { {1,2} , {0,3} , {0,3} , {1,2} };
		for (int ii = 0 ; ii < 4 ; ii++){
			//SUNDANCE_MSG3( verb , "test node = "<< ii << " , projectedEdgeCovered[ii]=" << projectedEdgeCovered[ii] <<
			//		" , N1=" << projectedEdgeCovered[ neighbourNodes[ii][0] ] << " , N2=" << projectedEdgeCovered[ neighbourNodes[ii][1] ]);
			// look for uncovered nodes, which have one neighbor which is covered
			if ( (projectedEdgeCovered[ii] < 0) &&
				 ((projectedEdgeCovered[ neighbourNodes[ii][0] ] > 0) || (projectedEdgeCovered[ neighbourNodes[ii][1] ] > 0)) ){
				//SUNDANCE_MSG3( verb , "found node = "<< ii );
				foundEdgeIndex[foundI] = ii;
				projectedEdgeCovered[ii] = 1;
				coveredNodes++;
				break;
			}
		}
		foundI = 0;
		// look for the nearest point, with respect to the found node
		dist_tmp = 1e+100;
		for (int ii = 0 ; ii < nrIntPoint ; ii++){
			Point tmpP( intersectPoints[ii][otherDims[0]] , intersectPoints[ii][otherDims[1]] );
			double dist_tmp_tmp = ::sqrt( (tmpP-quadProjectPoints[foundEdgeIndex[foundI]])*(tmpP-quadProjectPoints[foundEdgeIndex[foundI]]) );
			//SUNDANCE_MSG3( verb , "test projected Int Point = "<< tmpP << " , dist = " << dist_tmp_tmp);
			if ( dist_tmp_tmp < dist_tmp) { pF = ii; dist_tmp = dist_tmp_tmp;}
		}
		SUNDANCE_MSG3( verb , "pF = "<< pF );

		// loop as long all the neighbors of the neighbors are covered
		while ( (projectedEdgeCovered[ neighbourNodes[foundEdgeIndex[foundI]][0] ] < 0) ||
		      (projectedEdgeCovered[ neighbourNodes[foundEdgeIndex[foundI]][1] ] < 0) ) {
			int tmp_next_Neighbor = -1;
			// only one can be true
			if (projectedEdgeCovered[ neighbourNodes[foundEdgeIndex[foundI]][0] ] < 0) { tmp_next_Neighbor = neighbourNodes[foundEdgeIndex[foundI]][0]; }
			if (projectedEdgeCovered[ neighbourNodes[foundEdgeIndex[foundI]][1] ] < 0) { tmp_next_Neighbor = neighbourNodes[foundEdgeIndex[foundI]][1]; }
			//SUNDANCE_MSG3( verb , "tmp_next_Neighbor = "<< tmp_next_Neighbor );
			foundEdgeIndex[foundI+1] = tmp_next_Neighbor;
			projectedEdgeCovered[tmp_next_Neighbor] = 1;
			foundI++;
			coveredNodes++;
		}
		// once we finish the loop we look for a next intersection point which is closer to the last node and is different from "pF"
		dist_tmp = 1e+100;
		for (int ii = 0 ; ii < nrIntPoint ; ii++){
			Point tmpP( intersectPoints[ii][otherDims[0]] , intersectPoints[ii][otherDims[1]] );
			double dist_tmp_tmp = ::sqrt((tmpP-quadProjectPoints[foundEdgeIndex[foundI]])*(tmpP-quadProjectPoints[foundEdgeIndex[foundI]]));
			if ( (dist_tmp_tmp < dist_tmp) && (ii != pF ) ) { pL = ii; dist_tmp = dist_tmp_tmp; }
		}
		SUNDANCE_MSG3( verb , "pL = "<< pF );

		// - now we have a list of the nodes of the projected quad, which need to be added
		for (int nn = 0 ; nn <= foundI ; nn = nn + 1)
		{
			int firstI_Tri = pF, secondI_Tri = -1 , thirdI_Tri = -1;
			intersectPoints.resize( nrIntPoint + 1);
			nrIntPoint++;
			//  check if one point will not be added twice ... on the other side this is not a problem for the integration
			Point ptmp( 0 , 0 , 0 );
			ptmp[otherDims[0]] = quadProjectPoints[foundEdgeIndex[nn]][0];
			ptmp[otherDims[1]] = quadProjectPoints[foundEdgeIndex[nn]][1];
			ptmp[projDir] = intersectPoints[pF][projDir];
			//SUNDANCE_MSG3( verb , "Add Point 2D " << quadProjectPoints[foundEdgeIndex[nn]] << " , P:" << ptmp);
			intersectPoints[nrIntPoint-1] = ptmp;
			secondI_Tri = nrIntPoint -1;
			Point nextP( 0 , 0 , 0 );
			if (nn < foundI ){
				nextP[otherDims[0]] = quadProjectPoints[foundEdgeIndex[nn+1]][0];
				nextP[otherDims[1]] = quadProjectPoints[foundEdgeIndex[nn+1]][1];
				nextP[projDir] = intersectPoints[pF][projDir];
				//SUNDANCE_MSG3( verb , "Add Point 2D " << quadProjectPoints[foundEdgeIndex[nn+1]] );
				intersectPoints.resize( nrIntPoint + 1);
				nrIntPoint++;
				thirdI_Tri = nrIntPoint - 1;
				intersectPoints[nrIntPoint-1] = nextP;
			} else {
				// the last point is "pL"
				nextP = intersectPoints[pL];
				thirdI_Tri = pL;
			}
			int triagSize = triangleIndex.size();
			triangleIndex.resize(triagSize + 3);
			triangleIndex[triagSize] = firstI_Tri;
			triangleIndex[triagSize+1] = secondI_Tri;
			triangleIndex[triagSize+2] = thirdI_Tri;
			SUNDANCE_MSG3( verb , "Add triangle [ "<< firstI_Tri << " , " <<  secondI_Tri << " , " << thirdI_Tri << " ]");
			SUNDANCE_MSG3( verb , "Added triangle's points [ "<< intersectPoints[firstI_Tri] << " , "
					<<  intersectPoints[secondI_Tri] << " , " << intersectPoints[thirdI_Tri] << " ]");
		}

		SUNDANCE_MSG3( verb , "covered nodes = "<< coveredNodes );
	} //end from the while which makes sure that each nodes are covered


	// - first get the quadrature points for the whole quad
	Array<Point> wholeQuadPoints;
	Array<double> wholeQuadweights;
	getBrickRule(wholeQuadPoints,wholeQuadweights);
	Array<Point> linePoints;
	Array<double> lineWeights;
	getLineRule(linePoints,lineWeights);
	int nrLinePoints = linePoints.size() , nr3DPoints = wholeQuadPoints.size();
	// this will be the output of the quadrature
	Array<double> computedWeights(nr3DPoints,0.0);

	// compute the quadrature points and weights for the prism elements
	Array<Point> triagPoints;
	Array<double> triagWeights;
	getTriangleQuadPoints( triagPoints , triagWeights );
	int nrQuadTriagPoints = triagPoints.size();
	int nrQuadPrism = triagPoints.size() * linePoints.size();
	Array<Point> prismPoints(nrQuadPrism);
	Array<double> prismWeights(nrQuadPrism);

	// - integrate all the prism elements, if the heights are different then zero
	int nr_prism = triangleIndex.size()/3;
	for (int p = 0 ; p < nr_prism ; p++ ){
		Point p0 = intersectPoints[ triangleIndex[3*p] ];
		Point p1 = intersectPoints[ triangleIndex[3*p + 1] ];
		Point p2 = intersectPoints[ triangleIndex[3*p + 2] ];
		Point p0L = p0;
		Point p1L = p1;
		Point p2L = p2;
		p0L[projDir] = 0.0; p1L[projDir] = 0.0; p2L[projDir] = 0.0;
		//SUNDANCE_MSG3( verb , " Triangle [ " << p0 << " , "  << p1 << " , " << p2 << "]");
		//SUNDANCE_MSG3( verb , " Triangle bottom [ " << p0L << " , "  << p1L << " , " << p2L << "]");

		// only integrate the prism when there is a measurable height
		if ( (::fabs(p0[projDir]-p0L[projDir]) + ::fabs(p1[projDir]-p1L[projDir])+ ::fabs(p2[projDir]-p2L[projDir]) ) > 1e-12 ){
			Point normAreaVect = cross(p1L-p0L ,p2L-p0L);
			Point triagPointTmp;
			double areaTriag = 0.5*::sqrt(normAreaVect*normAreaVect) , Zc = 0.0;
			int quadPIndex = 0;
			// compute the quadrature points for this prism
			for (int triagQIndex = 0 ; triagQIndex < nrQuadTriagPoints ; triagQIndex++)
			{
				// Zc is the height of the prism at the quadrature position
			    Zc = (1.0 - triagPoints[triagQIndex][0] - triagPoints[triagQIndex][1]) * ( p0[projDir] - p0L[projDir] );
				Zc = Zc + triagPoints[triagQIndex][0]*( p1[projDir] - p1L[projDir] );
				Zc = Zc + triagPoints[triagQIndex][1]*( p2[projDir] - p2L[projDir] );
				triagPointTmp = p0L + triagPoints[triagQIndex][0]*(p1L-p0L) + triagPoints[triagQIndex][1]*(p2L-p0L);
				//SUNDANCE_MSG3( verb , " areaTriag = " << areaTriag << " , Zc = "  << Zc <<
				//		" , triagWeights[triagQIndex] = " << triagWeights[triagQIndex] << " , lineWeights[0]" << lineWeights[0]);
				for (int lineQIndex = 0 ; lineQIndex < nrLinePoints ; lineQIndex++ , quadPIndex++)
				{
					// take the tensor product of the triangle and the line quadrature points
					// it is important that we compute correctly the weights and the quadrature points positions
					prismPoints[quadPIndex] = triagPointTmp;
					prismPoints[quadPIndex][projDir] = Zc * linePoints[lineQIndex][0];
					prismWeights[quadPIndex] = Zc*areaTriag*triagWeights[triagQIndex]*lineWeights[lineQIndex];
				}
			}
			// call the method to integrate the Lagrange polynoms
			makeInterpolantBrick( nrLinePoints , nr3DPoints , linePoints ,
					prismPoints , prismWeights , computedWeights , 1.0);
		}
	}

	// - determine the integration coefficient , use "firstCoveredNode" with "0" and "1" offset
	Point tmpPointIntCoeff = intersectPoints[firstCoveredNode];
	tmpPointIntCoeff[projDir] = -0.001;
	Point tmpPointIntCoeff1 = Point( maxCellsPoints[0][0] + tmpPointIntCoeff[0]*(maxCellsPoints[7][0]-maxCellsPoints[0][0]) ,
			                         maxCellsPoints[0][1] + tmpPointIntCoeff[1]*(maxCellsPoints[7][1]-maxCellsPoints[0][1]) ,
			                         maxCellsPoints[0][2] + tmpPointIntCoeff[2]*(maxCellsPoints[7][2]-maxCellsPoints[0][2]) );
	double zeroIntCoeff = globalCurve.integrationParameter( tmpPointIntCoeff1 );
	tmpPointIntCoeff[projDir] = 1.001;
	tmpPointIntCoeff1 = Point( maxCellsPoints[0][0] + tmpPointIntCoeff[0]*(maxCellsPoints[7][0]-maxCellsPoints[0][0]) ,
				               maxCellsPoints[0][1] + tmpPointIntCoeff[1]*(maxCellsPoints[7][1]-maxCellsPoints[0][1]) ,
				               maxCellsPoints[0][2] + tmpPointIntCoeff[2]*(maxCellsPoints[7][2]-maxCellsPoints[0][2]) );
	double oneIntCoeff = globalCurve.integrationParameter( tmpPointIntCoeff1 );

	SUNDANCE_MSG3( verb , " alpha1 = " << zeroIntCoeff << " , alpha2 = "  << oneIntCoeff );

	// at the end just compute the weights, by simply taking the difference of the whole and the computed weights
	double sumWeights = 0.0;
	for (int q = 0 ; q < computedWeights.size() ; q++ ){
		quadWeights[q] = zeroIntCoeff * computedWeights[q] + oneIntCoeff*( wholeQuadweights[q]- computedWeights[q]);
		sumWeights = sumWeights + quadWeights[q];
	}

	SUNDANCE_MSG3(verb , " ---------- END METHOD -----------  sum weights = " << sumWeights);

}


// ======================================================================================================




void GaussLobattoQuadrature::getLineRule(
		Array<Point>& quadPointsL,
		Array<double>& quadWeights) const {

	int n = order()+1; //this is correct m points for a basis of order m-1

	Array<double> quadPoints;
	quadPoints.resize(n);
	quadPointsL.resize(n);
	quadWeights.resize(n);

	// ================== QUAD POINTS ========================
	switch (n) {
	case 2: { quadPoints[0]=0.0; quadPoints[1] = 1; break; }
	case 3: { quadPoints[0]=0.0; quadPoints[1] = 0.5; quadPoints[2] = 1.0; break; }
	case 4:	{ quadPoints[0]=0.0; quadPoints[1] = 0.5-0.5/::sqrt(5.0); quadPoints[2] = 0.5+0.5/::sqrt(5.0); quadPoints[3] = 1.0; break; }
	case 5:	{ quadPoints[0]=0.0; quadPoints[1] = 0.5-0.5*::sqrt(3.0/7.0); quadPoints[2] = 0.5;
	        quadPoints[3] = 0.5+0.5*::sqrt(3.0/7.0); quadPoints[4] = 1.0; break; }
	case 6: { double t0=::sqrt(7.0);
	        double t1=(7.0+2.0*t0)/21.0;
	        double t2=(7.0-2.0*t0)/21.0;
	        quadPoints[0] = 0; quadPoints[1] = 0.5-0.5*::sqrt(t1); quadPoints[2] = 0.5-0.5*::sqrt(t2);
	        quadPoints[3] = 0.5+0.5*::sqrt(t2); quadPoints[4] = 0.5+0.5*::sqrt(t1); quadPoints[5] = 1.0; break; }
	case 7: {
		quadPoints[0] = 0.00000000000000000000e+00;
		quadPoints[1] = 8.48880518607165179823e-02;
		quadPoints[2] = 2.65575603264642912116e-01;
		quadPoints[3] = 5.00000000000000000000e-01;
		quadPoints[4] = 7.34424396735357087884e-01;
		quadPoints[5] = 9.15111948139283537529e-01;
		quadPoints[6] = 1.00000000000000000000e+00;
	    break; }
	case 8: {
		quadPoints[0] = 0.00000000000000000000e+00;
		quadPoints[1] = 6.41299257451967141819e-02;
		quadPoints[2] = 2.04149909283428854234e-01;
		quadPoints[3] = 3.95350391048760574364e-01;
		quadPoints[4] = 6.04649608951239425636e-01;
		quadPoints[5] = 7.95850090716571090255e-01;
		quadPoints[6] = 9.35870074254803285818e-01;
		quadPoints[7] = 1.00000000000000000000e+00;
		break; }
	case 9: {
		quadPoints[0] = 0.00000000000000000000e+00;
		quadPoints[1] = 5.01210022942699118254e-02;
		quadPoints[2] = 1.61406860244631134016e-01;
		quadPoints[3] = 3.18441268086910922452e-01;
		quadPoints[4] = 5.00000000000000000000e-01;
		quadPoints[5] = 6.81558731913089133059e-01;
		quadPoints[6] = 8.38593139755368865984e-01;
		quadPoints[7] = 9.49878997705730032663e-01;
		quadPoints[8] = 1.00000000000000000000e+00;
	    break; }
	case 10: {
		quadPoints[0] = 0.00000000000000000000e+00;
		quadPoints[1] = 4.02330459167705711820e-02;
		quadPoints[2] = 1.30613067447247432895e-01;
		quadPoints[3] = 2.61037525094777733692e-01;
		quadPoints[4] = 4.17360521166806497373e-01;
		quadPoints[5] = 5.82639478833193447116e-01;
		quadPoints[6] = 7.38962474905222266308e-01;
		quadPoints[7] = 8.69386932552752567105e-01;
		quadPoints[8] = 9.59766954083229428818e-01;
		quadPoints[9] = 1.00000000000000000000e+00;
		break; }
	 case 11: {
		quadPoints[0] = 0.00000000000000000000e+00;
		quadPoints[1] = 3.29992847959704183047e-02;
		quadPoints[2] = 1.07758263168427792511e-01;
		quadPoints[3] = 2.17382336501897477365e-01;
		quadPoints[4] = 3.52120932206530290465e-01;
		quadPoints[5] = 5.00000000000000000000e-01;
		quadPoints[6] = 6.47879067793469709535e-01;
		quadPoints[7] = 7.82617663498102578146e-01;
		quadPoints[8] = 8.92241736831572263000e-01;
		quadPoints[9] = 9.67000715204029637206e-01;
		quadPoints[10] = 1.00000000000000000000e+00;
	    break; }
	 case 12: {
		quadPoints[0] = 0.00000000000000000000e+00;
		quadPoints[1] = 2.75503638885589152707e-02;
		quadPoints[2] = 9.03603391779966846897e-02;
		quadPoints[3] = 1.83561923484069688950e-01;
		quadPoints[4] = 3.00234529517325543502e-01;
		quadPoints[5] = 4.31723533572536233294e-01;
		quadPoints[6] = 5.68276466427463766706e-01;
		quadPoints[7] = 6.99765470482674456498e-01;
		quadPoints[8] = 8.16438076515930255539e-01;
		quadPoints[9] = 9.09639660822003315310e-01;
		quadPoints[10] = 9.72449636111441084729e-01;
		quadPoints[11] = 1.00000000000000000000e+00;
	    break; }
	}

	// transform the array of doubles into an array of Points
	for (int i = 0 ; i < n ; i++){
		quadPointsL[i] = Point(quadPoints[i]);
	}

	// ================== WEIGHTS ========================
	switch (n) {
	case 2:{
		quadWeights[0] = 0.5;
		quadWeights[1] = 0.5;
	    break; }
	case 3:{
		quadWeights[0] = 1.0/6.0; quadWeights[1] = 2.0/3.0; quadWeights[2] = 1.0/6.0;
		break;}
	case 4:{
		quadWeights[0] = 1.0/12.0; quadWeights[1] = 5.0/12.0; quadWeights[2] = 5.0/12.0; quadWeights[3] = 1.0/12.0;
		break;}
	case 5:{
		quadWeights[0] = 0.05; quadWeights[1] = 49.0/180.0; quadWeights[2] = 32.0/90.0; quadWeights[3] = 49.0/180.0; quadWeights[4] = 0.05;
		break;}
	case 6:{
	    double t0=::sqrt(7.0);
	    double t1=(7.0+2.0*t0)/21.0;
	    double t2=(7.0-2.0*t0)/21.0;
	    double k1=(1.0-t0)*(1.0-t0);
	    double k2=(1.0+t0)*(1.0+t0);
	    quadWeights[0] = 1.0/30.0; quadWeights[1] = 0.3/(t1*k1); quadWeights[2] = 0.3/(t2*k2); quadWeights[3] = 0.3/(t2*k2);
	    quadWeights[4] = 0.3/(t1*k1); quadWeights[5] = 1.0/30.0;
	    break;}
	case 7:{
		quadWeights[0] = 2.38095238095238082021e-02;
		quadWeights[1] = 1.38413023680782953928e-01;
		quadWeights[2] = 2.15872690604931305458e-01;
		quadWeights[3] = 2.43809523809523809312e-01;
		quadWeights[4] = 2.15872690604931305458e-01;
		quadWeights[5] = 1.38413023680782953928e-01;
		quadWeights[6] = 2.38095238095238082021e-02;
	    break;}
	case 8:{
		quadWeights[0] = 1.78571428571428561516e-02;
		quadWeights[1] = 1.05352113571753072674e-01;
		quadWeights[2] = 1.70561346241752204156e-01;
		quadWeights[3] = 2.06229397329351860080e-01;
		quadWeights[4] = 2.06229397329351860080e-01;
		quadWeights[5] = 1.70561346241752204156e-01;
		quadWeights[6] = 1.05352113571753072674e-01;
		quadWeights[7] = 1.78571428571428561516e-02;
	    break;}
	case 9:{
		quadWeights[0] = 1.38888888888888881179e-02;
		quadWeights[1] = 8.27476807804027880699e-02;
		quadWeights[2] = 1.37269356250080826198e-01;
		quadWeights[3] = 1.73214255486523083238e-01;
		quadWeights[4] = 1.85759637188208620584e-01;
		quadWeights[5] = 1.73214255486523083238e-01;
		quadWeights[6] = 1.37269356250080826198e-01;
		quadWeights[7] = 8.27476807804027880699e-02;
		quadWeights[8] = 1.38888888888888881179e-02;
	    break;}
	case 10:{
		quadWeights[0] = 1.11111111111111115352e-02;
		quadWeights[1] = 6.66529954255350304271e-02;
		quadWeights[2] = 1.12444671031563220298e-01;
		quadWeights[3] = 1.46021341839841889421e-01;
		quadWeights[4] = 1.63769880591948718829e-01;
		quadWeights[5] = 1.63769880591948718829e-01;
		quadWeights[6] = 1.46021341839841889421e-01;
		quadWeights[7] = 1.12444671031563220298e-01;
		quadWeights[8] = 6.66529954255350304271e-02;
		quadWeights[9] = 1.11111111111111115352e-02;
	    break;}
	case 11:{
		quadWeights[0] = 9.09090909090909046752e-03;
		quadWeights[1] = 5.48061366334974126024e-02;
		quadWeights[2] = 9.35849408901525958715e-02;
		quadWeights[3] = 1.24024052132014145355e-01;
		quadWeights[4] = 1.43439562389503921791e-01;
		quadWeights[5] = 1.50108797727845355574e-01;
		quadWeights[6] = 1.43439562389503921791e-01;
		quadWeights[7] = 1.24024052132014145355e-01;
		quadWeights[8] = 9.35849408901525958715e-02;
		quadWeights[9] = 5.48061366334974126024e-02;
		quadWeights[10] = 9.09090909090909046752e-03;
	    break;}
	case 12:{
		quadWeights[0] = 7.57575757575757596785e-03;
		quadWeights[1] = 4.58422587065981240739e-02;
		quadWeights[2] = 7.89873527821850218711e-02;
		quadWeights[3] = 1.06254208880510653268e-01;
		quadWeights[4] = 1.25637801599600640312e-01;
		quadWeights[5] = 1.35702620455348088591e-01;
		quadWeights[6] = 1.35702620455348088591e-01;
		quadWeights[7] = 1.25637801599600640312e-01;
		quadWeights[8] = 1.06254208880510653268e-01;
		quadWeights[9] = 7.89873527821850218711e-02;
		quadWeights[10] = 4.58422587065981240739e-02;
		quadWeights[11] = 7.57575757575757596785e-03;
	    break;}
   }
}


void GaussLobattoQuadrature::getTriangleQuadPoints(Array<Point>& pnt  ,Array<double>& weight ) const{
	// we took this directly from the Matlab code of Prof.Ulbrich
	int deg = nrPointin1D_ + nrPointin1D_ - 2;
	if (deg==2){
	  pnt.resize(3); weight.resize(3);
	  pnt[0] = Point(0.50000000000000000000 , 0.00000000000000000000);
	  pnt[1] = Point(0.50000000000000000000 , 0.50000000000000000000);
	  pnt[2] = Point(0.00000000000000000000 , 0.50000000000000000000);
	  weight[0] = 0.33333333333333333333;
	  weight[1] = 0.33333333333333333333;
	  weight[2] = 0.33333333333333333333;
	} else if (deg==3){
		pnt.resize(4); weight.resize(4);
		pnt[0] = Point(0.33333333333333333333 , 0.33333333333333333333);
		pnt[1] = Point(0.60000000000000000000 , 0.20000000000000000000);
		pnt[2] = Point(0.20000000000000000000 , 0.60000000000000000000);
		pnt[3] = Point(0.20000000000000000000 , 0.20000000000000000000);
		weight[0] = -0.56250000000000000000;
		weight[1] = 0.52083333333333333333;
		weight[2] = 0.52083333333333333333;
		weight[3] = 0.52083333333333333333;
	} else if (deg==4) {
		pnt.resize(6); weight.resize(6);
		pnt[0] = Point(0.816847572980459 , 0.091576213509771);
		pnt[1] = Point(0.091576213509771 , 0.816847572980459);
		pnt[2] = Point(0.091576213509771 , 0.091576213509771);
		pnt[3] = Point(0.108103018168070 , 0.445948490915965);
		pnt[4] = Point(0.445948490915965 , 0.108103018168070);
		pnt[5] = Point(0.445948490915965 , 0.445948490915965);
		weight[0] = 0.109951743655322;
		weight[1] = 0.109951743655322;
		weight[2] = 0.109951743655322;
		weight[3] = 0.223381589678011;
		weight[4] = 0.223381589678011;
		weight[5] = 0.223381589678011;
	} else if (deg==5) {
		pnt.resize(7); weight.resize(7);
		pnt[0] = Point(0.33333333333333333 , 0.33333333333333333);
		pnt[1] = Point(0.79742698535308720 , 0.10128650732345633);
		pnt[2] = Point(0.10128650732345633 , 0.79742698535308720);
		pnt[3] = Point(0.10128650732345633 , 0.10128650732345633);
		pnt[4] = Point(0.05971587178976981 , 0.47014206410511505);
		pnt[5] = Point(0.47014206410511505 , 0.05971587178976981);
		pnt[6] = Point(0.47014206410511505 , 0.47014206410511505);
		weight[0] = 0.22500000000000000;
		weight[1] = 0.12593918054482717;
		weight[2] = 0.12593918054482717;
		weight[3] = 0.12593918054482717;
		weight[4] = 0.13239415278850616;
		weight[5] = 0.13239415278850616;
		weight[6] = 0.13239415278850616;
	} else if (deg==6) {
		pnt.resize(9); weight.resize(9);
		pnt[0] = Point(0.124949503233232 , 0.437525248383384);
		pnt[1] = Point(0.437525248383384 , 0.124949503233232);
		pnt[2] = Point(0.437525248383384 , 0.437525248383384);
		pnt[3] = Point(0.797112651860071 , 0.165409927389841);
		pnt[4] = Point(0.797112651860071 , 0.037477420750088);
		pnt[5] = Point(0.165409927389841 , 0.797112651860071);
		pnt[6] = Point(0.165409927389841 , 0.037477420750088);
		pnt[7] = Point(0.037477420750088 , 0.797112651860071);
		pnt[8] = Point(0.037477420750088 , 0.165409927389841);
		weight[0] = 0.205950504760887;
		weight[1] = 0.205950504760887;
		weight[2] = 0.205950504760887;
		weight[3] = 0.063691414286223;
		weight[4] = 0.063691414286223;
		weight[5] = 0.063691414286223;
		weight[6] = 0.063691414286223;
		weight[7] = 0.063691414286223;
		weight[8] = 0.063691414286223;
	} else if (deg==7) {
		pnt.resize(13); weight.resize(13);
		pnt[0] = Point(0.333333333333333 , 0.333333333333333);
		pnt[1] = Point(0.479308067841923 , 0.260345966079038);
		pnt[2] = Point(0.260345966079038 , 0.479308067841923);
		pnt[3] = Point(0.260345966079038 , 0.260345966079038);
		pnt[4] = Point(0.869739794195568 , 0.065130102902216);
		pnt[5] = Point(0.065130102902216 , 0.869739794195568);
		pnt[6] = Point(0.065130102902216  ,0.065130102902216);
		pnt[7] = Point(0.638444188569809 , 0.312865496004875);
		pnt[8] = Point(0.638444188569809 , 0.048690315425316);
		pnt[9] = Point(0.312865496004875 , 0.638444188569809);
		pnt[10] = Point(0.312865496004875 , 0.048690315425316);
		pnt[11] = Point(0.048690315425316 , 0.638444188569809);
		pnt[12] = Point(0.048690315425316 , 0.312865496004875);
		weight[0] = -0.149570044467670;
		weight[1] = 0.175615257433204;
		weight[2] = 0.175615257433204;
		weight[3] = 0.175615257433204;
		weight[4] = 0.053347235608839;
		weight[5] = 0.053347235608839;
		weight[6] = 0.053347235608839;
		weight[7] = 0.077113760890257;
		weight[8] = 0.077113760890257;
		weight[9] = 0.077113760890257;
		weight[10] = 0.077113760890257;
		weight[11] = 0.077113760890257;
		weight[12] = 0.077113760890257;
    } else if (deg==8) {
		pnt.resize(19); weight.resize(19);
		pnt[0] = Point(0.3333333333333333 , 0.3333333333333333);
		pnt[1] = Point(0.7974269853530872 , 0.1012865073234563);
		pnt[2] = Point(0.1012865073234563 , 0.7974269853530872);
		pnt[3] = Point(0.1012865073234563 , 0.1012865073234563);
		pnt[4] = Point(0.0597158717897698 , 0.4701420641051151);
		pnt[5] = Point(0.4701420641051151 , 0.0597158717897698);
		pnt[6] = Point(0.4701420641051151 , 0.4701420641051151);
		pnt[7] = Point(0.5357953464498992 , 0.2321023267750504);
		pnt[8] = Point(0.2321023267750504 , 0.5357953464498992);
		pnt[9] = Point(0.2321023267750504 , 0.2321023267750504);
		pnt[10] = Point(0.9410382782311209 , 0.0294808608844396);
		pnt[11] = Point(0.0294808608844396 , 0.9410382782311209);
		pnt[12] = Point(0.0294808608844396 , 0.0294808608844396);
		pnt[13] = Point(0.7384168123405100 , 0.2321023267750504);
		pnt[14] = Point(0.7384168123405100 , 0.0294808608844396);
		pnt[15] = Point(0.2321023267750504 , 0.7384168123405100);
		pnt[16] = Point(0.2321023267750504 , 0.0294808608844396);
		pnt[17] = Point(0.0294808608844396 , 0.7384168123405100);
		pnt[18] = Point(0.0294808608844396 , 0.2321023267750504);
		weight[0] = 0.0378610912003147;
		weight[1] = 0.0376204254131829;
		weight[2] = 0.0376204254131829;
		weight[3] = 0.0376204254131829;
		weight[4] = 0.0783573522441174;
		weight[5] = 0.0783573522441174;
		weight[6] = 0.0783573522441174;
		weight[7] = 0.1162714796569659;
		weight[8] = 0.1162714796569659;
		weight[9] = 0.1162714796569659;
		weight[10] = 0.0134442673751655;
		weight[11] = 0.0134442673751655;
		weight[12] = 0.0134442673751655;
		weight[13] = 0.0375097224552317;
		weight[14] = 0.0375097224552317;
		weight[15] = 0.0375097224552317;
		weight[16] = 0.0375097224552317;
		weight[17] = 0.0375097224552317;
		weight[18] = 0.0375097224552317;
    } else if (deg == 9) {
		pnt.resize(19); weight.resize(19);
		pnt[0] = Point(0.33333333333333331     ,  0.33333333333333331);
		pnt[1] = Point(2.06349616025259287E-002,  0.48968251919873701);
		pnt[2] = Point(0.48968251919873701     ,  2.06349616025259287E-002);
		pnt[3] = Point(0.48968251919873701      , 0.48968251919873701);
		pnt[4] = Point(0.12582081701412900     ,  0.43708959149293553);
		pnt[5] = Point(0.43708959149293553     ,  0.12582081701412900);
		pnt[6] = Point(0.43708959149293553     ,  0.43708959149293553);
		pnt[7] = Point(0.62359292876193562     ,  0.18820353561903219);
		pnt[8] = Point(0.18820353561903219     ,  0.62359292876193562);
		pnt[9] = Point(0.18820353561903219     ,  0.18820353561903219);
		pnt[10] = Point(0.91054097321109406     ,  4.47295133944529688E-002);
		pnt[11] = Point(4.47295133944529688E-002,  0.91054097321109406);
		pnt[12] = Point(4.47295133944529688E-002,  4.47295133944529688E-002);
		pnt[13] = Point(0.74119859878449801     ,  3.68384120547362581E-002);
		pnt[14] = Point(0.74119859878449801     ,  0.22196298916076573);
		pnt[15] = Point(3.68384120547362581E-002,  0.74119859878449801);
		pnt[16] = Point(3.68384120547362581E-002,  0.22196298916076573);
		pnt[17] = Point(0.22196298916076573     ,  0.74119859878449801);
		pnt[18] = Point(0.22196298916076573     ,  3.68384120547362581E-002 );
		weight[0] = 9.71357962827961025E-002;
		weight[1] = 3.13347002271398278E-002;
		weight[2] = 3.13347002271398278E-002;
		weight[3] = 3.13347002271398278E-002;
		weight[4] = 7.78275410047754301E-002;
		weight[5] = 7.78275410047754301E-002;
		weight[6] = 7.78275410047754301E-002;
		weight[7] = 7.96477389272090969E-002;
		weight[8] = 7.96477389272090969E-002;
		weight[9] = 7.96477389272090969E-002;
		weight[10] = 2.55776756586981006E-002;
		weight[11] = 2.55776756586981006E-002;
		weight[12] = 2.55776756586981006E-002;
		weight[13] = 4.32835393772893970E-002;
		weight[14] = 4.32835393772893970E-002;
		weight[15] = 4.32835393772893970E-002;
		weight[16] = 4.32835393772893970E-002;
		weight[17] = 4.32835393772893970E-002;
		weight[18] = 4.32835393772893970E-002;
    } else if (deg<=11) {
		pnt.resize(28); weight.resize(28);
		pnt[0] = Point(0.33333333333333333 , 0.333333333333333333);
		pnt[1] = Point(0.9480217181434233  , 0.02598914092828833);
		pnt[2] = Point(0.02598914092828833 , 0.9480217181434233);
		pnt[3] = Point(0.02598914092828833 , 0.02598914092828833);
		pnt[4] = Point(0.8114249947041546  , 0.09428750264792270);
		pnt[5] = Point(0.09428750264792270 , 0.8114249947041546);
		pnt[6] = Point(0.09428750264792270 , 0.09428750264792270);
		pnt[7] = Point(0.01072644996557060 , 0.4946367750172147);
		pnt[8] = Point(0.4946367750172147  , 0.01072644996557060);
		pnt[9] = Point(0.4946367750172147  , 0.4946367750172147);
		pnt[10] = Point(0.5853132347709715  , 0.2073433826145142);
		pnt[11] = Point(0.2073433826145142  , 0.5853132347709715);
		pnt[12] = Point(0.2073433826145142  , 0.2073433826145142);
		pnt[13] = Point(0.1221843885990187  , 0.4389078057004907);
		pnt[14] = Point(0.4389078057004907  , 0.1221843885990187);
		pnt[15] = Point(0.4389078057004907  , 0.4389078057004907);
		pnt[16] = Point(0.6779376548825902  , 0.04484167758913055);
		pnt[17] = Point(0.6779376548825902  , 0.27722066752827925);
		pnt[18] = Point(0.04484167758913055 , 0.6779376548825902);
		pnt[19] = Point(0.04484167758913055 , 0.27722066752827925);
		pnt[20] = Point(0.27722066752827925 , 0.6779376548825902);
		pnt[21] = Point(0.27722066752827925 , 0.04484167758913055);
		pnt[22] = Point(0.8588702812826364  , 0.00000000000000000);
		pnt[23] = Point(0.8588702812826364  , 0.1411297187173636);
		pnt[24] = Point(0.0000000000000000  , 0.8588702812826364);
		pnt[25] = Point(0.0000000000000000  , 0.1411297187173636);
		pnt[26] = Point(0.1411297187173636  , 0.8588702812826364);
		pnt[27] = Point(0.1411297187173636  , 0.0000000000000000);
		weight[0] = 0.08797730116222190;
		weight[1] = 0.008744311553736190;
		weight[2] = 0.008744311553736190;
		weight[3] = 0.008744311553736190;
		weight[4] = 0.03808157199393533;
		weight[5] = 0.03808157199393533;
		weight[6] = 0.03808157199393533;
		weight[7] = 0.01885544805613125;
		weight[8] = 0.01885544805613125;
		weight[9] = 0.01885544805613125;
		weight[10] = 0.07215969754474100;
		weight[11] = 0.07215969754474100;
		weight[12] = 0.07215969754474100;
		weight[13] = 0.06932913870553720;
		weight[14] = 0.06932913870553720;
		weight[15] = 0.06932913870553720;
		weight[16] = 0.04105631542928860;
		weight[17] = 0.04105631542928860;
		weight[18] = 0.04105631542928860;
		weight[19] = 0.04105631542928860;
		weight[20] = 0.04105631542928860;
		weight[21] = 0.04105631542928860;
		weight[22] = 0.007362383783300573;
		weight[23] = 0.007362383783300573;
		weight[24] = 0.007362383783300573;
		weight[25] = 0.007362383783300573;
		weight[26] = 0.007362383783300573;
		weight[27] = 0.007362383783300573;
    } else if (deg<=13) {
    	pnt.resize(37); weight.resize(37);
    	pnt[0] = Point(0.333333333333333333333333333333,  0.333333333333333333333333333333);
    	pnt[1] = Point(0.950275662924105565450352089520,  0.024862168537947217274823955239);
    	pnt[2] = Point(0.024862168537947217274823955239,  0.950275662924105565450352089520);
    	pnt[3] = Point(0.024862168537947217274823955239,  0.024862168537947217274823955239);
    	pnt[4] = Point(0.171614914923835347556304795551,  0.414192542538082326221847602214);
    	pnt[5] = Point(0.414192542538082326221847602214,  0.171614914923835347556304795551);
    	pnt[6] = Point(0.414192542538082326221847602214,  0.414192542538082326221847602214);
    	pnt[7] = Point(0.539412243677190440263092985511,  0.230293878161404779868453507244);
    	pnt[8] = Point(0.230293878161404779868453507244,  0.539412243677190440263092985511);
    	pnt[9] = Point(0.230293878161404779868453507244,  0.230293878161404779868453507244);
    	pnt[10] = Point(0.772160036676532561750285570113,  0.113919981661733719124857214943);
    	pnt[11] = Point(0.113919981661733719124857214943,  0.772160036676532561750285570113);
    	pnt[12] = Point(0.113919981661733719124857214943,  0.113919981661733719124857214943);
    	pnt[13] = Point(0.009085399949835353883572964740,  0.495457300025082323058213517632);
    	pnt[14] = Point(0.495457300025082323058213517632,  0.009085399949835353883572964740);
    	pnt[15] = Point(0.495457300025082323058213517632,  0.495457300025082323058213517632);
    	pnt[16] = Point(0.062277290305886993497083640527,  0.468861354847056503251458179727);
    	pnt[17] = Point(0.468861354847056503251458179727,  0.062277290305886993497083640527);
    	pnt[18] = Point(0.468861354847056503251458179727,  0.468861354847056503251458179727);
    	pnt[19] = Point(0.022076289653624405142446876931,  0.851306504174348550389457672223);
    	pnt[20] = Point(0.022076289653624405142446876931,  0.126617206172027096933163647918);
    	pnt[21] = Point(0.851306504174348550389457672223,  0.022076289653624405142446876931);
    	pnt[22] = Point(0.851306504174348550389457672223,  0.126617206172027096933163647918);
    	pnt[23] = Point(0.126617206172027096933163647918,  0.022076289653624405142446876931);
    	pnt[24] = Point(0.126617206172027096933163647918,  0.851306504174348550389457672223);
    	pnt[25] = Point(0.018620522802520968955913511549,  0.689441970728591295496647976487);
    	pnt[26] = Point(0.018620522802520968955913511549,  0.291937506468887771754472382212);
    	pnt[27] = Point(0.689441970728591295496647976487,  0.018620522802520968955913511549);
    	pnt[28] = Point(0.689441970728591295496647976487,  0.291937506468887771754472382212);
    	pnt[29] = Point(0.291937506468887771754472382212,  0.018620522802520968955913511549);
    	pnt[30] = Point(0.291937506468887771754472382212,  0.689441970728591295496647976487);
    	pnt[31] = Point(0.096506481292159228736516560903,  0.635867859433872768286976979827);
    	pnt[32] = Point(0.096506481292159228736516560903,  0.267625659273967961282458816185);
    	pnt[33] = Point(0.635867859433872768286976979827,  0.096506481292159228736516560903);
    	pnt[34] = Point(0.635867859433872768286976979827,  0.267625659273967961282458816185);
    	pnt[35] = Point(0.267625659273967961282458816185,  0.096506481292159228736516560903);
    	pnt[36] = Point(0.267625659273967961282458816185,  0.635867859433872768286976979827);

        weight[0] = 0.051739766065744133555179145422;
        weight[1] = 0.008007799555564801597804123460;
        weight[2] = 0.008007799555564801597804123460;
        weight[3] = 0.008007799555564801597804123460;
        weight[4] = 0.046868898981821644823226732071;
        weight[5] = 0.046868898981821644823226732071;
        weight[6] = 0.046868898981821644823226732071;
        weight[7] = 0.046590940183976487960361770070;
        weight[8] = 0.046590940183976487960361770070;
        weight[9] = 0.046590940183976487960361770070;
        weight[10] = 0.031016943313796381407646220131;
        weight[11] = 0.031016943313796381407646220131;
        weight[12] = 0.031016943313796381407646220131;
        weight[13] = 0.010791612736631273623178240136;
        weight[14] = 0.010791612736631273623178240136;
        weight[15] = 0.010791612736631273623178240136;
        weight[16] = 0.032195534242431618819414482205;
        weight[17] = 0.032195534242431618819414482205;
        weight[18] = 0.032195534242431618819414482205;
        weight[19] = 0.015445834210701583817692900053;
        weight[20] = 0.015445834210701583817692900053;
        weight[21] = 0.015445834210701583817692900053;
        weight[22] = 0.015445834210701583817692900053;
        weight[23] = 0.015445834210701583817692900053;
        weight[24] = 0.015445834210701583817692900053;
        weight[25] = 0.017822989923178661888748319485;
        weight[26] = 0.017822989923178661888748319485;
        weight[27] = 0.017822989923178661888748319485;
        weight[28] = 0.017822989923178661888748319485;
        weight[29] = 0.017822989923178661888748319485;
        weight[30] = 0.017822989923178661888748319485;
        weight[31] = 0.037038683681384627918546472190;
        weight[32] = 0.037038683681384627918546472190;
        weight[33] = 0.037038683681384627918546472190;
        weight[34] = 0.037038683681384627918546472190;
        weight[35] = 0.037038683681384627918546472190;
        weight[36] = 0.037038683681384627918546472190;
    } else if (deg<=14) {
    	pnt.resize(46); weight.resize(46);
        pnt[0] = Point(0.33333333333333331000 , 0.33333333333333331000);
        pnt[1] = Point(0.00997976080645843200 , 0.00997976080645843200);
        pnt[2] = Point(0.00997976080645843200 , 0.98004047838708308000);
        pnt[3] = Point(0.98004047838708308000 , 0.00997976080645843200);
        pnt[4] = Point(0.47997789352118841000 , 0.47997789352118841000);
        pnt[5] = Point(0.47997789352118841000 , 0.04004421295762317100);
        pnt[6] = Point(0.04004421295762317100 , 0.47997789352118841000);
        pnt[7] = Point(0.15381195917696691000 , 0.15381195917696691000);
        pnt[8] = Point(0.15381195917696691000 , 0.69237608164606623000);
        pnt[9] = Point(0.69237608164606623000 , 0.15381195917696691000);
        pnt[10] = Point(0.07402347711698781300 , 0.07402347711698781300);
        pnt[11] = Point(0.07402347711698781300 , 0.85195304576602437000);
        pnt[12] = Point(0.85195304576602437000 , 0.07402347711698781300);
        pnt[13] = Point(0.13035468250332999000 , 0.13035468250332999000);
        pnt[14] = Point(0.13035468250332999000 , 0.73929063499334002000);
        pnt[15] = Point(0.73929063499334002000 , 0.13035468250332999000);
        pnt[16] = Point(0.23061722602665313000 , 0.23061722602665313000);
        pnt[17] = Point(0.23061722602665313000 , 0.53876554794669373000);
        pnt[18] = Point(0.53876554794669373000 , 0.23061722602665313000);
        pnt[19] = Point(0.42233208341914780000 , 0.42233208341914780000);
        pnt[20] = Point(0.42233208341914780000 , 0.15533583316170441000);
        pnt[21] = Point(0.15533583316170441000 , 0.42233208341914780000);
        pnt[22] = Point(0.78623738593466097000 , 0.19061636003190091000);
        pnt[23] = Point(0.78623738593466097000 , 0.02314625403343817400);
        pnt[24] = Point(0.19061636003190091000 , 0.78623738593466097000);
        pnt[25] = Point(0.19061636003190091000 , 0.02314625403343817400);
        pnt[26] = Point(0.02314625403343817400 , 0.78623738593466097000);
        pnt[27] = Point(0.02314625403343817400 , 0.19061636003190091000);
        pnt[28] = Point(0.63055214366060741000 , 0.36232313774354713000);
        pnt[29] = Point(0.63055214366060741000 , 0.00712471859584540290);
        pnt[30] = Point(0.36232313774354713000 , 0.63055214366060741000);
        pnt[31] = Point(0.36232313774354713000 , 0.00712471859584540290);
        pnt[32] = Point(0.00712471859584540290 , 0.63055214366060741000);
        pnt[33] = Point(0.00712471859584540290 , 0.36232313774354713000);
        pnt[34] = Point(0.62657732985630632000 , 0.29077120588366739000);
        pnt[35] = Point(0.62657732985630632000 , 0.08265146426002623100);
        pnt[36] = Point(0.29077120588366739000 , 0.62657732985630632000);
        pnt[37] = Point(0.29077120588366739000 , 0.08265146426002623100);
        pnt[38] = Point(0.08265146426002623100 , 0.62657732985630632000);
        pnt[39] = Point(0.08265146426002623100 , 0.29077120588366739000);
        pnt[40] = Point(0.91420998492962546000 , 0.07116571087775076800);
        pnt[41] = Point(0.91420998492962546000 , 0.01462430419262372700);
        pnt[42] = Point(0.07116571087775076800 , 0.91420998492962546000);
        pnt[43] = Point(0.07116571087775076800 , 0.01462430419262372700);
        pnt[44] = Point(0.01462430419262372700 , 0.91420998492962546000);
        pnt[45] = Point(0.01462430419262372700 , 0.07116571087775076800);

        weight[0] = 0.05859628522602859700;
        weight[1] = 0.00173515122972526760;
        weight[2] = 0.00173515122972526760;
        weight[3] = 0.00173515122972526760;
        weight[4] = 0.02616378255861452300;
        weight[5] = 0.02616378255861452300;
        weight[6] = 0.02616378255861452300;
        weight[7] = 0.00391972924240182890;
        weight[8] = 0.00391972924240182890;
        weight[9] = 0.00391972924240182890;
        weight[10] = 0.01224735975694086700;
        weight[11] = 0.01224735975694086700;
        weight[12] = 0.01224735975694086700;
        weight[13] = 0.02819962850325796000;
        weight[14] = 0.02819962850325796000;
        weight[15] = 0.02819962850325796000;
        weight[16] = 0.05088708718595948800;
        weight[17] = 0.05088708718595948800;
        weight[18] = 0.05088708718595948800;
        weight[19] = 0.05045343990160360000;
        weight[20] = 0.05045343990160360000;
        weight[21] = 0.05045343990160360000;
        weight[22] = 0.01706364421223345200;
        weight[23] = 0.01706364421223345200;
        weight[24] = 0.01706364421223345200;
        weight[25] = 0.01706364421223345200;
        weight[26] = 0.01706364421223345200;
        weight[27] = 0.01706364421223345200;
        weight[28] = 0.00968346642550660040;
        weight[29] = 0.00968346642550660040;
        weight[30] = 0.00968346642550660040;
        weight[31] = 0.00968346642550660040;
        weight[32] = 0.00968346642550660040;
        weight[33] = 0.00968346642550660040;
        weight[34] = 0.03638575592848500300;
        weight[35] = 0.03638575592848500300;
        weight[36] = 0.03638575592848500300;
        weight[37] = 0.03638575592848500300;
        weight[38] = 0.03638575592848500300;
        weight[39] = 0.03638575592848500300;
        weight[40] = 0.00696466337351841270;
        weight[41] = 0.00696466337351841270;
        weight[42] = 0.00696466337351841270;
        weight[43] = 0.00696466337351841270;
        weight[44] = 0.00696466337351841270;
        weight[45] = 0.00696466337351841270;
    } else if (deg<=20) {
    	pnt.resize(88); weight.resize(88);
        pnt[0] = Point(0.33333333333333331000 , 0.33333333333333331000);
        pnt[1] = Point(0.21587430593299198000 , 0.21587430593299198000);
        pnt[2] = Point(0.21587430593299198000 , 0.56825138813401610000);
        pnt[3] = Point(0.56825138813401610000 , 0.21587430593299198000);
        pnt[4] = Point(0.07537676652974727200 , 0.07537676652974727200);
        pnt[5] = Point(0.07537676652974727200 , 0.84924646694050543000);
        pnt[6] = Point(0.84924646694050543000 , 0.07537676652974727200);
        pnt[7] = Point(0.01030082813722179300 , 0.01030082813722179300);
        pnt[8] = Point(0.01030082813722179300 , 0.97939834372555645000);
        pnt[9] = Point(0.97939834372555645000 , 0.01030082813722179300);
        pnt[10] = Point(0.49360221129870019000 , 0.49360221129870019000);
        pnt[11] = Point(0.49360221129870019000 , 0.01279557740259962300);
        pnt[12] = Point(0.01279557740259962300 , 0.49360221129870019000);
        pnt[13] = Point(0.46155093810692532000 , 0.46155093810692532000);
        pnt[14] = Point(0.46155093810692532000 , 0.07689812378614935300);
        pnt[15] = Point(0.07689812378614935300 , 0.46155093810692532000);
        pnt[16] = Point(0.32862140642423698000 , 0.42934057025821037000);
        pnt[17] = Point(0.32862140642423698000 , 0.24203802331755264000);
        pnt[18] = Point(0.42934057025821037000 , 0.32862140642423698000);
        pnt[19] = Point(0.42934057025821037000 , 0.24203802331755264000);
        pnt[20] = Point(0.24203802331755264000 , 0.32862140642423698000);
        pnt[21] = Point(0.24203802331755264000 , 0.42934057025821037000);
        pnt[22] = Point(0.26048036178656875000 , 0.10157753428096944000);
        pnt[23] = Point(0.26048036178656875000 , 0.63794210393246176000);
        pnt[24] = Point(0.10157753428096944000 , 0.26048036178656875000);
        pnt[25] = Point(0.10157753428096944000 , 0.63794210393246176000);
        pnt[26] = Point(0.63794210393246176000 , 0.26048036178656875000);
        pnt[27] = Point(0.63794210393246176000 , 0.10157753428096944000);
        pnt[28] = Point(0.13707423584645531000 , 0.71006597300113017000);
        pnt[29] = Point(0.13707423584645531000 , 0.15285979115241455000);
        pnt[30] = Point(0.71006597300113017000 , 0.13707423584645531000);
        pnt[31] = Point(0.71006597300113017000 , 0.15285979115241455000);
        pnt[32] = Point(0.15285979115241455000 , 0.13707423584645531000);
        pnt[33] = Point(0.15285979115241455000 , 0.71006597300113017000);
        pnt[34] = Point(0.14672694587229979000 , 0.49854547767841484000);
        pnt[35] = Point(0.14672694587229979000 , 0.35472757644928543000);
        pnt[36] = Point(0.49854547767841484000 , 0.14672694587229979000);
        pnt[37] = Point(0.49854547767841484000 , 0.35472757644928543000);
        pnt[38] = Point(0.35472757644928543000 , 0.14672694587229979000);
        pnt[39] = Point(0.35472757644928543000 , 0.49854547767841484000);
        pnt[40] = Point(0.02699897774255329000 , 0.04918672267258199900);
        pnt[41] = Point(0.02699897774255329000 , 0.92381429958486472000);
        pnt[42] = Point(0.04918672267258199900 , 0.02699897774255329000);
        pnt[43] = Point(0.04918672267258199900 , 0.92381429958486472000);
        pnt[44] = Point(0.92381429958486472000 , 0.02699897774255329000);
        pnt[45] = Point(0.92381429958486472000 , 0.04918672267258199900);
        pnt[46] = Point(0.06187178593361702900 , 0.77966014654056937000);
        pnt[47] = Point(0.06187178593361702900 , 0.15846806752581366000);
        pnt[48] = Point(0.77966014654056937000 , 0.06187178593361702900);
        pnt[49] = Point(0.77966014654056937000 , 0.15846806752581366000);
        pnt[50] = Point(0.15846806752581366000 , 0.06187178593361702900);
        pnt[51] = Point(0.15846806752581366000 , 0.77966014654056937000);
        pnt[52] = Point(0.04772436742762199700 , 0.37049153914954763000);
        pnt[53] = Point(0.04772436742762199700 , 0.58178409342283044000);
        pnt[54] = Point(0.37049153914954763000 , 0.04772436742762199700);
        pnt[55] = Point(0.37049153914954763000 , 0.58178409342283044000);
        pnt[56] = Point(0.58178409342283044000 , 0.04772436742762199700);
        pnt[57] = Point(0.58178409342283044000 , 0.37049153914954763000);
        pnt[58] = Point(0.12060051518636437000 , 0.86334694875475260000);
        pnt[59] = Point(0.12060051518636437000 , 0.01605253605888301600);
        pnt[60] = Point(0.86334694875475260000 , 0.12060051518636437000);
        pnt[61] = Point(0.86334694875475260000 , 0.01605253605888301600);
        pnt[62] = Point(0.01605253605888301600 , 0.12060051518636437000);
        pnt[63] = Point(0.01605253605888301600 , 0.86334694875475260000);
        pnt[64] = Point(0.00269714779670978760 , 0.05619493818774550000);
        pnt[65] = Point(0.00269714779670978760 , 0.94110791401554472000);
        pnt[66] = Point(0.05619493818774550000 , 0.00269714779670978760);
        pnt[67] = Point(0.05619493818774550000 , 0.94110791401554472000);
        pnt[68] = Point(0.94110791401554472000 , 0.00269714779670978760);
        pnt[69] = Point(0.94110791401554472000 , 0.05619493818774550000);
        pnt[70] = Point(0.00301563327794236250 , 0.20867500674842135000);
        pnt[71] = Point(0.00301563327794236250 , 0.78830935997363627000);
        pnt[72] = Point(0.20867500674842135000 , 0.00301563327794236250);
        pnt[73] = Point(0.20867500674842135000 , 0.78830935997363627000);
        pnt[74] = Point(0.78830935997363627000 , 0.00301563327794236250);
        pnt[75] = Point(0.78830935997363627000 , 0.20867500674842135000);
        pnt[76] = Point(0.02990537578845702000 , 0.72115124091203409000);
        pnt[77] = Point(0.02990537578845702000 , 0.24894338329950894000);
        pnt[78] = Point(0.72115124091203409000 , 0.02990537578845702000);
        pnt[79] = Point(0.72115124091203409000 , 0.24894338329950894000);
        pnt[80] = Point(0.24894338329950894000 , 0.02990537578845702000);
        pnt[81] = Point(0.24894338329950894000 , 0.72115124091203409000);
        pnt[82] = Point(0.00675665422246098880 , 0.64005544194054187000);
        pnt[83] = Point(0.00675665422246098880 , 0.35318790383699716000);
        pnt[84] = Point(0.64005544194054187000 , 0.00675665422246098880);
        pnt[85] = Point(0.64005544194054187000 , 0.35318790383699716000);
        pnt[86] = Point(0.35318790383699716000 , 0.00675665422246098880);
        pnt[87] = Point(0.35318790383699716000 , 0.64005544194054187000);

        weight[0] = 0.01253760799449665600;
        weight[1] = 0.02747186987642421400;
        weight[2] = 0.02747186987642421400;
        weight[3] = 0.02747186987642421400;
        weight[4] = 0.00976527227705142370;
        weight[5] = 0.00976527227705142370;
        weight[6] = 0.00976527227705142370;
        weight[7] = 0.00139841953539182350;
        weight[8] = 0.00139841953539182350;
        weight[9] = 0.00139841953539182350;
        weight[10] = 0.00929210262518518310;
        weight[11] = 0.00929210262518518310;
        weight[12] = 0.00929210262518518310;
        weight[13] = 0.01657787603236692700;
        weight[14] = 0.01657787603236692700;
        weight[15] = 0.01657787603236692700;
        weight[16] = 0.02066776234866507900;
        weight[17] = 0.02066776234866507900;
        weight[18] = 0.02066776234866507900;
        weight[19] = 0.02066776234866507900;
        weight[20] = 0.02066776234866507900;
        weight[21] = 0.02066776234866507900;
        weight[22] = 0.02082223552115450600;
        weight[23] = 0.02082223552115450600;
        weight[24] = 0.02082223552115450600;
        weight[25] = 0.02082223552115450600;
        weight[26] = 0.02082223552115450600;
        weight[27] = 0.02082223552115450600;
        weight[28] = 0.00956863841984906090;
        weight[29] = 0.00956863841984906090;
        weight[30] = 0.00956863841984906090;
        weight[31] = 0.00956863841984906090;
        weight[32] = 0.00956863841984906090;
        weight[33] = 0.00956863841984906090;
        weight[34] = 0.02445277096897246300;
        weight[35] = 0.02445277096897246300;
        weight[36] = 0.02445277096897246300;
        weight[37] = 0.02445277096897246300;
        weight[38] = 0.02445277096897246300;
        weight[39] = 0.02445277096897246300;
        weight[40] = 0.00315573063063053420;
        weight[41] = 0.00315573063063053420;
        weight[42] = 0.00315573063063053420;
        weight[43] = 0.00315573063063053420;
        weight[44] = 0.00315573063063053420;
        weight[45] = 0.00315573063063053420;
        weight[46] = 0.01213679636532129800;
        weight[47] = 0.01213679636532129800;
        weight[48] = 0.01213679636532129800;
        weight[49] = 0.01213679636532129800;
        weight[50] = 0.01213679636532129800;
        weight[51] = 0.01213679636532129800;
        weight[52] = 0.01496648014388644900;
        weight[53] = 0.01496648014388644900;
        weight[54] = 0.01496648014388644900;
        weight[55] = 0.01496648014388644900;
        weight[56] = 0.01496648014388644900;
        weight[57] = 0.01496648014388644900;
        weight[58] = 0.00632759332177773930;
        weight[59] = 0.00632759332177773930;
        weight[60] = 0.00632759332177773930;
        weight[61] = 0.00632759332177773930;
        weight[62] = 0.00632759332177773930;
        weight[63] = 0.00632759332177773930;
        weight[64] = 0.00134256031206369590;
        weight[65] = 0.00134256031206369590;
        weight[66] = 0.00134256031206369590;
        weight[67] = 0.00134256031206369590;
        weight[68] = 0.00134256031206369590;
        weight[69] = 0.00134256031206369590;
        weight[70] = 0.00277607691634755400;
        weight[71] = 0.00277607691634755400;
        weight[72] = 0.00277607691634755400;
        weight[73] = 0.00277607691634755400;
        weight[74] = 0.00277607691634755400;
        weight[75] = 0.00277607691634755400;
        weight[76] = 0.01073984447418494100;
        weight[77] = 0.01073984447418494100;
        weight[78] = 0.01073984447418494100;
        weight[79] = 0.01073984447418494100;
        weight[80] = 0.01073984447418494100;
        weight[81] = 0.01073984447418494100;
        weight[82] = 0.00536780573818745280;
        weight[83] = 0.00536780573818745280;
        weight[84] = 0.00536780573818745280;
        weight[85] = 0.00536780573818745280;
        weight[86] = 0.00536780573818745280;
        weight[87] = 0.00536780573818745280;
    } else if (deg<=25) {
    	pnt.resize(126); weight.resize(126);
        pnt[0] = Point(0.02794648307316999900 , 0.48602675846340998000);
        pnt[1] = Point(0.48602675846340998000 , 0.48602675846340998000);
        pnt[2] = Point(0.48602675846340998000 , 0.02794648307316999900);
        pnt[3] = Point(0.13117860132765000000 , 0.43441069933616999000);
        pnt[4] = Point(0.43441069933616999000 , 0.43441069933616999000);
        pnt[5] = Point(0.43441069933616999000 , 0.13117860132765000000);
        pnt[6] = Point(0.22022172951207000000 , 0.38988913524396002000);
        pnt[7] = Point(0.38988913524396002000 , 0.38988913524396002000);
        pnt[8] = Point(0.38988913524396002000 , 0.22022172951207000000);
        pnt[9] = Point(0.40311353196039001000 , 0.29844323401980000000);
        pnt[10] = Point(0.29844323401980000000 , 0.29844323401980000000);
        pnt[11] = Point(0.29844323401980000000 , 0.40311353196039001000);
        pnt[12] = Point(0.53191165532525997000 , 0.23404417233736999000);
        pnt[13] = Point(0.23404417233736999000 , 0.23404417233736999000);
        pnt[14] = Point(0.23404417233736999000 , 0.53191165532525997000);
        pnt[15] = Point(0.69706333078196003000 , 0.15146833460902001000);
        pnt[16] = Point(0.15146833460902001000 , 0.15146833460902001000);
        pnt[17] = Point(0.15146833460902001000 , 0.69706333078196003000);
        pnt[18] = Point(0.77453221290801000000 , 0.11273389354599000000);
        pnt[19] = Point(0.11273389354599000000 , 0.11273389354599000000);
        pnt[20] = Point(0.11273389354599000000 , 0.77453221290801000000);
        pnt[21] = Point(0.84456861581694997000 , 0.07771569209152999500);
        pnt[22] = Point(0.07771569209152999500 , 0.07771569209152999500);
        pnt[23] = Point(0.07771569209152999500 , 0.84456861581694997000);
        pnt[24] = Point(0.93021381277141002000 , 0.03489309361430000000);
        pnt[25] = Point(0.03489309361430000000 , 0.03489309361430000000);
        pnt[26] = Point(0.03489309361430000000 , 0.93021381277141002000);
        pnt[27] = Point(0.98548363075812995000 , 0.00725818462093000010);
        pnt[28] = Point(0.00725818462093000010 , 0.00725818462093000010);
        pnt[29] = Point(0.00725818462093000010 , 0.98548363075812995000);
        pnt[30] = Point(0.00129235270443999990 , 0.22721445215336000000);
        pnt[31] = Point(0.22721445215336000000 , 0.77149319514218995000);
        pnt[32] = Point(0.77149319514218995000 , 0.00129235270443999990);
        pnt[33] = Point(0.22721445215336000000 , 0.00129235270443999990);
        pnt[34] = Point(0.77149319514218995000 , 0.22721445215336000000);
        pnt[35] = Point(0.00129235270443999990 , 0.77149319514218995000);
        pnt[36] = Point(0.00539970127212000000 , 0.43501055485356999000);
        pnt[37] = Point(0.43501055485356999000 , 0.55958974387431004000);
        pnt[38] = Point(0.55958974387431004000 , 0.00539970127212000000);
        pnt[39] = Point(0.43501055485356999000 , 0.00539970127212000000);
        pnt[40] = Point(0.55958974387431004000 , 0.43501055485356999000);
        pnt[41] = Point(0.00539970127212000000 , 0.55958974387431004000);
        pnt[42] = Point(0.00638400303398000030 , 0.32030959927219999000);
        pnt[43] = Point(0.32030959927219999000 , 0.67330639769381995000);
        pnt[44] = Point(0.67330639769381995000 , 0.00638400303398000030);
        pnt[45] = Point(0.32030959927219999000 , 0.00638400303398000030);
        pnt[46] = Point(0.67330639769381995000 , 0.32030959927219999000);
        pnt[47] = Point(0.00638400303398000030 , 0.67330639769381995000);
        pnt[48] = Point(0.00502821150199000020 , 0.09175032228000999700);
        pnt[49] = Point(0.09175032228000999700 , 0.90322146621800004000);
        pnt[50] = Point(0.90322146621800004000 , 0.00502821150199000020);
        pnt[51] = Point(0.09175032228000999700 , 0.00502821150199000020);
        pnt[52] = Point(0.90322146621800004000 , 0.09175032228000999700);
        pnt[53] = Point(0.00502821150199000020 , 0.90322146621800004000);
        pnt[54] = Point(0.00682675862178000040 , 0.03801083585871999800);
        pnt[55] = Point(0.03801083585871999800 , 0.95516240551949005000);
        pnt[56] = Point(0.95516240551949005000 , 0.00682675862178000040);
        pnt[57] = Point(0.03801083585871999800 , 0.00682675862178000040);
        pnt[58] = Point(0.95516240551949005000 , 0.03801083585871999800);
        pnt[59] = Point(0.00682675862178000040 , 0.95516240551949005000);
        pnt[60] = Point(0.01001619963993000000 , 0.15742521848530999000);
        pnt[61] = Point(0.15742521848530999000 , 0.83255858187475995000);
        pnt[62] = Point(0.83255858187475995000 , 0.01001619963993000000);
        pnt[63] = Point(0.15742521848530999000 , 0.01001619963993000000);
        pnt[64] = Point(0.83255858187475995000 , 0.15742521848530999000);
        pnt[65] = Point(0.01001619963993000000 , 0.83255858187475995000);
        pnt[66] = Point(0.02575781317339000100 , 0.23988965977853000000);
        pnt[67] = Point(0.23988965977853000000 , 0.73435252704807996000);
        pnt[68] = Point(0.73435252704807996000 , 0.02575781317339000100);
        pnt[69] = Point(0.23988965977853000000 , 0.02575781317339000100);
        pnt[70] = Point(0.73435252704807996000 , 0.23988965977853000000);
        pnt[71] = Point(0.02575781317339000100 , 0.73435252704807996000);
        pnt[72] = Point(0.03022789811992000100 , 0.36194311812606000000);
        pnt[73] = Point(0.36194311812606000000 , 0.60782898375401995000);
        pnt[74] = Point(0.60782898375401995000 , 0.03022789811992000100);
        pnt[75] = Point(0.36194311812606000000 , 0.03022789811992000100);
        pnt[76] = Point(0.60782898375401995000 , 0.36194311812606000000);
        pnt[77] = Point(0.03022789811992000100 , 0.60782898375401995000);
        pnt[78] = Point(0.03050499010715999900 , 0.08355196095483000100);
        pnt[79] = Point(0.08355196095483000100 , 0.88594304893801001000);
        pnt[80] = Point(0.88594304893801001000 , 0.03050499010715999900);
        pnt[81] = Point(0.08355196095483000100 , 0.03050499010715999900);
        pnt[82] = Point(0.88594304893801001000 , 0.08355196095483000100);
        pnt[83] = Point(0.03050499010715999900 , 0.88594304893801001000);
        pnt[84] = Point(0.04595654736257000200 , 0.14844322073242000000);
        pnt[85] = Point(0.14844322073242000000 , 0.80560023190500996000);
        pnt[86] = Point(0.80560023190500996000 , 0.04595654736257000200);
        pnt[87] = Point(0.14844322073242000000 , 0.04595654736257000200);
        pnt[88] = Point(0.80560023190500996000 , 0.14844322073242000000);
        pnt[89] = Point(0.04595654736257000200 , 0.80560023190500996000);
        pnt[90] = Point(0.06744280054027999800 , 0.28373970872753002000);
        pnt[91] = Point(0.28373970872753002000 , 0.64881749073218997000);
        pnt[92] = Point(0.64881749073218997000 , 0.06744280054027999800);
        pnt[93] = Point(0.28373970872753002000 , 0.06744280054027999800);
        pnt[94] = Point(0.64881749073218997000 , 0.28373970872753002000);
        pnt[95] = Point(0.06744280054027999800 , 0.64881749073218997000);
        pnt[96] = Point(0.07004509141591000500 , 0.40689937511878999000);
        pnt[97] = Point(0.40689937511878999000 , 0.52305553346529998000);
        pnt[98] = Point(0.52305553346529998000 , 0.07004509141591000500);
        pnt[99] = Point(0.40689937511878999000 , 0.07004509141591000500);
        pnt[100] = Point(0.52305553346529998000 , 0.40689937511878999000);
        pnt[101] = Point(0.07004509141591000500 , 0.52305553346529998000);
        pnt[102] = Point(0.08391152464012000000 , 0.19411398702488999000);
        pnt[103] = Point(0.19411398702488999000 , 0.72197448833499001000);
        pnt[104] = Point(0.72197448833499001000 , 0.08391152464012000000);
        pnt[105] = Point(0.19411398702488999000 , 0.08391152464012000000);
        pnt[106] = Point(0.72197448833499001000 , 0.19411398702488999000);
        pnt[107] = Point(0.08391152464012000000 , 0.72197448833499001000);
        pnt[108] = Point(0.12037553567714999000 , 0.32413434700069998000);
        pnt[109] = Point(0.32413434700069998000 , 0.55549011732214004000);
        pnt[110] = Point(0.55549011732214004000 , 0.12037553567714999000);
        pnt[111] = Point(0.32413434700069998000 , 0.12037553567714999000);
        pnt[112] = Point(0.55549011732214004000 , 0.32413434700069998000);
        pnt[113] = Point(0.12037553567714999000 , 0.55549011732214004000);
        pnt[114] = Point(0.14806689915737001000 , 0.22927748355597999000);
        pnt[115] = Point(0.22927748355597999000 , 0.62265561728664998000);
        pnt[116] = Point(0.62265561728664998000 , 0.14806689915737001000);
        pnt[117] = Point(0.22927748355597999000 , 0.14806689915737001000);
        pnt[118] = Point(0.62265561728664998000 , 0.22927748355597999000);
        pnt[119] = Point(0.14806689915737001000 , 0.62265561728664998000);
        pnt[120] = Point(0.19177186586733000000 , 0.32561812259598000000);
        pnt[121] = Point(0.32561812259598000000 , 0.48261001153668998000);
        pnt[122] = Point(0.48261001153668998000 , 0.19177186586733000000);
        pnt[123] = Point(0.32561812259598000000 , 0.19177186586733000000);
        pnt[124] = Point(0.48261001153668998000 , 0.32561812259598000000);
        pnt[125] = Point(0.19177186586733000000 , 0.48261001153668998000);
        weight[0] = 0.00800558188002041710;
        weight[1] = 0.00800558188002041710;
        weight[2] = 0.00800558188002041710;
        weight[3] = 0.01594707683239050100;
        weight[4] = 0.01594707683239050100;
        weight[5] = 0.01594707683239050100;
        weight[6] = 0.01310914123079553000;
        weight[7] = 0.01310914123079553000;
        weight[8] = 0.01310914123079553000;
        weight[9] = 0.01958300096563562000;
        weight[10] = 0.01958300096563562000;
        weight[11] = 0.01958300096563562000;
        weight[12] = 0.01647088544153727000;
        weight[13] = 0.01647088544153727000;
        weight[14] = 0.01647088544153727000;
        weight[15] = 0.00854727907409210020;
        weight[16] = 0.00854727907409210020;
        weight[17] = 0.00854727907409210020;
        weight[18] = 0.00816188585722649180;
        weight[19] = 0.00816188585722649180;
        weight[20] = 0.00816188585722649180;
        weight[21] = 0.00612114653998377910;
        weight[22] = 0.00612114653998377910;
        weight[23] = 0.00612114653998377910;
        weight[24] = 0.00290849826493666490;
        weight[25] = 0.00290849826493666490;
        weight[26] = 0.00290849826493666490;
        weight[27] = 0.00069227524566199629;
        weight[28] = 0.00069227524566199629;
        weight[29] = 0.00069227524566199629;
        weight[30] = 0.00124828919927739700;
        weight[31] = 0.00124828919927739700;
        weight[32] = 0.00124828919927739700;
        weight[33] = 0.00124828919927739700;
        weight[34] = 0.00124828919927739700;
        weight[35] = 0.00124828919927739700;
        weight[36] = 0.00340475290880302200;
        weight[37] = 0.00340475290880302200;
        weight[38] = 0.00340475290880302200;
        weight[39] = 0.00340475290880302200;
        weight[40] = 0.00340475290880302200;
        weight[41] = 0.00340475290880302200;
        weight[42] = 0.00335965432606405090;
        weight[43] = 0.00335965432606405090;
        weight[44] = 0.00335965432606405090;
        weight[45] = 0.00335965432606405090;
        weight[46] = 0.00335965432606405090;
        weight[47] = 0.00335965432606405090;
        weight[48] = 0.00171615653949675410;
        weight[49] = 0.00171615653949675410;
        weight[50] = 0.00171615653949675410;
        weight[51] = 0.00171615653949675410;
        weight[52] = 0.00171615653949675410;
        weight[53] = 0.00171615653949675410;
        weight[54] = 0.00148085631671560600;
        weight[55] = 0.00148085631671560600;
        weight[56] = 0.00148085631671560600;
        weight[57] = 0.00148085631671560600;
        weight[58] = 0.00148085631671560600;
        weight[59] = 0.00148085631671560600;
        weight[60] = 0.00351131261072868500;
        weight[61] = 0.00351131261072868500;
        weight[62] = 0.00351131261072868500;
        weight[63] = 0.00351131261072868500;
        weight[64] = 0.00351131261072868500;
        weight[65] = 0.00351131261072868500;
        weight[66] = 0.00739355014970648380;
        weight[67] = 0.00739355014970648380;
        weight[68] = 0.00739355014970648380;
        weight[69] = 0.00739355014970648380;
        weight[70] = 0.00739355014970648380;
        weight[71] = 0.00739355014970648380;
        weight[72] = 0.00798308747737655820;
        weight[73] = 0.00798308747737655820;
        weight[74] = 0.00798308747737655820;
        weight[75] = 0.00798308747737655820;
        weight[76] = 0.00798308747737655820;
        weight[77] = 0.00798308747737655820;
        weight[78] = 0.00435596261315804140;
        weight[79] = 0.00435596261315804140;
        weight[80] = 0.00435596261315804140;
        weight[81] = 0.00435596261315804140;
        weight[82] = 0.00435596261315804140;
        weight[83] = 0.00435596261315804140;
        weight[84] = 0.00736505670141783180;
        weight[85] = 0.00736505670141783180;
        weight[86] = 0.00736505670141783180;
        weight[87] = 0.00736505670141783180;
        weight[88] = 0.00736505670141783180;
        weight[89] = 0.00736505670141783180;
        weight[90] = 0.01096357284641955000;
        weight[91] = 0.01096357284641955000;
        weight[92] = 0.01096357284641955000;
        weight[93] = 0.01096357284641955000;
        weight[94] = 0.01096357284641955000;
        weight[95] = 0.01096357284641955000;
        weight[96] = 0.01174996174354112100;
        weight[97] = 0.01174996174354112100;
        weight[98] = 0.01174996174354112100;
        weight[99] = 0.01174996174354112100;
        weight[100] = 0.01174996174354112100;
        weight[101] = 0.01174996174354112100;
        weight[102] = 0.01001560071379857000;
        weight[103] = 0.01001560071379857000;
        weight[104] = 0.01001560071379857000;
        weight[105] = 0.01001560071379857000;
        weight[106] = 0.01001560071379857000;
        weight[107] = 0.01001560071379857000;
        weight[108] = 0.01330964078762868000;
        weight[109] = 0.01330964078762868000;
        weight[110] = 0.01330964078762868000;
        weight[111] = 0.01330964078762868000;
        weight[112] = 0.01330964078762868000;
        weight[113] = 0.01330964078762868000;
        weight[114] = 0.01415444650522614000;
        weight[115] = 0.01415444650522614000;
        weight[116] = 0.01415444650522614000;
        weight[117] = 0.01415444650522614000;
        weight[118] = 0.01415444650522614000;
        weight[119] = 0.01415444650522614000;
        weight[120] = 0.01488137956116801000;
        weight[121] = 0.01488137956116801000;
        weight[122] = 0.01488137956116801000;
        weight[123] = 0.01488137956116801000;
        weight[124] = 0.01488137956116801000;
        weight[125] = 0.01488137956116801000;
    } else if (deg<=30) {
    	pnt.resize(175); weight.resize(175);
        pnt[0] = Point(0.33333333333332998000 , 0.33333333333332998000);
        pnt[1] = Point(0.00733011643276999980 , 0.49633494178361998000);
        pnt[2] = Point(0.49633494178361998000 , 0.49633494178361998000);
        pnt[3] = Point(0.49633494178361998000 , 0.00733011643276999980);
        pnt[4] = Point(0.08299567580295999500 , 0.45850216209852002000);
        pnt[5] = Point(0.45850216209852002000 , 0.45850216209852002000);
        pnt[6] = Point(0.45850216209852002000 , 0.08299567580295999500);
        pnt[7] = Point(0.15098095612540999000 , 0.42450952193729002000);
        pnt[8] = Point(0.42450952193729002000 , 0.42450952193729002000);
        pnt[9] = Point(0.42450952193729002000 , 0.15098095612540999000);
        pnt[10] = Point(0.23590585989217000000 , 0.38204707005392002000);
        pnt[11] = Point(0.38204707005392002000 , 0.38204707005392002000);
        pnt[12] = Point(0.38204707005392002000 , 0.23590585989217000000);
        pnt[13] = Point(0.43802430840785000000 , 0.28098784579607999000);
        pnt[14] = Point(0.28098784579607999000 , 0.28098784579607999000);
        pnt[15] = Point(0.28098784579607999000 , 0.43802430840785000000);
        pnt[16] = Point(0.54530204829192996000 , 0.22734897585402999000);
        pnt[17] = Point(0.22734897585402999000 , 0.22734897585402999000);
        pnt[18] = Point(0.22734897585402999000 , 0.54530204829192996000);
        pnt[19] = Point(0.65088177698254002000 , 0.17455911150872999000);
        pnt[20] = Point(0.17455911150872999000 , 0.17455911150872999000);
        pnt[21] = Point(0.17455911150872999000 , 0.65088177698254002000);
        pnt[22] = Point(0.75348314559713003000 , 0.12325842720143999000);
        pnt[23] = Point(0.12325842720143999000 , 0.12325842720143999000);
        pnt[24] = Point(0.12325842720143999000 , 0.75348314559713003000);
        pnt[25] = Point(0.83983154221560996000 , 0.08008422889220000200);
        pnt[26] = Point(0.08008422889220000200 , 0.08008422889220000200);
        pnt[27] = Point(0.08008422889220000200 , 0.83983154221560996000);
        pnt[28] = Point(0.90445106518420004000 , 0.04777446740790000000);
        pnt[29] = Point(0.04777446740790000000 , 0.04777446740790000000);
        pnt[30] = Point(0.04777446740790000000 , 0.90445106518420004000);
        pnt[31] = Point(0.95655897063971995000 , 0.02172051468014000000);
        pnt[32] = Point(0.02172051468014000000 , 0.02172051468014000000);
        pnt[33] = Point(0.02172051468014000000 , 0.95655897063971995000);
        pnt[34] = Point(0.99047064476913005000 , 0.00476467761544000030);
        pnt[35] = Point(0.00476467761544000030 , 0.00476467761544000030);
        pnt[36] = Point(0.00476467761544000030 , 0.99047064476913005000);
        pnt[37] = Point(0.00092537119334999999 , 0.41529527091330998000);
        pnt[38] = Point(0.41529527091330998000 , 0.58377935789334001000);
        pnt[39] = Point(0.58377935789334001000 , 0.00092537119334999999);
        pnt[40] = Point(0.41529527091330998000 , 0.00092537119334999999);
        pnt[41] = Point(0.58377935789334001000 , 0.41529527091330998000);
        pnt[42] = Point(0.00092537119334999999 , 0.58377935789334001000);
        pnt[43] = Point(0.00138592585556000010 , 0.06118990978535000100);
        pnt[44] = Point(0.06118990978535000100 , 0.93742416435909004000);
        pnt[45] = Point(0.93742416435909004000 , 0.00138592585556000010);
        pnt[46] = Point(0.06118990978535000100 , 0.00138592585556000010);
        pnt[47] = Point(0.93742416435909004000 , 0.06118990978535000100);
        pnt[48] = Point(0.00138592585556000010 , 0.93742416435909004000);
        pnt[49] = Point(0.00368241545590999990 , 0.16490869013691001000);
        pnt[50] = Point(0.16490869013691001000 , 0.83140889440718002000);
        pnt[51] = Point(0.83140889440718002000 , 0.00368241545590999990);
        pnt[52] = Point(0.16490869013691001000 , 0.00368241545590999990);
        pnt[53] = Point(0.83140889440718002000 , 0.16490869013691001000);
        pnt[54] = Point(0.00368241545590999990 , 0.83140889440718002000);
        pnt[55] = Point(0.00390322342416000000 , 0.02503506223199999900);
        pnt[56] = Point(0.02503506223199999900 , 0.97106171434384003000);
        pnt[57] = Point(0.97106171434384003000 , 0.00390322342416000000);
        pnt[58] = Point(0.02503506223199999900 , 0.00390322342416000000);
        pnt[59] = Point(0.97106171434384003000 , 0.02503506223199999900);
        pnt[60] = Point(0.00390322342416000000 , 0.97106171434384003000);
        pnt[61] = Point(0.00323324815500999980 , 0.30606446515109997000);
        pnt[62] = Point(0.30606446515109997000 , 0.69070228669389000000);
        pnt[63] = Point(0.69070228669389000000 , 0.00323324815500999980);
        pnt[64] = Point(0.30606446515109997000 , 0.00323324815500999980);
        pnt[65] = Point(0.69070228669389000000 , 0.30606446515109997000);
        pnt[66] = Point(0.00323324815500999980 , 0.69070228669389000000);
        pnt[67] = Point(0.00646743211223999980 , 0.10707328373022000000);
        pnt[68] = Point(0.10707328373022000000 , 0.88645928415754005000);
        pnt[69] = Point(0.88645928415754005000 , 0.00646743211223999980);
        pnt[70] = Point(0.10707328373022000000 , 0.00646743211223999980);
        pnt[71] = Point(0.88645928415754005000 , 0.10707328373022000000);
        pnt[72] = Point(0.00646743211223999980 , 0.88645928415754005000);
        pnt[73] = Point(0.00324747549132999980 , 0.22995754934557999000);
        pnt[74] = Point(0.22995754934557999000 , 0.76679497516308004000);
        pnt[75] = Point(0.76679497516308004000 , 0.00324747549132999980);
        pnt[76] = Point(0.22995754934557999000 , 0.00324747549132999980);
        pnt[77] = Point(0.76679497516308004000 , 0.22995754934557999000);
        pnt[78] = Point(0.00324747549132999980 , 0.76679497516308004000);
        pnt[79] = Point(0.00867509080675000000 , 0.33703663330577999000);
        pnt[80] = Point(0.33703663330577999000 , 0.65428827588745997000);
        pnt[81] = Point(0.65428827588745997000 , 0.00867509080675000000);
        pnt[82] = Point(0.33703663330577999000 , 0.00867509080675000000);
        pnt[83] = Point(0.65428827588745997000 , 0.33703663330577999000);
        pnt[84] = Point(0.00867509080675000000 , 0.65428827588745997000);
        pnt[85] = Point(0.01559702646731000000 , 0.05625657618206000200);
        pnt[86] = Point(0.05625657618206000200 , 0.92814639735062998000);
        pnt[87] = Point(0.92814639735062998000 , 0.01559702646731000000);
        pnt[88] = Point(0.05625657618206000200 , 0.01559702646731000000);
        pnt[89] = Point(0.92814639735062998000 , 0.05625657618206000200);
        pnt[90] = Point(0.01559702646731000000 , 0.92814639735062998000);
        pnt[91] = Point(0.01797672125369000100 , 0.40245137521239999000);
        pnt[92] = Point(0.40245137521239999000 , 0.57957190353390997000);
        pnt[93] = Point(0.57957190353390997000 , 0.01797672125369000100);
        pnt[94] = Point(0.40245137521239999000 , 0.01797672125369000100);
        pnt[95] = Point(0.57957190353390997000 , 0.40245137521239999000);
        pnt[96] = Point(0.01797672125369000100 , 0.57957190353390997000);
        pnt[97] = Point(0.01712424535389000000 , 0.24365470201083000000);
        pnt[98] = Point(0.24365470201083000000 , 0.73922105263528004000);
        pnt[99] = Point(0.73922105263528004000 , 0.01712424535389000000);
        pnt[100] = Point(0.24365470201083000000 , 0.01712424535389000000);
        pnt[101] = Point(0.73922105263528004000 , 0.24365470201083000000);
        pnt[102] = Point(0.01712424535389000000 , 0.73922105263528004000);
        pnt[103] = Point(0.02288340534658000000 , 0.16538958561452999000);
        pnt[104] = Point(0.16538958561452999000 , 0.81172700903887995000);
        pnt[105] = Point(0.81172700903887995000 , 0.02288340534658000000);
        pnt[106] = Point(0.16538958561452999000 , 0.02288340534658000000);
        pnt[107] = Point(0.81172700903887995000 , 0.16538958561452999000);
        pnt[108] = Point(0.02288340534658000000 , 0.81172700903887995000);
        pnt[109] = Point(0.03273759728777000200 , 0.09930187449584999800);
        pnt[110] = Point(0.09930187449584999800 , 0.86796052821639003000);
        pnt[111] = Point(0.86796052821639003000 , 0.03273759728777000200);
        pnt[112] = Point(0.09930187449584999800 , 0.03273759728777000200);
        pnt[113] = Point(0.86796052821639003000 , 0.09930187449584999800);
        pnt[114] = Point(0.03273759728777000200 , 0.86796052821639003000);
        pnt[115] = Point(0.03382101234234000100 , 0.30847833306904998000);
        pnt[116] = Point(0.30847833306904998000 , 0.65770065458860005000);
        pnt[117] = Point(0.65770065458860005000 , 0.03382101234234000100);
        pnt[118] = Point(0.30847833306904998000 , 0.03382101234234000100);
        pnt[119] = Point(0.65770065458860005000 , 0.30847833306904998000);
        pnt[120] = Point(0.03382101234234000100 , 0.65770065458860005000);
        pnt[121] = Point(0.03554761446001999900 , 0.46066831859210999000);
        pnt[122] = Point(0.46066831859210999000 , 0.50378406694787004000);
        pnt[123] = Point(0.50378406694787004000 , 0.03554761446001999900);
        pnt[124] = Point(0.46066831859210999000 , 0.03554761446001999900);
        pnt[125] = Point(0.50378406694787004000 , 0.46066831859210999000);
        pnt[126] = Point(0.03554761446001999900 , 0.50378406694787004000);
        pnt[127] = Point(0.05053979030687000300 , 0.21881529945393000000);
        pnt[128] = Point(0.21881529945393000000 , 0.73064491023919997000);
        pnt[129] = Point(0.73064491023919997000 , 0.05053979030687000300);
        pnt[130] = Point(0.21881529945393000000 , 0.05053979030687000300);
        pnt[131] = Point(0.73064491023919997000 , 0.21881529945393000000);
        pnt[132] = Point(0.05053979030687000300 , 0.73064491023919997000);
        pnt[133] = Point(0.05701471491573000000 , 0.37920955156026998000);
        pnt[134] = Point(0.37920955156026998000 , 0.56377573352399002000);
        pnt[135] = Point(0.56377573352399002000 , 0.05701471491573000000);
        pnt[136] = Point(0.37920955156026998000 , 0.05701471491573000000);
        pnt[137] = Point(0.56377573352399002000 , 0.37920955156026998000);
        pnt[138] = Point(0.05701471491573000000 , 0.56377573352399002000);
        pnt[139] = Point(0.06415280642119999800 , 0.14296081941819000000);
        pnt[140] = Point(0.14296081941819000000 , 0.79288637416061003000);
        pnt[141] = Point(0.79288637416061003000 , 0.06415280642119999800);
        pnt[142] = Point(0.14296081941819000000 , 0.06415280642119999800);
        pnt[143] = Point(0.79288637416061003000 , 0.14296081941819000000);
        pnt[144] = Point(0.06415280642119999800 , 0.79288637416061003000);
        pnt[145] = Point(0.08050114828762999800 , 0.28373128210592002000);
        pnt[146] = Point(0.28373128210592002000 , 0.63576756960644998000);
        pnt[147] = Point(0.63576756960644998000 , 0.08050114828762999800);
        pnt[148] = Point(0.28373128210592002000 , 0.08050114828762999800);
        pnt[149] = Point(0.63576756960644998000 , 0.28373128210592002000);
        pnt[150] = Point(0.08050114828762999800 , 0.63576756960644998000);
        pnt[151] = Point(0.10436706813453001000 , 0.19673744100443999000);
        pnt[152] = Point(0.19673744100443999000 , 0.69889549086102998000);
        pnt[153] = Point(0.69889549086102998000 , 0.10436706813453001000);
        pnt[154] = Point(0.19673744100443999000 , 0.10436706813453001000);
        pnt[155] = Point(0.69889549086102998000 , 0.19673744100443999000);
        pnt[156] = Point(0.10436706813453001000 , 0.69889549086102998000);
        pnt[157] = Point(0.11384489442875000000 , 0.35588914121165999000);
        pnt[158] = Point(0.35588914121165999000 , 0.53026596435958995000);
        pnt[159] = Point(0.53026596435958995000 , 0.11384489442875000000);
        pnt[160] = Point(0.35588914121165999000 , 0.11384489442875000000);
        pnt[161] = Point(0.53026596435958995000 , 0.35588914121165999000);
        pnt[162] = Point(0.11384489442875000000 , 0.53026596435958995000);
        pnt[163] = Point(0.14536348771551999000 , 0.25981868535190999000);
        pnt[164] = Point(0.25981868535190999000 , 0.59481782693256002000);
        pnt[165] = Point(0.59481782693256002000 , 0.14536348771551999000);
        pnt[166] = Point(0.25981868535190999000 , 0.14536348771551999000);
        pnt[167] = Point(0.59481782693256002000 , 0.25981868535190999000);
        pnt[168] = Point(0.14536348771551999000 , 0.59481782693256002000);
        pnt[169] = Point(0.18994565282198000000 , 0.32192318123129998000);
        pnt[170] = Point(0.32192318123129998000 , 0.48813116594672001000);
        pnt[171] = Point(0.48813116594672001000 , 0.18994565282198000000);
        pnt[172] = Point(0.32192318123129998000 , 0.18994565282198000000);
        pnt[173] = Point(0.48813116594672001000 , 0.32192318123129998000);
        pnt[174] = Point(0.18994565282198000000 , 0.48813116594672001000);

        weight[0] = 0.01557996020289920000;
        weight[1] = 0.00317723370053413400;
        weight[2] = 0.00317723370053413400;
        weight[3] = 0.00317723370053413400;
        weight[4] = 0.01048342663573077100;
        weight[5] = 0.01048342663573077100;
        weight[6] = 0.01048342663573077100;
        weight[7] = 0.01320945957774363000;
        weight[8] = 0.01320945957774363000;
        weight[9] = 0.01320945957774363000;
        weight[10] = 0.01497500696627149900;
        weight[11] = 0.01497500696627149900;
        weight[12] = 0.01497500696627149900;
        weight[13] = 0.01498790444338419000;
        weight[14] = 0.01498790444338419000;
        weight[15] = 0.01498790444338419000;
        weight[16] = 0.01333886474102166000;
        weight[17] = 0.01333886474102166000;
        weight[18] = 0.01333886474102166000;
        weight[19] = 0.01088917111390201100;
        weight[20] = 0.01088917111390201100;
        weight[21] = 0.01088917111390201100;
        weight[22] = 0.00818944066089346070;
        weight[23] = 0.00818944066089346070;
        weight[24] = 0.00818944066089346070;
        weight[25] = 0.00557538758860778510;
        weight[26] = 0.00557538758860778510;
        weight[27] = 0.00557538758860778510;
        weight[28] = 0.00319121647341197600;
        weight[29] = 0.00319121647341197600;
        weight[30] = 0.00319121647341197600;
        weight[31] = 0.00129671514432704500;
        weight[32] = 0.00129671514432704500;
        weight[33] = 0.00129671514432704500;
        weight[34] = 0.00029826282613491719;
        weight[35] = 0.00029826282613491719;
        weight[36] = 0.00029826282613491719;
        weight[37] = 0.00099890568507889641;
        weight[38] = 0.00099890568507889641;
        weight[39] = 0.00099890568507889641;
        weight[40] = 0.00099890568507889641;
        weight[41] = 0.00099890568507889641;
        weight[42] = 0.00099890568507889641;
        weight[43] = 0.00046285084917325331;
        weight[44] = 0.00046285084917325331;
        weight[45] = 0.00046285084917325331;
        weight[46] = 0.00046285084917325331;
        weight[47] = 0.00046285084917325331;
        weight[48] = 0.00046285084917325331;
        weight[49] = 0.00123445133638241290;
        weight[50] = 0.00123445133638241290;
        weight[51] = 0.00123445133638241290;
        weight[52] = 0.00123445133638241290;
        weight[53] = 0.00123445133638241290;
        weight[54] = 0.00123445133638241290;
        weight[55] = 0.00057071985224320615;
        weight[56] = 0.00057071985224320615;
        weight[57] = 0.00057071985224320615;
        weight[58] = 0.00057071985224320615;
        weight[59] = 0.00057071985224320615;
        weight[60] = 0.00057071985224320615;
        weight[61] = 0.00112694612587762410;
        weight[62] = 0.00112694612587762410;
        weight[63] = 0.00112694612587762410;
        weight[64] = 0.00112694612587762410;
        weight[65] = 0.00112694612587762410;
        weight[66] = 0.00112694612587762410;
        weight[67] = 0.00174786694940733710;
        weight[68] = 0.00174786694940733710;
        weight[69] = 0.00174786694940733710;
        weight[70] = 0.00174786694940733710;
        weight[71] = 0.00174786694940733710;
        weight[72] = 0.00174786694940733710;
        weight[73] = 0.00118281881503165690;
        weight[74] = 0.00118281881503165690;
        weight[75] = 0.00118281881503165690;
        weight[76] = 0.00118281881503165690;
        weight[77] = 0.00118281881503165690;
        weight[78] = 0.00118281881503165690;
        weight[79] = 0.00199083929467503380;
        weight[80] = 0.00199083929467503380;
        weight[81] = 0.00199083929467503380;
        weight[82] = 0.00199083929467503380;
        weight[83] = 0.00199083929467503380;
        weight[84] = 0.00199083929467503380;
        weight[85] = 0.00190041279503598000;
        weight[86] = 0.00190041279503598000;
        weight[87] = 0.00190041279503598000;
        weight[88] = 0.00190041279503598000;
        weight[89] = 0.00190041279503598000;
        weight[90] = 0.00190041279503598000;
        weight[91] = 0.00449836580881745120;
        weight[92] = 0.00449836580881745120;
        weight[93] = 0.00449836580881745120;
        weight[94] = 0.00449836580881745120;
        weight[95] = 0.00449836580881745120;
        weight[96] = 0.00449836580881745120;
        weight[97] = 0.00347871946027471890;
        weight[98] = 0.00347871946027471890;
        weight[99] = 0.00347871946027471890;
        weight[100] = 0.00347871946027471890;
        weight[101] = 0.00347871946027471890;
        weight[102] = 0.00347871946027471890;
        weight[103] = 0.00410239903672395340;
        weight[104] = 0.00410239903672395340;
        weight[105] = 0.00410239903672395340;
        weight[106] = 0.00410239903672395340;
        weight[107] = 0.00410239903672395340;
        weight[108] = 0.00410239903672395340;
        weight[109] = 0.00402176154974416210;
        weight[110] = 0.00402176154974416210;
        weight[111] = 0.00402176154974416210;
        weight[112] = 0.00402176154974416210;
        weight[113] = 0.00402176154974416210;
        weight[114] = 0.00402176154974416210;
        weight[115] = 0.00603316466079506590;
        weight[116] = 0.00603316466079506590;
        weight[117] = 0.00603316466079506590;
        weight[118] = 0.00603316466079506590;
        weight[119] = 0.00603316466079506590;
        weight[120] = 0.00603316466079506590;
        weight[121] = 0.00394629030212959810;
        weight[122] = 0.00394629030212959810;
        weight[123] = 0.00394629030212959810;
        weight[124] = 0.00394629030212959810;
        weight[125] = 0.00394629030212959810;
        weight[126] = 0.00394629030212959810;
        weight[127] = 0.00664404453768026840;
        weight[128] = 0.00664404453768026840;
        weight[129] = 0.00664404453768026840;
        weight[130] = 0.00664404453768026840;
        weight[131] = 0.00664404453768026840;
        weight[132] = 0.00664404453768026840;
        weight[133] = 0.00825430585607845810;
        weight[134] = 0.00825430585607845810;
        weight[135] = 0.00825430585607845810;
        weight[136] = 0.00825430585607845810;
        weight[137] = 0.00825430585607845810;
        weight[138] = 0.00825430585607845810;
        weight[139] = 0.00649605663340641070;
        weight[140] = 0.00649605663340641070;
        weight[141] = 0.00649605663340641070;
        weight[142] = 0.00649605663340641070;
        weight[143] = 0.00649605663340641070;
        weight[144] = 0.00649605663340641070;
        weight[145] = 0.00925277814414660230;
        weight[146] = 0.00925277814414660230;
        weight[147] = 0.00925277814414660230;
        weight[148] = 0.00925277814414660230;
        weight[149] = 0.00925277814414660230;
        weight[150] = 0.00925277814414660230;
        weight[151] = 0.00916492072629427990;
        weight[152] = 0.00916492072629427990;
        weight[153] = 0.00916492072629427990;
        weight[154] = 0.00916492072629427990;
        weight[155] = 0.00916492072629427990;
        weight[156] = 0.00916492072629427990;
        weight[157] = 0.01156952462809767100;
        weight[158] = 0.01156952462809767100;
        weight[159] = 0.01156952462809767100;
        weight[160] = 0.01156952462809767100;
        weight[161] = 0.01156952462809767100;
        weight[162] = 0.01156952462809767100;
        weight[163] = 0.01176111646760917000;
        weight[164] = 0.01176111646760917000;
        weight[165] = 0.01176111646760917000;
        weight[166] = 0.01176111646760917000;
        weight[167] = 0.01176111646760917000;
        weight[168] = 0.01176111646760917000;
        weight[169] = 0.01382470218216540000;
        weight[170] = 0.01382470218216540000;
        weight[171] = 0.01382470218216540000;
        weight[172] = 0.01382470218216540000;
        weight[173] = 0.01382470218216540000;
        weight[174] = 0.01382470218216540000;
    }
    else {
    	// throw error , we do not have have the required order
    	TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error  , " GaussLobattoQuadrature::getTriangleQuadPoints order of the quadrature to high !!! ");
    }
}
