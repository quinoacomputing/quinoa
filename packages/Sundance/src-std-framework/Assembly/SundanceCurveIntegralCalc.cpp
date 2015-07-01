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
 * SundanceCurveIntagralCalc.cpp
 *
 *  Created on: Sep 2, 2010
 *      Author: benk
 */

#include "SundanceCurveIntegralCalc.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundancePolygon2D.hpp"
#include "SundancePolygonQuadrature.hpp"
#include "SundanceSurf3DCalc.hpp"
#include "SundanceSurfQuadrature.hpp"

using namespace Sundance;

CurveIntegralCalc::CurveIntegralCalc() {
	//nothing to do
}

void CurveIntegralCalc::getCurveQuadPoints(CellType  maxCellType ,
							   int maxCellLID ,
							   const Mesh& mesh ,
							   const ParametrizedCurve& paramCurve,
							   const QuadratureFamily& quad ,
							   Array<Point>& curvePoints ,
							   Array<Point>& curveDerivs ,
							   Array<Point>& curveNormals ){

	int verb = 0;
	Tabs tabs;

    // get all the points from the maxDimCell
	int nr_point = mesh.numFacets( mesh.spatialDim() , 0 , 0);
	Array<Point> maxCellsPoints(nr_point);
	int tmp_i , point_LID;

	SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints nr points per cell: " << nr_point)
	for (int jj = 0 ; jj < nr_point ; jj++){
		point_LID = mesh.facetLID( mesh.spatialDim() , maxCellLID , 0 , jj , tmp_i );
		maxCellsPoints[jj] = mesh.nodePosition(point_LID);
		SUNDANCE_MSG3(verb, tabs << " max cell point p[" << jj << "]:"<< maxCellsPoints[jj]);
	}

	// if we have a polygon in 2D then use specialized curve integrals
	if (usePolygoneCurveQuadrature( maxCellType , mesh , paramCurve)){
		// specialized curve integral
		CurveIntegralCalc::getCurveQuadPoints_polygon( maxCellLID , maxCellType , mesh , maxCellsPoints , paramCurve,
								    quad , curvePoints , curveDerivs , curveNormals);
	} else {
		// general curve integral
		if ( use3DSurfQuadrature( maxCellType , mesh ) ){
			CurveIntegralCalc::getSurfQuadPoints(maxCellLID , maxCellType , mesh , maxCellsPoints , paramCurve,
					quad , curvePoints , curveDerivs , curveNormals);
			 //CurveIntegralCalc::getCurveQuadPoints_line( maxCellLID , maxCellType , mesh , maxCellsPoints , paramCurve,
				//								quad , curvePoints , curveDerivs , curveNormals);
		} else {
		    CurveIntegralCalc::getCurveQuadPoints_line( maxCellLID , maxCellType , mesh , maxCellsPoints , paramCurve,
									quad , curvePoints , curveDerivs , curveNormals);
		}
	}

}

void CurveIntegralCalc::getCurveQuadPoints_line(
                                   int  maxCellLID ,
                                   CellType  maxCellType ,
                                   const Mesh& mesh ,
                                   const Array<Point>& cellPoints,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals ){

	int verb = 0;
	Tabs tabs;

    // select the dimension of the curve
	switch (paramCurve.getCurveDim()){
	  case 1: {

		  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_line, curve has ONE dimnesion")
		  // 1D curve integral in 2D or 3D

		  // first determine the two intersection points with the edges
		  Point startPoint(0.0,0.0);
		  Point endPoint(0.0,0.0);
	      int nrPoints , total_points = 0 ;

		  switch (maxCellType){
		    case QuadCell:{
		    	// loop over each edge and detect intersection point
		    	// there can be only one
		    	TEUCHOS_TEST_FOR_EXCEPTION( cellPoints.size() != 4 ,
		    			std::runtime_error ," CurveIntegralCalc::getCurveQuadPoints , QuadCell must have 4 points" );
				SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_line, on QuadCell")
		    	Array<Point> result(0);
		    	int edegIndex[4][2] = { {0,1} , {0,2} , {1,3} , {2,3} };

		    	// loop over the edges
		    	for (int jj = 0 ; jj < 4 ; jj++ ){
					paramCurve.returnIntersectPoints(cellPoints[edegIndex[jj][0]], cellPoints[edegIndex[jj][1]], nrPoints , result);
					// test if we have intersection point
			    	if (nrPoints == 1){
			    		if (total_points == 0) startPoint = result[0];
			    		else endPoint = result[0];
			    		SUNDANCE_MSG3(verb, tabs << "found Int. point:" << result[0]);
						total_points += nrPoints;
			    	}
					SUNDANCE_MSG3(verb, tabs << "ind:" << jj << ", nr Int points :" << nrPoints
							<< " , start:" << startPoint << ", end:"<< endPoint);
					// if we have more than one intersection point then just ignore that edge
					TEUCHOS_TEST_FOR_EXCEPTION( nrPoints > 1 ,
							std::runtime_error , " CurveIntegralCalc::getCurveQuadPoints , QuadCell one edge " << jj << " , can have only one intersection point"
							<< " edgeP0:" << cellPoints[edegIndex[jj][0]] << " edgeP1:" << cellPoints[edegIndex[jj][1]] );
		    	}
		    	// test if we have too much intersection points
		    	TEUCHOS_TEST_FOR_EXCEPTION( total_points > 2 ,std::runtime_error , " CurveIntegralCalc::getCurveQuadPoints total_points > 2 : " << total_points );
				SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_line, end finding intersection points")
		    } break;
		    case TriangleCell:{
		    	TEUCHOS_TEST_FOR_EXCEPTION( cellPoints.size() != 3 , std::runtime_error , " CurveIntegralCalc::getCurveQuadPoints , TriangleCell must have 3 points" );
		    	Array<Point> result;
		    	int edegIndex[3][2] = { {0,1} , {0,2} , {1,2} };

		    	// loop over the edges
		    	for (int jj = 0 ; jj < 3 ; jj++ ){
					paramCurve.returnIntersectPoints(cellPoints[edegIndex[jj][0]], cellPoints[edegIndex[jj][1]], nrPoints , result);
					// test if we have intersection point
			    	if (nrPoints == 1){
			    		if (total_points == 0) startPoint = result[0];
			    		else endPoint = result[0];
			    		SUNDANCE_MSG3(verb, tabs << "found Int. point:" << result[0]);
						total_points += nrPoints;
			    	}
			    	// if we have more than one intersection point then just ignore that edge
			    	TEUCHOS_TEST_FOR_EXCEPTION( nrPoints > 1 ,
							std::runtime_error , " CurveIntegralCalc::getCurveQuadPoints , TriangleCell one edge " << jj << " , can have only one intersection point"
							<< " edgeP0:" << cellPoints[edegIndex[jj][0]] << " edgeP1:" << cellPoints[edegIndex[jj][1]] );
		    	}
		    	// test if we have too much intersection points
		    	TEUCHOS_TEST_FOR_EXCEPTION( total_points > 2 ,std::runtime_error , " CurveIntegralCalc::getCurveQuadPoints total_points > 2 : " << total_points );
		    } break;
		    default : {
		    	TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error , "CurveIntegralCalc::getCurveQuadPoints , Unknown Cell in 2D" );
		    }
		  }
		  // test for to many intersection points
		  TEUCHOS_TEST_FOR_EXCEPTION( total_points > 2 ,std::runtime_error , " CurveIntegralCalc::getCurveQuadPoints , no much intersection points: " << total_points);

		  // -> having the intersection points now we have to determine the:
		  // quad points and the gradients at the points(which norm will be in the integration)
		  // and the normalized normal vector (normalized since we need the sin and cos for Nitsche )
		  // -> as a very simple approach we consider the curve as a line between intersection points "startPoint -> endPoint"

		  // in the X direction the line should always have increasing values
		  if (startPoint[0] > endPoint[0]){
			  Point tmp = startPoint;
			  startPoint = endPoint;
			  endPoint = tmp;
		  }
		  SUNDANCE_MSG3(verb, tabs << "start end and points , start:" << startPoint << ", end:"<< endPoint)

		  // get the quadrature points for the line
		  Array<Point> quadPoints;
		  Array<double> quadWeights;
		  quad.getPoints( LineCell , quadPoints , quadWeights );
		  int nr_quad_points = quadPoints.size();

		  // The intersection points , we distribute the points along the line
		  curvePoints.resize(nr_quad_points);
		  SUNDANCE_MSG3(verb, tabs << " setting reference quadrature points" );
		  for (int ii = 0 ; ii < nr_quad_points ; ii++) {
			  SUNDANCE_MSG3(verb, tabs << " nr:" << ii << " Line quad point:" << quadPoints[ii] << ", line:" << (endPoint - startPoint));
			  curvePoints[ii] = startPoint + quadPoints[ii][0]*(endPoint - startPoint);

			  // we transform the intersection points to the reference cell
			  // this should work for unstructured case
			  // - get the Jacobian of the cells
			  // - get the inverse of this Jacobian
			  // - apply this Jacobian to the physical points (pos0 = (0,0))of the cell

			  Array<int> cellLID(1); cellLID[0] = maxCellLID;
		      CellJacobianBatch JBatch;
		      JBatch.resize(1, 2, 2);
			  mesh.getJacobians(2, cellLID , JBatch);
			  double pointc[2] = { curvePoints[ii][0] - cellPoints[0][0] , curvePoints[ii][1] - cellPoints[0][1]};
			  JBatch.applyInvJ(0, pointc , 1 , false);
			  curvePoints[ii][0] = pointc[0];
			  curvePoints[ii][1] = pointc[1];

			  SUNDANCE_MSG3(verb, tabs << " quad point nr:" << ii << " = " << curvePoints[ii]);
		  }


		  // The derivatives at points, for this simple method are simple and constant over the whole quadrature
		  curveDerivs.resize(nr_quad_points);
		  Point dist_point(endPoint[0] - startPoint[0] , endPoint[1] - startPoint[1]);
		  SUNDANCE_MSG3(verb, tabs << " setting derivative values points" );
		  for (int ii = 0 ; ii < nr_quad_points ; ii++) {
			  curveDerivs[ii] = endPoint - startPoint;
			  SUNDANCE_MSG3(verb, tabs << " curve Derivatives point nr:" << ii << " = " << curveDerivs[ii]);
		  }


		  // calculate the norms to the curve
		  curveNormals.resize(nr_quad_points);

		  // this is a shorter variant
		  SUNDANCE_MSG3(verb, tabs << " setting normal values at points" );
		  for (int ii = 0 ; ii < nr_quad_points ; ii++) {
			  Point norm_vec = Point(-curveDerivs[ii][1], curveDerivs[ii][0]);
			  norm_vec = norm_vec / ::sqrt(norm_vec*norm_vec);
			  double relative_length = 1e-2 * ::sqrt((endPoint - startPoint)*(endPoint - startPoint));
			  if ( paramCurve.curveEquation( startPoint + relative_length*norm_vec ) > 0 ) {
				  curveNormals[ii] = -norm_vec;
			  }
			  else {
				  curveNormals[ii] = norm_vec;
			  }
			  SUNDANCE_MSG3(verb, tabs << " curve Normals point nr:" << ii << " = " << curveNormals[ii]);
		  }


	  } break;

//========================== END 1D curve in 2D context ==========================



	  case 2: {
// ====================2D curve integral in 3D context =============================

		  // take the intersection with each edge,
		  // from these intersection points construct one 3D (convex) polygon (from triangle to hexadron)
		  //  - maybe divide the polygons into many triangles
		  // map the quadrature points somehow inside this surface

		  // here we consider a simple surface, as it is in bilinear interpolation
		  Tabs tabs;
		  switch (maxCellType){
		    case BrickCell:{
		    	int edegIndex[12][2] = { {0,1} , {0,2} , {0,4} , {1,3} , {1,5} , {2,3} , {2,6} , {3,7} ,
		    			                 {4,5} , {4,6} , {5,7} , {6,7} };
		    	int faceEdges[6][4]  = { {0,1,3,5} , {0,2,4,8} , {1,2,6,9} , {3,4,7,10} , {5,6,7,11} , {8,9,10,11} };
		    	// loop over the edges
		    	int t_inters_points = 0;
		    	int nrPoints = 0;
		    	Array<Point> edgeIntersectPoint(t_inters_points);
		    	Array<int> edgeIndex(t_inters_points);
		    	Array<int> edgePointI(12,-1);

		    	// loop over the edges and get the intersection points
		    	for (int jj = 0 ; jj < 12 ; jj++ ){
			    	Array<Point> result(0);
					paramCurve.returnIntersectPoints(cellPoints[edegIndex[jj][0]], cellPoints[edegIndex[jj][1]], nrPoints , result);
					// test if we have intersection point
			    	if (nrPoints > 0){
			    		edgeIntersectPoint.resize( t_inters_points + 1 );
			    		edgeIndex.resize( t_inters_points + 1 );
			    		edgeIntersectPoint[t_inters_points] = result[0];
			    		edgeIndex[t_inters_points] = jj;
			    		edgePointI[jj] = t_inters_points;
			    		SUNDANCE_MSG3(verb, tabs << "found Int. point [" << t_inters_points <<"]=" << result[0]);
			    		t_inters_points++;
			    	}
					SUNDANCE_MSG3(verb, tabs << "ind:" << jj << ", nr Int points :" << nrPoints );
					TEUCHOS_TEST_FOR_EXCEPTION( nrPoints > 1 ,
							std::runtime_error , " CurveIntegralCalc::getCurveQuadPoints , QuadCell one edge " << jj << " , can have only one intersection point" );
		    	} // getting intersection points


		    	// -------- get the point which is nearest to the other points ---------
		    	Array<double> totDist(t_inters_points);
		    	// calculate the distance from each point to each point
		    	SUNDANCE_MSG3(verb, tabs << " Calculating distances ");
		    	for (int i = 0 ; i < t_inters_points ; i++){
		    		double sum = 0.0;
		    		for (int j = 0 ; j < t_inters_points ; j++){
		    			sum = sum + ::sqrt(edgeIntersectPoint[i]*edgeIntersectPoint[j]);
		    		}
		    		totDist[i] = sum;
		    	}
		    	// find the point which is mostly in the middle
		    	int minIndex = 0;
		    	for (int i = 1 ; i < t_inters_points ; i++){
		    		if (totDist[i] < totDist[minIndex]) minIndex = i;
		    	}
		    	SUNDANCE_MSG3(verb, tabs << " minIndex = " << minIndex );


		    	// --------- discover the polygon and the triangles-----------
		    	// todo: what if there is no or only too few intersection points
		    	int intersTriangleInd = 0;
		    	int nrTriangles = t_inters_points-2;
		    	Array<int> triangles(3*nrTriangles,-1);
		    	Array<double> triangle_area(nrTriangles,0.0);
		    	double total_area = 0.0;
		    	SUNDANCE_MSG3(verb, tabs << " nrTriangles = " << nrTriangles );
		    	int intPointProc = -1;
		    	int intPointProc_lastNeigh = -1;

// -------------------------------------
		    	// get one neighboring intersection point of the "minIndex", this will be the starting direction of the polygon
		    	for (int jj = 0 ; jj < 6 ; jj++ )
		    	{
		    		int nrIpointPerFace = 0;
		    		int firstPI = -1 , secondPI = -1;
		    		// loop over each edge in one face
		    		for (int ii = 0 ; ii < 4 ; ii++){
		    			// if this edge has one intersection point
		    			if ( edgePointI[faceEdges[jj][ii]] >= 0){
		    				if (nrIpointPerFace > 0){
		    					secondPI = edgePointI[faceEdges[jj][ii]];
		    				}else{
		    					firstPI = edgePointI[faceEdges[jj][ii]];
		    				}
		    				TEUCHOS_TEST_FOR_EXCEPTION( nrIpointPerFace > 2 , std::runtime_error,"getCurveQuadPoints , nrIpointPerFace > 2 , " << nrIpointPerFace );
		    				nrIpointPerFace++;
		    			}
		    		}
		    		TEUCHOS_TEST_FOR_EXCEPTION( nrIpointPerFace == 1 , std::runtime_error,"getCurveQuadPoints , nrIpointPerFace == 1 , " << nrIpointPerFace );
		    		// Found the first neighboring point of "minIndex"
		    		if (( nrIpointPerFace > 1) && ((minIndex == firstPI) || (minIndex == secondPI) )){
			    			SUNDANCE_MSG3(verb, tabs << " Found  starting line:" << firstPI << " -> " << secondPI);
			    			intPointProc_lastNeigh = minIndex;
			    			if (minIndex == firstPI){
			    				intPointProc = secondPI;
			    			}else{
			    				intPointProc = firstPI;
			    			}
			    			// once we found then break the current for loop
			    			break;
		    		}
		    	}
		    	TEUCHOS_TEST_FOR_EXCEPTION( intPointProc < 0 , std::runtime_error,"getCurveQuadPoints , intPointProc < 0 , " << intPointProc );
		    	// store this as act. point and get the next one(which is not neighbor to "minIndex") which is not the one which we already had
		        // here we create all triangles, by traversing the polygon in one direction, so that the mapping will be continous
		        for (int pI = 0 ; pI < nrTriangles; pI++)
		        {
		        	// for each new
			    	for (int jj = 0 ; jj < 6 ; jj++ )
			    	{
			    		int nrIpointPerFace = 0;
			    		int firstPI = -1 , secondPI = -1;
			    		// loop over each edge in one face
			    		for (int ii = 0 ; ii < 4 ; ii++){
			    			// if this edge has one intersection point
			    			if ( edgePointI[faceEdges[jj][ii]] >= 0){
			    				if (nrIpointPerFace > 0){
			    					secondPI = edgePointI[faceEdges[jj][ii]];
			    				}else{
			    					firstPI = edgePointI[faceEdges[jj][ii]];
			    				}
			    				TEUCHOS_TEST_FOR_EXCEPTION( nrIpointPerFace > 2 , std::runtime_error,"getCurveQuadPoints , nrIpointPerFace > 2 , " << nrIpointPerFace );
			    				nrIpointPerFace++;
			    			}
			    		}
			    		TEUCHOS_TEST_FOR_EXCEPTION( nrIpointPerFace == 1 , std::runtime_error,"getCurveQuadPoints , nrIpointPerFace == 1 , " << nrIpointPerFace );
			    		// condition to find the neighbor of "intPointProc" but not the one which has been already found
			    		if (( nrIpointPerFace > 1) && ((intPointProc == firstPI) || (intPointProc == secondPI))
			    				&& ((intPointProc_lastNeigh != firstPI) && (intPointProc_lastNeigh != secondPI)))
			    		{
				    			SUNDANCE_MSG3(verb, tabs << " Found next line:" << firstPI << " -> " << secondPI);
				    			if (intPointProc == firstPI){
				    				intPointProc = secondPI;
					    			intPointProc_lastNeigh = firstPI;
				    			}else{
				    				intPointProc = firstPI;
					    			intPointProc_lastNeigh = secondPI;
				    			}
				    			// add triangle
			    				triangles[intersTriangleInd] = minIndex;
			    				triangles[intersTriangleInd+1] = intPointProc_lastNeigh;
			    				triangles[intersTriangleInd+2] = intPointProc;
			    				SUNDANCE_MSG3(verb, tabs << " Found triangle:" << minIndex << " , " << firstPI << " , " << secondPI);
			    				Point v1 = edgeIntersectPoint[firstPI] - edgeIntersectPoint[minIndex];
			    				Point v2 = edgeIntersectPoint[secondPI] - edgeIntersectPoint[minIndex];
			    				Point areaV = cross(v1,v2);
			    				SUNDANCE_MSG3(verb, tabs << " TriangINdex = " << intersTriangleInd << " , area = " << ::sqrt(areaV*areaV));
			    				triangle_area[ intersTriangleInd/3 ] = ::sqrt(areaV*areaV) / (double)2.0 ; // div by 2 to get the are of the triangle
			    				total_area = total_area + triangle_area[ intersTriangleInd/3 ];
			    				intersTriangleInd = intersTriangleInd + 3;
				    			// once we found the next triangle then break the current for loop
				    			break;
			    		}
			    	}
			    }
// ---------------------------
		    	SUNDANCE_MSG3(verb, tabs << " END Found triangle  total_area=" << total_area);


		    	// if we have only 1 triangle (3 intersection points) then add to the longest edge one middle point
		    	// so that we will have at least 2 triangles
		    	if (t_inters_points < 4){
		    		// 0,1,2
		    		SUNDANCE_MSG3(verb, tabs << " Add one additional point to have at least 2 triangles ");
		    		int longEdI1 = -1 , longEdI2 = -1;
		    		for (int ii = 0 ; ii < t_inters_points ; ii++)
		    			if ( (ii != minIndex) ){
		    				longEdI1 = ii;
		    				break;
		    			}
		    		for (int ii = 0 ; ii < t_inters_points ; ii++)
		    			if ( (ii != minIndex) && (ii != longEdI1) ){
		    				longEdI2 = ii;
		    				break;
		    			}
		    		TEUCHOS_TEST_FOR_EXCEPTION( (longEdI1 < 0)||(longEdI2 < 0), std::runtime_error," (longEdI1 < 0)||(longEdI2 < 0):" << longEdI1 << "," <<longEdI2);

		    		//add the middle point
		    		edgeIntersectPoint.resize(t_inters_points+1);
		    		edgeIntersectPoint[t_inters_points] = edgeIntersectPoint[longEdI1]/2.0 + edgeIntersectPoint[longEdI2]/2.0;

		    		int new_p = t_inters_points;
		    		t_inters_points += 1;
		    		nrTriangles += 1;

		    		// we'll have two new triangles
    				triangles.resize(3*(nrTriangles));
    				triangles[0] = minIndex;
    				triangles[1] = longEdI1;
    				triangles[2] = new_p;
    				triangles[3] = minIndex;
    				triangles[4] = new_p;
    				triangles[5] = longEdI2;
                    // the are will be half of the original one
    				triangle_area[0] = triangle_area[0] / 2.0;
    				triangle_area.resize(2);
    				triangle_area[1] = triangle_area[0];
		    	}

		    	// now map the quadrature points to the real cell (triangles)
				Array<Point> quadPoints;
				Array<double> quadWeights;
				quad.getPoints( QuadCell , quadPoints , quadWeights );
				int nr_quad_points = quadPoints.size();

				// --------- the quadrature points per cell ------------
				curvePoints.resize(nr_quad_points);
				Array<int> triangleInd(nr_quad_points,-1);
				SUNDANCE_MSG3(verb, tabs << " setting reference quadrature points" );
				for (int ii = 0 ; ii < nr_quad_points ; ii++) {
					//SUNDANCE_MSG3(verb, tabs << " nr:" << ii << " Quad quad point:" << quadPoints[ii] << ", line:" << (endPoint - startPoint));
					int tInd = -1;
					Point tmp(0.0,0.0,0.0);
					get3DRealCoordsOnSurf( quadPoints[ii] , cellPoints , triangles ,
					        nrTriangles , edgeIntersectPoint, tInd , tmp );
					triangleInd[ii] = tInd;
					curvePoints[ii] = tmp;
					//set the real quad points
					SUNDANCE_MSG3(verb, tabs << " quad point nr:" << ii << " = " << curvePoints[ii]);
				}

				// ------- The derivatives at points, for this simple method are simple and constant over the whole quadrature -----
				curveDerivs.resize(nr_quad_points);
				SUNDANCE_MSG3(verb, tabs << " setting surface values points" );
				for (int ii = 0 ; ii < nr_quad_points ; ii++) {
					// set the vector according to formula (only the vector noem will be taken in account)
					int tInd = triangleInd[ii] ;
					curveDerivs[ii] = Point( nrTriangles*triangle_area[tInd] , 0.0, 0.0 );
					//curveDerivs[ii] = Point( total_area , 0.0, 0.0 ); //total_area
					SUNDANCE_MSG3(verb, tabs << " curve Derivatives point nr:" << ii << " = " << curveDerivs[ii]);
				}


				// calculate the norms to the curve
				curveNormals.resize(nr_quad_points);
				SUNDANCE_MSG3(verb, tabs << " setting normal values at points" );
				for (int ii = 0 ; ii < nr_quad_points ; ii++) {
					int tInd = triangleInd[ii];
					// get the normal to the triangle
				    Point norm_vec = cross( edgeIntersectPoint[triangles[3*tInd+1]] - edgeIntersectPoint[triangles[3*tInd]] ,
				    		edgeIntersectPoint[triangles[3*tInd+2]] - edgeIntersectPoint[triangles[3*tInd]]);
				    double relative_length = 1e-2 * ::sqrt(norm_vec*norm_vec);
				    norm_vec = norm_vec / ::sqrt(norm_vec*norm_vec);
				    if ( paramCurve.curveEquation( edgeIntersectPoint[triangles[3*tInd]] + relative_length*norm_vec ) > 0 ) {
					    curveNormals[ii] = -norm_vec;
				    }
				    else {
					    curveNormals[ii] = norm_vec;
				    }
				    SUNDANCE_MSG3(verb, tabs << " curve Normals point nr:" << ii << " = " << curveNormals[ii]);
				}

                break;
		    } // ----------- end brick cell --------------
		    case TetCell:{
		    	// todo: later when this will be needed
		    	TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error,"CurveIntagralCalc::getCurveQuadPoints , not implemented for " << maxCellType );
		    	break;
		    }
		    default:{
		    	TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error,"CurveIntagralCalc::getCurveQuadPoints , not implemented for " << maxCellType );
		    }
		  }
		  SUNDANCE_MSG3(verb, tabs << " END CurveIntagralCalc::getCurveQuadPoints");
	  } break;

//========================== END 2D curve in 3D context ==========================

	  default: {
        // throw exception
		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error,"CurveIntagralCalc::getCurveQuadPoints , curve dimension must be 1 or two ");
	  }
	}
}


void CurveIntegralCalc::get3DRealCoordsOnSurf(const Point &refP ,
		const Array<Point>& cellPoints,
        const Array<int> &triangles ,
        const int nrTriag ,
        const Array<Point> edgeIntersectPoint,
        int &triagIndex ,
        Point &realPoint){

	Tabs tabs;
	int verb = 0;
	realPoint[0] = 0.0; realPoint[1] = 0.0; realPoint[2] = 0.0;
	SUNDANCE_MSG3(verb, tabs << " start get3DRealCoordsOnSurf refP="<<refP );
	Point v1;
	Point v2;
	// these matrixes were calculated in octave
	double case1[2][4] = {{1 ,-1  ,0   ,1},{ 1 , 0 ,  -1  , 1}};
	double case2[3][4] = {{1 ,-2, 0, 2} , { 2 , -2, -1, 2},{ 1 , 0, -1, 1}};
	double case3[4][4] = {{1, -2,0 , 2} , { 2 , -2, -1,2},{ 2, -1,-2, 2},{ 2 ,0, -2,  1}};
	double *refInv = 0;
	switch (nrTriag){
		case 2:{
			if (refP[0] >= refP[1]){ // first triangle
				triagIndex = 0; v1 = Point(1.0,0.0); v2 = Point(1.0,1.0);
			}else{// second triangle
				triagIndex = 1; v1 = Point(1.0,1.0); v2 = Point(0.0,1.0);
			}
			refInv = case1[triagIndex];
			break;
		}
		case 3:{
			if (refP[0] >= 2*refP[1]){ // first triangle
				triagIndex = 0; v1 = Point(1.0,0.0); v2 = Point(1.0,0.5);
			}else if ((refP[0] <= 2*refP[1]) && (refP[0] >= refP[1])){ // second triangle
				triagIndex = 1; v1 = Point(1.0,0.5); v2 = Point(1.0,1.0);
			}else{ // third triangle
				triagIndex = 2; v1 = Point(1.0,1.0); v2 = Point(0.0,1.0);
			}
			refInv = case2[triagIndex];
			break;
		}
		case 4:{
			if (refP[0] >= 2*refP[1]){ // first triangle
				triagIndex = 0; v1 = Point(1.0,0.0); v2 = Point(1.0,0.5);
			}else if ((refP[0] <= 2*refP[1]) && (refP[0] >= refP[1])){// second triangle
				triagIndex = 1; v1 = Point(1.0,0.5); v2 = Point(1.0,1.0);
			}else if ((refP[0] <= refP[1]) && (2*refP[0] >= refP[1])){// third triangle
				triagIndex = 2; v1 = Point(1.0,1.0); v2 = Point(0.5,1.0);
			}else{ // fourth triangle
				triagIndex = 3;v1 = Point(0.5,1.0); v2 = Point(0.0,1.0);
			}
			refInv = case3[triagIndex];
			break;
		}
		default:{
			TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error,"get3DRealCoordsOnSurf , too many triangles " << nrTriag);
		}
	}
	// get the three points of the triangle
	Point p0 = edgeIntersectPoint[triangles[3*triagIndex]];
	Point p1 = edgeIntersectPoint[triangles[3*triagIndex+1]];
	Point p2 = edgeIntersectPoint[triangles[3*triagIndex+2]];
	SUNDANCE_MSG3(verb, tabs << "get3DRealCoordsOnSurf p0="<<p0<<" , p1="<<p1<<" , p2="<<p2);
	SUNDANCE_MSG3(verb, tabs << "get3DRealCoordsOnSurf v1=" << v1 << " , v2=" << v2 << " , triagIndex=" << triagIndex);

	// transfor to the reference 2D coordinates
	double ref_x = refInv[0] * refP[0] + refInv[1]*refP[1];
	double ref_y = refInv[2] * refP[0] + refInv[3]*refP[1];
	SUNDANCE_MSG3(verb, tabs << " after traf to ref_x ="<<ref_x<<" , ref_y="<<ref_y );

	// transform from 2D ref to 3D triangle real coordinate
	realPoint = p0 + ref_x*(p1-p0) + ref_y*(p2-p0);
	SUNDANCE_MSG3(verb, tabs << " end get3DRealCoordsOnSurf , realPoint="<<realPoint );
	SUNDANCE_MSG3(verb, tabs << " end get3DRealCoordsOnSurf cellPoints[0]="<<cellPoints[0]<<" , cellPoints[7]="<<cellPoints[7] );

	// transform the point in real coordinates back to reference 3D coordinates
	// this works only for cubes (structured)
	realPoint[0] = (realPoint[0] - cellPoints[0][0]) / (cellPoints[7][0] - cellPoints[0][0]);
	realPoint[1] = (realPoint[1] - cellPoints[0][1]) / (cellPoints[7][1] - cellPoints[0][1]);
	realPoint[2] = (realPoint[2] - cellPoints[0][2]) / (cellPoints[7][2] - cellPoints[0][2]);
	SUNDANCE_MSG3(verb, tabs << " end get3DRealCoordsOnSurf triagIndex="<<triagIndex<<" , realPoint="<<realPoint );

}


bool CurveIntegralCalc::usePolygoneCurveQuadrature(
		const CellType& cellType ,
		const Mesh& mesh ,
		const ParametrizedCurve& paramCurve){

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
		{ break; }
		default: {}
		}
		break;
	}
	case QuadCell:
	{
		switch(cellType)
		{
		case QuadCell:
		{
			// this works only in 2D for polygons, it is important that we ask the underlying object and not the object itself
			const Polygon2D* poly = dynamic_cast<const Polygon2D* > ( paramCurve.ptr().get()->getUnderlyingCurveObj() );
			if ( poly != 0 ){
				return true;
			}
			break;
		}
		default: {}
		}
		break;
	}
	default: { }
	}
	//
	return false;
}


void CurveIntegralCalc::getCurveQuadPoints_polygon(
                               int  maxCellLID ,
                               CellType  maxCellType ,
                               const Mesh& mesh ,
                               const Array<Point>& cellPoints,
							   const ParametrizedCurve& paramCurve,
							   const QuadratureFamily& quad ,
							   Array<Point>& curvePoints ,
							   Array<Point>& curveDerivs ,
							   Array<Point>& curveNormals ) {

	// this works only in 2D for polygons, it is important that we ask the underlying object and not the object itself
	const Polygon2D* polygon = dynamic_cast<const Polygon2D* > ( paramCurve.ptr().get()->getUnderlyingCurveObj() );
	if ( (polygon == 0) || (maxCellType != QuadCell) ) {
		TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error , " getCurveQuadPoints_polygon , paramCurve must be a Polygon and the cell must be a quad" );
		return;
	}

	// the algorithm is to create separate line segments
	//(P1,P2,d1),(P2,P3,d2),(P4,P5,d3) ... d1,d2,d3 are normalized lengths
	// see these line segments as one continuous line and maps the quadrature points on the segments

	Array<Point> allIntPoints(0);
	Array<int> allIntPointsEdge(0);
	Array<Point> allPolygonPoints(0);
	int nrPoints , total_points , polyNrPoint ;
	Tabs tabs;
	int verb = 0;
	double eps = 1e-14;

	SUNDANCE_MSG3(verb, " ---------- START METHOD -----------  ");

	//
	switch (maxCellType){
	  case QuadCell:{
	  	// loop over each edge and detect intersection point
	   	// there can be only one
		TEUCHOS_TEST_FOR_EXCEPTION( cellPoints.size() != 4 ,
	   			std::runtime_error ," CurveIntegralCalc::getCurveQuadPoints , QuadCell must have 4 points" );
		//SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_line, on QuadCell")
	   	Array<Point> result(5);
	   	int edegIndex[4][2] = { {0,1} , {0,2} , {1,3} , {2,3} };
	   	bool pointExists;

	   	total_points = 0;
	   	// loop over the edges
	   	for (int jj = 0 ; jj < 4 ; jj++ ){
			paramCurve.returnIntersectPoints(cellPoints[edegIndex[jj][0]], cellPoints[edegIndex[jj][1]], nrPoints , result);
			// add the intersection points to the list of intersection points
	    	for (int ptmp = 0 ; ptmp < nrPoints ; ptmp++){
	    		pointExists = false;
	    		//SUNDANCE_MSG3(verb, tabs << "getCurveQuadPoints_line  INVESTIGATE Intersect point:" << result[ptmp]);
	    		for (int tmp_r = 0 ; tmp_r < total_points ; tmp_r++){
	    			Point tmpPoint = result[ptmp];
	    			tmpPoint = tmpPoint - allIntPoints[tmp_r];
	    			if ( ::sqrt(tmpPoint * tmpPoint) < eps ) pointExists = true;
	    		}
	    		// add only if this point does not exist in the list
	    		if (!pointExists){
	    			//SUNDANCE_MSG3(verb, tabs << "getCurveQuadPoints_line  ADD  Intersect point:" << result[ptmp]);
	    			allIntPoints.resize( total_points + 1 );
	    			allIntPointsEdge.resize( total_points + 1 );
	    			allIntPoints[total_points] = result[ptmp];
	    			allIntPointsEdge[total_points] = jj;
	    			total_points = total_points + 1;
	    		}
	    	}
	   	}
	    // test if we have too much intersection points
		//SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_line, end finding intersection points")
     } break;
     case TriangleCell:{
    	TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error , " CurveIntegralCalc::getCurveQuadPoints_polygon , TriangleCell not implemented yet" );
     } break;
     default : {
    	TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error , "CurveIntegralCalc::getCurveQuadPoints_polygon , Unknown Cell in 2D" );
	 }
	}

	// now we have all the intersection points which are not the same (duplicates are eliminated)
	// we get the points from the polygon
	polygon->getCellsPolygonesPoint( maxCellLID ,allPolygonPoints);
	polyNrPoint = allPolygonPoints.size();

	// now we should merge the intersection points and the points from the polygon
	// into a line segments

	nrPoints = polyNrPoint + total_points;
	Array<Point> segmentStart(0);
	Array<Point> segmentEnd(0);
	int segmentNr = 0 , flagMode , edgeActual = 0;
	Array<bool> usedIntPoint( total_points , false );

	// todo:
	if (total_points < 2) {
		std::cout << "total_points = " << total_points << " , P0:" << allIntPoints[0] << std::endl;
	}
	if (total_points == 3) {
		std::cout << "total_points = " << total_points << "P0:" << allIntPoints[0] <<
				" , P1:" << allIntPoints[1] << " , P2:" << allIntPoints[2] << std::endl;
		total_points = 2;
	}
	if (total_points > 4) {
		std::cout << "total_points = " << total_points << "P0:" << allIntPoints[0] << " , P1:" << allIntPoints[1]
		       << " , P2:" << allIntPoints[2] << " , P3:" << allIntPoints[3] << std::endl;
		total_points = 4;
	}

	TEUCHOS_TEST_FOR_EXCEPTION( (total_points != 2) && (total_points != 4) , std::runtime_error , "nr Intersect point can be eighter 2 or 4 but is " << total_points);

	// ---- create total_points/2 different segments ----
	flagMode = 0;
	for (int intP = 0 ; intP < total_points ; intP++ ){
		if (flagMode == 0){
			// get a point which is not used
			int Pind = -1;
			for (int stmp = 0 ; stmp < total_points ; stmp++){
				// this point was not used then take it
				if (!usedIntPoint[stmp]) { Pind = stmp; break;}
			}
			//SUNDANCE_MSG3(verb, tabs << "getCurveQuadPoints_line flag mode 0 , found Index Pind:" << Pind << " , start:"<<allIntPoints[Pind]);
			segmentStart.resize( segmentNr + 1);
			segmentEnd.resize( segmentNr + 1 );
			edgeActual = allIntPointsEdge[Pind];
			segmentStart[ segmentNr ] = allIntPoints[Pind];
			usedIntPoint[ Pind ] = true;
			flagMode = 1;
		} else {
			// here we have a start point we are looking for an end of the segment line
			double dist = 1e+100;
			int Pind = -1;
			for (int stmp = 0 ; stmp < total_points ; stmp++){
				// the condition when to test one intersection point
				if ( !usedIntPoint[stmp] && ( (edgeActual != allIntPointsEdge[stmp]) || total_points < 3 ) ){
					Point tmpPointC = (segmentStart[segmentNr] - allIntPoints[stmp]);
					double thisDist = ::sqrt(tmpPointC*tmpPointC);
					if ( thisDist < dist ) {
						dist = thisDist;
						Pind = stmp;
					}
				}
			}
			//SUNDANCE_MSG3(verb, tabs << "getCurveQuadPoints_line flag mode 1 , found Index Pind:" << Pind << " , end:"<<allIntPoints[Pind]);
			segmentEnd[ segmentNr ] = allIntPoints[Pind];
			usedIntPoint[ Pind ] = true;
			flagMode = 0;
			segmentNr = segmentNr + 1;
		}
	}

	// ---- at this point we have the start and end points of the segments
	// ---- now we should insert the polygon points
	for (int polyPoint = 0 ; polyPoint < polyNrPoint ; polyPoint++ ){
		// find the line segment which will be split
		double totalDist = 1e+100;
		int segInd = -1;
		for (int pl = 0; pl < segmentNr ; pl++ ){
			Point d1p = segmentStart[pl] - allPolygonPoints[polyPoint];
			double d1 = ::sqrt(d1p*d1p);
			Point d2p = segmentEnd[pl] - allPolygonPoints[polyPoint];
			double d2 = ::sqrt(d2p*d2p);
			Point dist_prev = segmentEnd[pl] - segmentStart[pl];
			double d_p = ::sqrt(dist_prev*dist_prev);
			//SUNDANCE_MSG3(verb, tabs << " Probe actual segment START :" << segmentStart[pl] << " , END:" << segmentEnd[pl] );
			//SUNDANCE_MSG3(verb, tabs << " Probe actual point to pl:" << pl << " actual dist:" << (d1+d2-d_p) << " , totalDist:" << totalDist);
			// test if this is the shortest segment till now where to insert the polygon point
			if ( (d1 + d2 - d_p < totalDist) && (d1 > eps) && (d2 > eps) ) {
				totalDist = d1 + d2 - d_p;
				segInd = pl;
			}
		}

		//SUNDANCE_MSG3(verb, tabs << " Add polygon point to segInd:" << segInd << " , polyPoint:" << polyPoint
		//		<< " , allPolygonPoints[polyPoint]:" << allPolygonPoints[polyPoint]);
		// if the point is already in or no segment found
		if ( segInd < 0) continue;

		// "segInd" shows the segment which should be split
		segmentStart.resize( segmentNr + 1);
		segmentEnd.resize( segmentNr + 1 );
		// enlarge first the vector
		for (int tmpP = segmentNr-1; tmpP >= segInd ; tmpP--){
			segmentStart[tmpP+1] = segmentStart[tmpP];
			segmentEnd[tmpP+1] = segmentEnd[tmpP];
		}

		// insert the new line segment
		segmentEnd[segInd] = allPolygonPoints[polyPoint];
		segmentStart[segInd+1] = allPolygonPoints[polyPoint];
		segmentNr = segmentNr + 1;
	}

	// --- now we have all line segments now calculate each length
	Array<double> segmentLengthRel(0);
	Array<double> segmentStartLengthRel(0);
	segmentLengthRel.resize(segmentNr);
	segmentStartLengthRel.resize(segmentNr);

	// calculate the length of each line, so that we can find
	double totalLength = 0.0 , actLen = 0.0;
	for (int seg = 0 ; seg < segmentNr ; seg++){
		totalLength = totalLength + ::sqrt( (segmentEnd[seg]-segmentStart[seg])*(segmentEnd[seg]-segmentStart[seg]) );
		//SUNDANCE_MSG3(verb, tabs << " total length seg:" << seg << " , totalLength = " << totalLength);
		//SUNDANCE_MSG3(verb, tabs << "  segmentStart[seg] = " << segmentStart[seg] <<" segmentEnd[seg]:" << segmentEnd[seg] );
	}

    // update the lengths of each segment
	for (int seg = 0 ; seg < segmentNr ; seg++){
		double thisLength = ::sqrt( (segmentEnd[seg]-segmentStart[seg])*(segmentEnd[seg]-segmentStart[seg]) );
		segmentLengthRel[seg] = thisLength/totalLength;
		segmentStartLengthRel[seg] = actLen/totalLength;
		//SUNDANCE_MSG3(verb, tabs << " Absolute length seg:" << seg << " , segmentLengthRel[seg] = " << segmentLengthRel[seg]
		//                                << " , segmentStartLengthRel[seg]=" << segmentStartLengthRel[seg]);
		//SUNDANCE_MSG3( verb , tabs << " segmentStart[seg]: " << segmentStart[seg] << " , segmentEnd[seg]:" << segmentEnd[seg]);
		actLen = actLen + thisLength;
	}

	//try to convert the quadrature to a polygon curve integral class
	// this quadrature class is special for polygon curve integrals
	const PolygonQuadrature* polyquad = dynamic_cast<const PolygonQuadrature*> (quad.ptr().get());
	if ( polyquad == 0 ) {
		TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error ,
			" getCurveQuadPoints_polygon , Quadrature for a polygon curve must be PolygonQuadrature" );
		return;
	}

	// get the quadrature points for line segments
	Array<Point> quadPoints;
	Array<double> quadWeights;
	polyquad->getPoints( LineCell , quadPoints , quadWeights );
	int nr_quad_points = quadPoints.size();
	Point endPoint , startPoint;
	int maxLineSegments = PolygonQuadrature::getNrMaxLinePerCell();
	Array<int> cellLID(1); cellLID[0] = maxCellLID;
	CellJacobianBatch JBatch;
	JBatch.resize(1, 2, 2);
	mesh.getJacobians(2, cellLID , JBatch);

	TEUCHOS_TEST_FOR_EXCEPTION( maxLineSegments < segmentNr , std::runtime_error , " getCurveQuadPoints_polygon , nr of polygon lines inside one cell too high" <<
			        " Use PolygonQuadrature::getNrMaxLinePerCell(" << segmentNr << ") , now it is only: " << maxLineSegments );

	// ---- at this stage we have the line segments,
	// we just need to map the quadrature points to the line and calculate the normal vector
	curvePoints.resize(nr_quad_points);
	curveDerivs.resize(nr_quad_points);
	curveNormals.resize(nr_quad_points);

	// get the line segment which contains this quadrature point
	int ii = 0;
	for (int seg = 0; seg < segmentNr ; seg++ )
	{
		// the start and end points of the line segment
		startPoint = segmentStart[seg];
		endPoint = segmentEnd[seg];

		// the calculation of the normal vector is usual, we look to the middle of the segment and see if it is inside or outside
		// since this is a polygon this must be determined
		Point curvderi = endPoint - startPoint;
		Point norm_vec = Point( -curvderi[1] , curvderi[0] );
		norm_vec = norm_vec / ::sqrt(norm_vec*norm_vec);
		double relative_length = 1e-2 * ::sqrt(curvderi*curvderi);
		bool invNormVect = false;
		// this segment must be a line, and we test which part is inside
		if ( paramCurve.curveEquation( (0.5*endPoint + 0.5*startPoint) + relative_length*norm_vec ) > 0 ) {
			invNormVect = true;
		}
		else {
			invNormVect = false;
		}

		// loop over each quadrature point for one line segment
		for (int pointNr = 0 ; pointNr < nr_quad_points/maxLineSegments ; pointNr++ , ii++ )
		{
			// now we have
			//SUNDANCE_MSG3(verb, tabs << "pointNr:" << pointNr << " , nr:" << ii << " Line quad point:" << quadPoints[ii] << ", line:" << (endPoint - startPoint)
			//	<< " , segmentsLength: " << segmentStartLengthRel[seg] );
			//SUNDANCE_MSG3(verb, tabs << " nr:" << ii << " startPoint:" <<  startPoint << " , endPoint:" <<  endPoint );
			// here we map the quadrature point on the line segment
			curvePoints[ii] = startPoint + quadPoints[ii][0]*(endPoint - startPoint);

			// here we transform the real coordinate back to the
			double pointc[2] = { curvePoints[ii][0] - cellPoints[0][0] , curvePoints[ii][1] - cellPoints[0][1] };
			JBatch.applyInvJ(0, pointc , 1 , false);
			curvePoints[ii][0] = pointc[0];
			curvePoints[ii][1] = pointc[1];
			//SUNDANCE_MSG3(verb, tabs << " quad point nr:" << ii << " = " << curvePoints[ii]);

			// the derivative we calculate by simply take the total length
			// in the curve integral we use only the norm of this point
			curveDerivs[ii] = Point( segmentLengthRel[seg]*totalLength*((double)maxLineSegments) , 0.0 );
			//SUNDANCE_MSG3(verb, tabs << " curve Derivatives point nr:" << ii << " = " << curveDerivs[ii]);

			// set the normal vector
			if ( invNormVect ) {
				curveNormals[ii] = -norm_vec;
			}
			else {
				curveNormals[ii] = norm_vec;
			}
			//SUNDANCE_MSG3(verb, tabs << " curve Normals point nr:" << ii << " = " << curveNormals[ii]);
		}
	}

	// for the rest of the points just set zero values , these point should not have influence
	for  ( ; ii < nr_quad_points ; ii++ ){
		curvePoints[ii] = Point( 0.0 , 0.0 );
		curveDerivs[ii] = Point( 0.0 , 0.0 );
		curveNormals[ii] = Point( 0.0 , 0.0 );
	}

	SUNDANCE_MSG3(verb, " ---------- END METHOD -----------  ");
}


void CurveIntegralCalc::getSurfQuadPoints(
                               int  maxCellLID ,
                               CellType  maxCellType ,
                               const Mesh& mesh ,
                               const Array<Point>& cellPoints,
							   const ParametrizedCurve& paramCurve,
							   const QuadratureFamily& quad ,
							   Array<Point>& curvePoints ,
							   Array<Point>& curveDerivs ,
							   Array<Point>& curveNormals ) {

	const SurfQuadrature* surfquad = dynamic_cast<const SurfQuadrature*> (quad.ptr().get());
	if ( surfquad == 0 ) {
		TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error ,
			" getSurfQuadPoints , Surface Quadrature for a 3D must be SurfQuadrature" );
		return;
	}

	Array<Point> quadPoints;
	Array<double> quadWeights;
	surfquad->getPoints( TriangleCell , quadPoints , quadWeights );
	int nr_quad_points = quadPoints.size();
	int maxTriangles = SurfQuadrature::getNrMaxTrianglePerCell();
	int nrQuadPointPerTriag = nr_quad_points/maxTriangles;
	int verb = 0;
	curvePoints.resize(nr_quad_points);
	curveDerivs.resize(nr_quad_points);
	curveNormals.resize(nr_quad_points);

	Array<Point> intersectPoints;
	Array<int> edgeIndex;
	Array<int> triangleIndex;
    // the results in this form will be in physical coordinates
	SundanceSurf3DCalc::getCurveQuadPoints( maxCellType , maxCellLID , mesh , paramCurve, cellPoints,
						 intersectPoints ,  edgeIndex , triangleIndex );

	int  nr_triag = triangleIndex.size()/3  ,  qind = 0;
	double totalArea = 0.0 , area;
	// first compute the total area of all triangles
	for (int t = 0 ; t < nr_triag ; t++){
		// for each triangle
		Point p0 = intersectPoints[triangleIndex[3*t]];
		Point p1 = intersectPoints[triangleIndex[3*t + 1]];
		Point p2 = intersectPoints[triangleIndex[3*t + 2]];
		Point n = cross(p1-p0,p2-p0);
		totalArea = totalArea + ::sqrt(n*n);
	}


	SUNDANCE_MSG3(verb, " nr_triag = " << nr_triag);
	// for each resulted triangles compute the quadrature points and the normal
	for (int t = 0 ; t < nr_triag ; t++){
		// for each triangle
		Point p0 = intersectPoints[triangleIndex[3*t]];
		Point p1 = intersectPoints[triangleIndex[3*t + 1]];
		Point p2 = intersectPoints[triangleIndex[3*t + 2]];

		Point n = cross(p1-p0,p2-p0);
		// calculate area
		area = ::sqrt(n*n); // / 2.0; -> no division by two
		// invert the normal vector if necessary
		if ( paramCurve.curveEquation( p0 + 5e-3*n ) > 0 ) { n = (-1)*n; }
		n = (1/::sqrt(n*n))*n; // normalize the normal vector

		// for each quadrature point on this triangle
		//SUNDANCE_MSG3(verb, " p0 = " << p0 << " , p1 = " << p1 << " , p2 = " << p2 );
		//SUNDANCE_MSG3(verb, " area = " << area << " , n = " << n );
		for (int q = 0 ; q < nrQuadPointPerTriag ; q++ , qind++){
			// this is the point in real coordinates
			Point pTmp = p0 + quadPoints[qind][0] * (p1-p0) + quadPoints[qind][1] * (p2-p0);
			//SUNDANCE_MSG3(verb, "REAL curvePoints["<<qind<<"]=" << pTmp );
			// here we transform it to reference coordinates
			curvePoints[qind] = Point( (pTmp[0] -  cellPoints[0][0])/(cellPoints[7][0] - cellPoints[0][0]),
					                   (pTmp[1] -  cellPoints[0][1])/(cellPoints[7][1] - cellPoints[0][1]),
					                   (pTmp[2] -  cellPoints[0][2])/(cellPoints[7][2] - cellPoints[0][2]) );
			//SUNDANCE_MSG3(verb, "REF curvePoints["<<qind<<"]=" << curvePoints[qind] );
			curveDerivs[qind] = Point( area * ((double)maxTriangles) , 0.0 , 0.0 );
			curveNormals[qind] = n;
		}
	}

	// the rest of the quadrature points we set to zero
	for ( ; qind < nr_quad_points ; qind++){
		curvePoints[qind] = Point( 0 , 0 , 0 );
		curveDerivs[qind] = Point( 0 , 0 , 0 );
		curveNormals[qind] = Point( 0 , 0 , 0 );
	}
	SUNDANCE_MSG3(verb, " ---------- END METHOD -----------  ");
}
