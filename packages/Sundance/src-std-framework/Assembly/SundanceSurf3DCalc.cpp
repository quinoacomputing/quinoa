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
 * SundanceSurf3DCalc.cpp
 *
 *  Created on: Oct 25, 2011
 *      Author: benk
 */

#include "SundanceSurf3DCalc.hpp"

using namespace Sundance;

const int SundanceSurf3DCalc::edegIndex[12][2] = { {0,1} , {0,2} , {0,4} , {1,3} , {1,5} , {2,3} , {2,6} , {3,7} ,
		    			                 {4,5} , {4,6} , {5,7} , {6,7} };
                                              //       0           1           2           3            4             5
const int SundanceSurf3DCalc::faceEdges[6][4]  = { {0,1,3,5} , {0,2,4,8} , {1,2,6,9} , {3,4,7,10} , {5,6,7,11} , {8,9,10,11} };

                                              //     0       1       2       3       4       5       6       7      8       9       10      11
const int SundanceSurf3DCalc::edgeFaces[12][2] = { {0,1} , {0,2} , {1,2} , {0,3} , {1,3} , {0,4} , {2,4} , {3,4} , {1,5} , {2,5} , {3,5} , {4,5} };

//                                                            0   1   2   3   4   5   6   7   8   9   10  11
const int SundanceSurf3DCalc::edgeNeighboredges[12][12] = { { 1 , 1 , 1 , 1 , 1 , 1 , 0 , 0 , 1 , 0 , 0 , 0} ,   /* 0 */
		                                                    { 1 , 1 , 1 , 1 , 0 , 1 , 1 , 0 , 0 , 1 , 0 , 0} ,   /* 1 */
		                                                    { 1 , 1 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 1 , 0 , 0} ,   /* 2 */
		                                                    { 1 , 1 , 0 , 1 , 1 , 1 , 0 , 1 , 0 , 0 , 1 , 0} ,   /* 3 */
		                                                    { 1 , 0 , 1 , 1 , 1 , 0 , 0 , 1 , 1 , 0 , 1 , 0} ,   /* 4 */
		                                                    { 1 , 1 , 0 , 1 , 0 , 1 , 1 , 1 , 0 , 0 , 0 , 1} ,   /* 5 */
		                                                    { 0 , 1 , 1 , 0 , 0 , 1 , 1 , 1 , 0 , 1 , 0 , 1} ,   /* 6 */
		                                                    { 0 , 0 , 0 , 1 , 1 , 1 , 1 , 1 , 0 , 0 , 1 , 1} ,   /* 7 */
		                                                    { 1 , 0 , 1 , 0 , 1 , 0 , 0 , 0 , 1 , 1 , 1 , 1} ,   /* 8 */
		                                                    { 0 , 1 , 1 , 0 , 0 , 0 , 1 , 0 , 1 , 1 , 1 , 1} ,   /* 9 */
		                                                    { 0 , 0 , 0 , 1 , 1 , 0 , 0 , 1 , 1 , 1 , 1 , 1} ,   /* 10 */
		                                                    { 0 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1 , 1 , 1 , 1} , };/* 11 */

void SundanceSurf3DCalc::getCurveQuadPoints(CellType  maxCellType ,
    		                       int maxCellLID ,
    		                       const Mesh& mesh ,
								   const ParametrizedCurve& paramCurve,
								   const Array<Point>& brickPoints,
								   Array<Point>& intersectPoints ,
								   Array<int>& edgeIndex,
								   Array<int>& triangleIndex){
	TEUCHOS_TEST_FOR_EXCEPTION( maxCellType != BrickCell ,
					std::runtime_error , " maxCellType must be a BrickCell: " << BrickCell << " but it is:" << maxCellType);

	int nrPoints , edge , t_inters_points = 0;
	int verb = 2;
	int tmp_i;
	int faceIntPoint[6] = { 0 , 0 , 0 , 0 , 0 , 0 };
	int edgePointI[12] = { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 };

// -------------------- get the intersection point according to the policy -------------------

	// get the intersection points between the edges of the brick and the surface
	for (edge = 0 ; edge < 12 ; edge++){
    	Array<Point> result(0);
		paramCurve.returnIntersectPoints(brickPoints[edegIndex[edge][0]], brickPoints[edegIndex[edge][1]], nrPoints , result);
		// test if we have intersection point, and take this only if there is one
    	if (nrPoints == 1){
    		// we test if the impacted faces have already two intersection points, if yes then, ignore the intersection point
    		// ignore this intersection point
    		if ( (faceIntPoint[ edgeFaces[edge][0]] < 2) &&  (faceIntPoint[ edgeFaces[edge][1] ] < 2 )) {
    		   faceIntPoint[ edgeFaces[edge][0] ]++; faceIntPoint[ edgeFaces[edge][1] ]++;
    		   intersectPoints.resize( t_inters_points + 1 );
    		   edgeIndex.resize( t_inters_points + 1 );
    		   intersectPoints[t_inters_points] = result[0];
    		   edgeIndex[t_inters_points] = edge;
    		   edgePointI[edge] = t_inters_points;
    		   //SUNDANCE_MSG3(verb, "found Int. point [" << t_inters_points <<"]=" << result[0]);
    		   t_inters_points++;
    		} else {
    			SUNDANCE_MSG1(verb, "WARNING One face had already two intersection points :"
    					<< faceIntPoint[ edgeFaces[edge][0]] <<" , " << faceIntPoint[ edgeFaces[edge][1]]);
    		}
    	}
    	// if there is more than one point than we ignore all of them, and just print a warning
    	if (nrPoints > 1){
    		SUNDANCE_MSG1(verb, " WARNING: nr of intersection points: " << nrPoints << " for edge [ " << brickPoints[edegIndex[edge][0]] << " , "
    				<< brickPoints[edegIndex[edge][1]] << "] , the first two points i1:" << result[0] << " , i2: " << result[1]);
    	}
	}

	// this means no surface point has found
	if ( t_inters_points < 3){
		intersectPoints.resize(0);
		edgeIndex.resize(0);
		triangleIndex.resize(0);
		SUNDANCE_MSG1(verb, " WARNING: less than 3 int point found : " << t_inters_points << " for maxCellLID=" << maxCellLID
				        << " P[0]=" << brickPoints[0] );
		return;
	}

// ------------------------ discover the triangles which form the surf -------------------------
	// get the polygon
	Array<int> polygonList(t_inters_points,-1);
	Array<int> pointUsed(t_inters_points,-1);
	int startIndex = -1 ;
	polygonList[0] = 0; pointUsed[0] = 1;
	SUNDANCE_MSG3(verb, " polygonList[ 0 ] = 0 , P: " << intersectPoints[ 0 ]);
	for (int jj = 1 ; jj < t_inters_points ; jj++){
		// look for a point which is connected to polygonList[jj] ad has not been used yet
		tmp_i = -1;
		for (int ii = 0 ; ii < t_inters_points ; ii++){
			// test if this edge is neighboring
			if ( (pointUsed[ii] < 0) && (edgeNeighboredges[ edgeIndex[polygonList[jj-1]] ][ edgeIndex[ii] ] > 0) ) {
				// we found one neighboring point
				tmp_i = ii;
				polygonList[jj] = ii;
				pointUsed[ii] = 1;
				SUNDANCE_MSG3(verb, " polygonList[ " << jj << " ] = " << ii << " , P: " << intersectPoints[ ii ]);
				break;
			}
		}
		TEUCHOS_TEST_FOR_EXCEPTION( tmp_i < 0 , std::runtime_error , " error by polygon detecting tmp_i > -1 , might be that the surface is not KONVEX");
		// look for an intersection point which is in the same face as the current intersection point
	}

	//get the starting point by getting the point which is closest
	Array<double> dist(t_inters_points);
	double dist_tmp = 1e+100;
	for (int jj = 0 ; jj < t_inters_points ; jj++){
		double sum = 0.0;
		for (int ii = 0 ; ii < t_inters_points ; ii++){
			sum = sum + ::sqrt( (intersectPoints[ polygonList[ii] ] - intersectPoints[ polygonList[jj] ])
					*(intersectPoints[ polygonList[ii] ] - intersectPoints[ polygonList[jj] ]) );
		}
		// store the minimum distance and the index
		if (sum < dist_tmp) { startIndex = jj; dist_tmp = sum; }
	}
	SUNDANCE_MSG3(verb, " central point found: " << startIndex << " P = " << intersectPoints[ polygonList[startIndex] ]);

// ----------------- having the polygon and the starting point define the triangles ---------------
	int fI = 0 , actIF = startIndex , actIL = startIndex;
	triangleIndex.resize(3*(t_inters_points-2),-1);
	for (int jj = 0 ; jj < t_inters_points-2 ; jj++){
		actIF = (actIF+1) % t_inters_points;
		actIL = (actIF+1) % t_inters_points;
		SUNDANCE_MSG3(verb, " startIndex = " << startIndex << " , actIF = " << actIF << " , actIL = " << actIL);
		triangleIndex[fI]   = polygonList[startIndex];
		triangleIndex[fI+1] = polygonList[actIF];
		triangleIndex[fI+2] = polygonList[actIL];
		SUNDANCE_MSG3(verb, "add triangle [ " << polygonList[startIndex] << " , " << polygonList[actIF] << " , " << polygonList[actIL] << " ]");
		SUNDANCE_MSG3(verb, "add triangle points = [ " << intersectPoints[polygonList[startIndex]] << " , " << intersectPoints[polygonList[actIF]]
		                      << " , " << intersectPoints[ polygonList[actIL] ] << " ]");
		fI = fI + 3;
	}
// end the triangle is set
}
