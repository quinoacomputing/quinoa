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
 * SundanceTriangleSurf3D.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: benk
 */

#include "SundanceTriangleSurf3D.hpp"
#include "SundancePoint.hpp"
#include "SundanceMesh.hpp"
#include "SundanceDefs.hpp"
#include "PlayaMPIComm.hpp"

#include <iostream>
#include <fstream>

using namespace Sundance;
using Playa::MPIOp;
using Playa::MPIDataType;

TriangleSurf3D::TriangleSurf3D(const Mesh& mesh , const Array<Point>& points ,
		const Array<int>& triags , double a1, double a2, bool flipD ) :
	CurveBase(2, a1, a2, flipD), hasMesh_(true), mesh_(&(mesh)) ,
	triagSurfSpaceValues_() , triagSurfSpaceNames_() , nrScalarField_(0) , isRecording_(false), twinTriangleSurf_(0) // here we twist the alpha parameters for the polygon convention
{
	int verb = 0;
	// just store the input points
	triagPoints_.resize(points.size());
	SUNDANCE_MSG3( verb , "TriangleSurf3D() Ctor nrPoints=" << points.size() );
	maxX_ = minX_ = points[0][0];
	maxY_ = minY_ = points[0][1];
	maxZ_ = minZ_ = points[0][2];
	for (int j = 0 ; j < points.size() ; j++ ){
		triagPoints_[j] = points[j];
		SUNDANCE_MSG3( verb , " point[" << j << "] = " << triagPoints_[j]);
		maxX_ = (maxX_ > points[j][0]) ? maxX_ : points[j][0];
		minX_ = (minX_ < points[j][0]) ? minX_ : points[j][0];
		maxY_ = (maxY_ > points[j][1]) ? maxY_ : points[j][1];
		minY_ = (minY_ < points[j][1]) ? minY_ : points[j][1];
		maxZ_ = (maxZ_ > points[j][2]) ? maxZ_ : points[j][2];
		minZ_ = (minZ_ < points[j][2]) ? minZ_ : points[j][2];
	}
	pointMaxCellLID_.resize(triagPoints_.size(),-1);

	// store the triangles points
	triagIndexes_.resize(triags.size());
	for (int j = 0 ; j < triags.size() ; j++ ){
		triagIndexes_[j] = triags[j];
	}
	triagIDs_.resize(triagIndexes_.size()/3);
	for (int j = 0 ; j < triagIDs_.size() ; j++ ){ triagIDs_[j] = j; }
	computeMaxCellLIDs();
	// set up space tree
	setupSpaceTree();
}

TriangleSurf3D::TriangleSurf3D(const Array<Point>& points , const Array<int>& triags ,
		double a1, double a2, bool flipD ) :
	CurveBase(2, a1, a2, flipD), hasMesh_(false), mesh_(0),
	triagSurfSpaceValues_() , triagSurfSpaceNames_() , nrScalarField_(0) , isRecording_(false) ,twinTriangleSurf_(0) // here we twist the alpha parameters for the polygon convention
{
	int verb = 0;
	// just store the input points
	triagPoints_.resize(points.size());
	SUNDANCE_MSG3( verb , "TriangleSurf3D() Ctor nrPoints=" << points.size() );
	maxX_ = minX_ = points[0][0];
	maxY_ = minY_ = points[0][1];
	maxZ_ = minZ_ = points[0][2];
	for (int j = 0 ; j < points.size() ; j++ ){
		triagPoints_[j] = points[j];
		SUNDANCE_MSG3( verb , " point[" << j << "] = " << triagPoints_[j]);
		maxX_ = (maxX_ > points[j][0]) ? maxX_ : points[j][0];
		minX_ = (minX_ < points[j][0]) ? minX_ : points[j][0];
		maxY_ = (maxY_ > points[j][1]) ? maxY_ : points[j][1];
		minY_ = (minY_ < points[j][1]) ? minY_ : points[j][1];
		maxZ_ = (maxZ_ > points[j][2]) ? maxZ_ : points[j][2];
		minZ_ = (minZ_ < points[j][2]) ? minZ_ : points[j][2];
	}
	pointMaxCellLID_.resize(triagPoints_.size(),-1);

	// store the triangles points
	triagIndexes_.resize(triags.size());
	for (int j = 0 ; j < triags.size() ; j++ ){
		triagIndexes_[j] = triags[j];
	}
	triagIDs_.resize(triagIndexes_.size()/3);
	for (int j = 0 ; j < triagIDs_.size() ; j++ ){ triagIDs_[j] = j; }
	// set up space tree
	setupSpaceTree();
}

TriangleSurf3D::TriangleSurf3D(const Mesh& mesh , const std::string& filename , double a1, double a2, bool flipD ) :
CurveBase(2, a1, a2, flipD), hasMesh_(true), mesh_(&(mesh)),
triagSurfSpaceValues_() , triagSurfSpaceNames_() , nrScalarField_(0) , isRecording_(false), twinTriangleSurf_(0) // here we twist the alpha parameters for the polygon convention
{
	std::string str_tmp;
	std::ifstream myfile;
	std::vector<std::string> elems;
	int verb = 6;

	char *pEnd , *pEnd1;
	double d1 , d2, d3;
	int    i1 , i2, i3;
	int nr_points = 0 , nr_triags;

	SUNDANCE_MSG3( verb , "TriangleSurf3D() read from file ");

	myfile.open( filename.c_str() , std::ios::in );

	// the first line should contain the number of points
	getline( myfile , str_tmp );
	nr_points = round(strtod (str_tmp.c_str() , &pEnd));

	// each line of the file should contain a point
	triagPoints_.resize(nr_points);
    for (int p = 0 ; p < nr_points ; p++)
    {
    	// read the line
    	getline( myfile , str_tmp );
    	d1 = strtod (str_tmp.c_str() , &pEnd);
    	d2 = strtod (pEnd , &pEnd1);
    	d3 = strtod (pEnd1 , NULL);
		// we add the point to the polygon
		Point tmpP(d1,d2,d3);
		triagPoints_[p] = tmpP;
		SUNDANCE_MSG3( verb , " point[" << p << "] = " << triagPoints_[p]);
		if (p == 0){
			maxX_ = minX_ = tmpP[0];
			maxY_ = minY_ = tmpP[1];
			maxZ_ = minZ_ = tmpP[2];
		}else{
			maxX_ = (maxX_ > tmpP[0]) ? maxX_ : tmpP[0];
			minX_ = (minX_ < tmpP[0]) ? minX_ : tmpP[0];
			maxY_ = (maxY_ > tmpP[1]) ? maxY_ : tmpP[1];
			minY_ = (minY_ < tmpP[1]) ? minY_ : tmpP[1];
			maxZ_ = (maxZ_ > tmpP[2]) ? maxZ_ : tmpP[2];
			minZ_ = (minZ_ < tmpP[2]) ? minZ_ : tmpP[2];
		}
    }
	pointMaxCellLID_.resize(triagPoints_.size(),-1);

    // this line should contain the nu
	getline( myfile , str_tmp );
    nr_triags = round(strtod (str_tmp.c_str() , &pEnd));
    triagIndexes_.resize(3*nr_triags);
    // read in all the triangles
    for (int p = 0 ; p < 3*nr_triags ; p = p + 3)
    {
    	getline( myfile , str_tmp );
    	i1 = round(strtod (str_tmp.c_str() , &pEnd));
    	i2 = round(strtod (pEnd , &pEnd1));
    	i3 = round(strtod (pEnd1 , NULL));
    	triagIndexes_[p] = i1; triagIndexes_[p+1] = i2; triagIndexes_[p+2] = i3;
    }
	triagIDs_.resize(triagIndexes_.size()/3);
	for (int j = 0 ; j < triagIDs_.size() ; j++ ){ triagIDs_[j] = j; }
	// close the file
	myfile.close();
	computeMaxCellLIDs();
	// set up space tree
	setupSpaceTree();
}

TriangleSurf3D::TriangleSurf3D( const std::string& filename , double a1, double a2 ,bool flipD ) :
CurveBase(2, a1, a2, flipD), hasMesh_(false), mesh_(0),
triagSurfSpaceValues_() , triagSurfSpaceNames_() , nrScalarField_(0) , isRecording_(false), twinTriangleSurf_(0)// here we twist the alpha parameters for the polygon convention
{
	std::string str_tmp;
	std::ifstream myfile;
	std::vector<std::string> elems;
	int verb = 6;

	char *pEnd , *pEnd1;
	double d1 , d2, d3;
	int    i1 , i2, i3;
	int nr_points = 0 , nr_triags;

	SUNDANCE_MSG3( verb , "TriangleSurf3D() read from file ");

	myfile.open( filename.c_str() , std::ios::in );

	// the first line should contain the number of points
	getline( myfile , str_tmp );
	nr_points = round(strtod (str_tmp.c_str() , &pEnd));

	// each line of the file should contain a point
	triagPoints_.resize(nr_points);
    for (int p = 0 ; p < nr_points ; p++)
    {
    	// read the line
    	getline( myfile , str_tmp );
    	d1 = strtod (str_tmp.c_str() , &pEnd);
    	d2 = strtod (pEnd , &pEnd1);
    	d3 = strtod (pEnd1 , NULL);
		// we add the point to the polygon
		Point tmpP(d1,d2,d3);
		triagPoints_[p] = tmpP;
		SUNDANCE_MSG3( verb , " point[" << p << "] = " << triagPoints_[p]);
		if (p == 0){
			maxX_ = minX_ = tmpP[0];
			maxY_ = minY_ = tmpP[1];
			maxZ_ = minZ_ = tmpP[2];
		}else{
			maxX_ = (maxX_ > tmpP[0]) ? maxX_ : tmpP[0];
			minX_ = (minX_ < tmpP[0]) ? minX_ : tmpP[0];
			maxY_ = (maxY_ > tmpP[1]) ? maxY_ : tmpP[1];
			minY_ = (minY_ < tmpP[1]) ? minY_ : tmpP[1];
			maxZ_ = (maxZ_ > tmpP[2]) ? maxZ_ : tmpP[2];
			minZ_ = (minZ_ < tmpP[2]) ? minZ_ : tmpP[2];
		}
    }
	pointMaxCellLID_.resize(triagPoints_.size(),-1);

    // this line should contain the nu
	getline( myfile , str_tmp );
    nr_triags = round(strtod (str_tmp.c_str() , &pEnd));
    triagIndexes_.resize(3*nr_triags);
    // read in all the triangles
    for (int p = 0 ; p < 3*nr_triags ; p = p + 3)
    {
    	getline( myfile , str_tmp );
    	i1 = round(strtod (str_tmp.c_str() , &pEnd));
    	i2 = round(strtod (pEnd , &pEnd1));
    	i3 = round(strtod (pEnd1 , NULL));
    	triagIndexes_[p] = i1; triagIndexes_[p+1] = i2; triagIndexes_[p+2] = i3;
    }
	triagIDs_.resize(triagIndexes_.size()/3);
	for (int j = 0 ; j < triagIDs_.size() ; j++ ){ triagIDs_[j] = j; }
	// close the file
	myfile.close();

	// set up space tree
	setupSpaceTree();
}

void TriangleSurf3D::setMesh(const Mesh& mesh){
	mesh_ = &(mesh);
	hasMesh_ = true;
	computeMaxCellLIDs();
}

// update the state of the curve of the control point were changed
void TriangleSurf3D::update() {
	maxX_ = minX_ = triagPoints_[0][0];
	maxY_ = minY_ = triagPoints_[0][1];
	maxZ_ = minZ_ = triagPoints_[0][2];
	// get the maximum and minimum values
	for (int j = 0 ; j < triagPoints_.size() ; j++ ){
		maxX_ = (maxX_ > triagPoints_[j][0]) ? maxX_ : triagPoints_[j][0];
		minX_ = (minX_ < triagPoints_[j][0]) ? minX_ : triagPoints_[j][0];
		maxY_ = (maxY_ > triagPoints_[j][1]) ? maxY_ : triagPoints_[j][1];
		minY_ = (minY_ < triagPoints_[j][1]) ? minY_ : triagPoints_[j][1];
		maxZ_ = (maxZ_ > triagPoints_[j][2]) ? maxZ_ : triagPoints_[j][2];
		minZ_ = (minZ_ < triagPoints_[j][2]) ? minZ_ : triagPoints_[j][2];
	}
	computeMaxCellLIDs();
	setupSpaceTree();
}

void TriangleSurf3D::computeMaxCellLIDs(){

	// test if we have a mesh, otherwise
	if (hasMesh_ == false) return;

	int meshDim = mesh_->spatialDim();
	int nrLID = mesh_->numCells(meshDim);
	double eps = 1e-10;
	int verb = 6;

	SUNDANCE_MSG3( verb , " TriangleSurf3D::computeMaxCellLIDs " << triagPoints_.size() );

	// initially set everything to -1, so in parallel case we'll know which do not belong to this processor
	pointMaxCellLID_.resize(triagPoints_.size());
	for (int j = 0 ; j < triagPoints_.size() ; j++){ pointMaxCellLID_[j] = -1; }

	for (int cellLID = 0 ; cellLID < nrLID ; cellLID++){
				int tmp1;
				Point p0 = mesh_->nodePosition( mesh_->facetLID(meshDim,cellLID,0,0,tmp1) );
				Point p7 = mesh_->nodePosition( mesh_->facetLID(meshDim,cellLID,0,7,tmp1) );
				//SUNDANCE_MSG3( verb , " Polygon2D::computeMaxCellLIDs cellLID =" << cellLID << " ,p0:" << p0 << " ,p7:" << p7 );
				for (int j = 0 ; j < triagPoints_.size() ; j++)
				{
	               //SUNDANCE_MSG3( verb , " test j = " << j );
                   Point& tmp = triagPoints_[j];
                   if ( (tmp[0] >= p0[0]-eps) &&  (tmp[0] <= p7[0]+eps) &&
                		(tmp[1] >= p0[1]-eps) &&  (tmp[1] <= p7[1]+eps) &&
                        (tmp[2] >= p0[2]-eps) &&  (tmp[2] <= p7[2]+eps)
                        && (pointMaxCellLID_[j] < 0))
                   {
                	   // store the maxCellLID
                	   //SUNDANCE_MSG3( verb , " TriangleSurf3D::computeMaxCellLIDs j=" << j << " maxCellLID=" << cellLID );
                	   pointMaxCellLID_[j] = cellLID;
                   }
				}
	}
	for (int j = 0 ; j < triagPoints_.size() ; j++){
		if (pointMaxCellLID_[j] < 0){
			SUNDANCE_MSG3( verb , " TriangleSurf3D::computeMaxCellLIDs NO CELL j=" << j << " triagPoints_[j]=" << triagPoints_[j] );
		}
	}
}

Expr TriangleSurf3D::getParams() const
{
	//todo: return the points of the polygon, only later for the optimization
	// if necessary at all ...
	return Expr(List(0.0));
}


bool TriangleSurf3D::shortDistanceCalculation(const Point& p, double &res) const{

	double epsx = 1e-2*(maxX_ - minX_);
	double epsy = 1e-2*(maxY_ - minY_);
	double epsz = 1e-2*(maxZ_ - minZ_);

	//SUNDANCE_MSG3( 6 , "shortDistanceCalculation ,  p = " << p );
	//SUNDANCE_MSG3( 6 , "[ " << minX_ << "," << maxX_ << " ] " << "[ " << minY_ << "," << maxY_ << " ] " << "[ " << minZ_ << "," << maxZ_ << " ] ");
	if ( (p[0] < minX_- epsx ) || (p[0] > maxX_ + epsx) ||
		 (p[1] < minY_- epsy ) || (p[1] > maxY_ + epsy) ||
		 (p[2] < minZ_- epsz ) || (p[2] > maxZ_ + epsz) ) {
        // if the point is outside then just set the distance to a positive high number
		res = +1;
		return true;
	}

	return false;
}

double TriangleSurf3D::curveEquation_intern(const Point& evalPoint) const
{

	int verb = 0;
	double dist = 1e+100 , distTmp , distTmp1 , sign_tmp , tseg = 0.0;

	// make a shorter version of distance calculation
	SUNDANCE_MSG3( verb , "TriangleSurf3D::curveEquation_intern ,  evalPoint = " << evalPoint );
	if ( shortDistanceCalculation(evalPoint, dist) == true ) { return dist; }

	Point p0 , p1, p2 , n , intP , tmpP1 , tmpP0 = evalPoint , P;
	int spaceTreeCellID = -1;
	bool isIntersected , insideTriangle;
	const Array<int>& triangleIDs = getSpaceTree( evalPoint , spaceTreeCellID);

	for (int t = 0; t < triangleIDs.size() ; t++){
		//SUNDANCE_MSG3( verb , "triangleIDs[t]=" << triangleIDs[t] << " , t1=" << triagIndexes_[3*triangleIDs[t]] <<
		//		" , t2 = " << triagIndexes_[3*triangleIDs[t]+1] << " t3 = " << triagIndexes_[3*triangleIDs[t]+2]);
		p0 = triagPoints_[triagIndexes_[3*triangleIDs[t]]];
		p1 = triagPoints_[triagIndexes_[3*triangleIDs[t]+1]];
		p2 = triagPoints_[triagIndexes_[3*triangleIDs[t]+2]];
		//SUNDANCE_MSG3( verb , "p0=" << p0 << " , p1=" << p1 << " , p2=" << p2 );
		n = cross( p1-p0 , p2 - p0 );
		n = (1/::sqrt(n*n)) * n;
		tmpP1 = tmpP0 + n;
       	//SUNDANCE_MSG3( verb , "p0=" << p0 << " , p1=" << p1 << " , p2=" << p2 <<" , tmpP0=" << tmpP0 << ", n = " << n << " tmpP1 = " << tmpP1);
		// get the intersection point
		intP = triangleLineIntersect( p0 , p1 , p2 , tmpP0 , tmpP1 , isIntersected , tseg);
		if (isIntersected) {
       	   //SUNDANCE_MSG3( verb , "intP = " << intP );
		   // since the directions normal there should be always an intersection point
		   insideTriangle = pointInTriangle( p0 , p1, p2 , intP );
		   if (!insideTriangle) {
			  // todo: alternatively one could just take the average point and then the difference or other less computing intensive variants
			  projectPointToLine(p0 , p1 , tmpP0 , P ); tmpP1 = P; distTmp = distTmp1 = ::sqrt( (P-tmpP0)*(P-tmpP0));
			  projectPointToLine(p1 , p2 , tmpP0 , P ); distTmp1 = ::sqrt( (P-tmpP0)*(P-tmpP0));
			  if ( distTmp > distTmp1 )
			  { distTmp = distTmp1; tmpP1 = P; }
			  projectPointToLine(p0 , p2 , tmpP0 , P ); distTmp1 = ::sqrt( (P-tmpP0)*(P-tmpP0));
			  if ( distTmp > distTmp1 )
			  { distTmp = distTmp1; tmpP1 = P; }
			  // we take the projected point and take the distance from this point
			  distTmp = ::sqrt( (tmpP1-tmpP0)*(tmpP1-tmpP0) );
			  //SUNDANCE_MSG3( verb , "outside triag tmpP1 = " << tmpP1 );
		   }
		   else
		   {
			  distTmp = ::sqrt( (intP-tmpP0)*(intP-tmpP0) );
			  //SUNDANCE_MSG3( verb , "inside triag tmpP1 = " << tmpP1 );
		   }
      	   sign_tmp = (tmpP0-intP)*n;
       	   // if the vector product is positive then the point is inside
       	   if ( sign_tmp > 0.0){
       		  distTmp = -distTmp;
       	   }
       	   //SUNDANCE_MSG3( verb , "curveEquation() dist_tmp=" << distTmp << "  sign_tmp= " << sign_tmp);
       	   // we store the absolute minimal distance
       	   if (::fabs(dist) > ::fabs(distTmp)){
       		  //SUNDANCE_MSG3( verb , "curveEquation() store distance " << distTmp << " prev dist: " << dist << " , sign_tmp:" << sign_tmp << " , n=" << n);
       		  //SUNDANCE_MSG3( verb , "p0=" << p0 << " , p1=" << p1 << " , p2=" << p2 <<" , tmpP0=" << tmpP0 << ", n = " << n << " tmpP1 = " << tmpP1);
       		  dist = distTmp;
       	   }
		}
	}

    SUNDANCE_MSG3( verb , "curveEquation() valPoint=" << evalPoint << " return distance = " << dist);
	return dist;
}

bool TriangleSurf3D::shortIntersectCalculation(const Point& st, const Point& end) const{

	double epsx = 1e-2*(maxX_ - minX_);
	double epsy = 1e-2*(maxY_ - minY_);
	double epsz = 1e-2*(maxZ_ - minZ_);

	bool st_outside = ( (st[0] < minX_- epsx ) || (st[0] > maxX_ + epsx) ||
			            (st[1] < minY_- epsy ) || (st[1] > maxY_ + epsy) ||
			            (st[2] < minZ_- epsz ) || (st[2] > maxZ_ + epsz) );
	bool end_outside = ( (end[0] < minX_- epsx ) || (end[0] > maxX_ + epsx) ||
			            (end[1] < minY_- epsy ) || (end[1] > maxY_ + epsy) ||
			            (end[2] < minZ_- epsz ) || (end[2] > maxZ_ + epsz) );

	// if both points are outside then no intersection is possible
	// specially because we have only lines which are parallel to axes
	if ( st_outside && end_outside ) {
		return true;
	}
	return false;
}


void TriangleSurf3D::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
		Array<Point>& result) const
{
	int verb = 0;
	result.resize(0); nrPoints = 0;
	double eps = 1e-8 , tseg;

    //SUNDANCE_MSG3( verb , "returnIntersectPoints start=" << start << " , end= " << end );
	// make a shorter version of distance calculation
	if ( shortIntersectCalculation(start,end) == true ) {
		return;
	}

	Point p0 , p1, p2 , n , intP , tmpP0 = start , tmpP1 = end;
	int spaceTreeCellIDstart = -1 , spaceTreeCellIDend = -1;
	bool isIntersected , insideTriangle;
	const Array<int>& triangleIDStart = getSpaceTree( start , spaceTreeCellIDstart);
	const Array<int>& triangleIDEnd = getSpaceTree( end , spaceTreeCellIDend);

	// loop over each triangle and look for the intersection points
	for (int t = 0; t < triangleIDStart.size() ; t++)
	{
		p0 = triagPoints_[triagIndexes_[3*triangleIDStart[t]]];
		p1 = triagPoints_[triagIndexes_[3*triangleIDStart[t]+1]];
		p2 = triagPoints_[triagIndexes_[3*triangleIDStart[t]+2]];
		// get the intersection point
		intP = triangleLineIntersect( p0 , p1 , p2 , tmpP0 , tmpP1 , isIntersected , tseg);
		//SUNDANCE_MSG3( verb , " intersection intP =" << intP << " , isIntersected:" << isIntersected);
		if (isIntersected) {
			//SUNDANCE_MSG3( verb , " int P , intP =" << intP );
			insideTriangle = pointInTriangle( p0 , p1, p2 , intP );
			//if (insideTriangle) { SUNDANCE_MSG3( verb , "intP=" << intP ); }
			// test if the intersection point is inside the start-end line
			insideTriangle = insideTriangle && ( (intP[0] <= end[0]+eps) && (intP[0] >= start[0]-eps) );
			insideTriangle = insideTriangle && ( (intP[1] <= end[1]+eps) && (intP[1] >= start[1]-eps) );
			insideTriangle = insideTriangle && ( (intP[2] <= end[2]+eps) && (intP[2] >= start[2]-eps) );
			if (insideTriangle){
				// add to the list
				double minDist = 1e+100, distP;
				for (int pp = 0 ; pp < nrPoints ; pp++ ){
					distP = ::sqrt((intP-result[pp])*(intP-result[pp]));
					minDist = ::fmin(minDist,distP);
				}
				// ---
				if ( minDist > 1e-10 )  {
				   result.resize(nrPoints+1);
				   result[nrPoints] = intP;
				   //SUNDANCE_MSG3( verb , "result[" << nrPoints<< "]=" << result[nrPoints] );
				   nrPoints++;
				}
			}
		}
	}
	if (spaceTreeCellIDend != spaceTreeCellIDstart){
		// if the space tree cells are different then loop through there triangles
		for (int t = 0; t < triangleIDEnd.size() ; t++)
		{
			p0 = triagPoints_[triagIndexes_[3*triangleIDEnd[t]]];
			p1 = triagPoints_[triagIndexes_[3*triangleIDEnd[t]+1]];
			p2 = triagPoints_[triagIndexes_[3*triangleIDEnd[t]+2]];
			// get the intersection point
			intP = triangleLineIntersect( p0 , p1 , p2 , tmpP0 , tmpP1 , isIntersected , tseg);
			if (isIntersected) {
				insideTriangle = pointInTriangle( p0 , p1, p2 , intP );
				// test if the intersection point is inside the start-end line
				insideTriangle = insideTriangle && ( (intP[0] <= end[0]+eps) && (intP[0] >= start[0]-eps) );
				insideTriangle = insideTriangle && ( (intP[1] <= end[1]+eps) && (intP[1] >= start[1]-eps) );
				insideTriangle = insideTriangle && ( (intP[2] <= end[2]+eps) && (intP[2] >= start[2]-eps) );
				if (insideTriangle) {
					double minDist = 1e+100, distP;
					for (int pp = 0 ; pp < nrPoints ; pp++ ){
						distP = ::sqrt((intP-result[pp])*(intP-result[pp]));
						minDist = ::fmin(minDist,distP);
					}
					// ---
					if ( minDist > 1e-10 )  {
					   result.resize(nrPoints+1);
					   result[nrPoints] = intP;
					   //SUNDANCE_MSG3( verb , "result[" << nrPoints<< "]=" << result[nrPoints] );
					   nrPoints++;
					}
				}
			}
		}
	}
	//SUNDANCE_MSG3( 4 , "returnIntersectPoints nrPoints=" << nrPoints << " , start=" << start << " , end= " << end );
}


void TriangleSurf3D::returnIntersect(const Point& start, const Point& end,
		int& nrPoints, Array<double>& result) const
{
	Array<Point> t;
	returnIntersectPoints(start, end, nrPoints, t);

	result.resize(nrPoints);

	// Return coordinates instead of t values
	for (int i = 0; i < nrPoints; i++)
	{
		Point tmp( end[0]-t[i][0]/(end[0]-start[0]) , end[1]-t[i][1]/(end[1]-start[1]) , end[2]-t[i][2]/(end[2]-start[2]));
		result[i] =  sqrt( tmp*tmp );
	}
}


void TriangleSurf3D::writeToVTK(const std::string& filename) const
{
     std::ofstream myfile;
     myfile.open(filename.c_str());
     myfile << "# vtk DataFile Version 2.0 \n";
     myfile << "Generated by Sundance::TriangleSurf3D Author: Janos Benk \n";
     myfile << "ASCII\n";
     myfile << "\n";
     myfile << "DATASET UNSTRUCTURED_GRID\n";
     myfile << "POINTS "<< triagPoints_.size() << " float \n";
     myfile << "\n";
     for (int ii = 0 ; ii < triagPoints_.size() ; ii++){
    	 myfile << triagPoints_[ii][0] << " " << triagPoints_[ii][1] << " " << triagPoints_[ii][2] << "\n";
     }
     myfile << "\n \n";
     myfile << "CELLS " << triagIndexes_.size()/3 << " " << 4*(triagIndexes_.size()/3) << "\n";
     myfile << "\n";
     for (int ii = 0 ; ii < triagIndexes_.size()-1 ; ii = ii + 3){
    	 myfile << "3 "<< triagIndexes_[ii] << " " << triagIndexes_[ii+1] << " " << triagIndexes_[ii+2] << "\n";
     }
     myfile << "\n \n";
     myfile << "CELL_TYPES " << triagIndexes_.size()/3 << "\n";
     myfile << "\n";
     for (int ii = 0 ; ii < triagIndexes_.size()/3 ; ii++){
    	 myfile << "5\n";
     }
     myfile << "\n";

     if ( nrScalarField_ > 0) {
    	 myfile << "POINT_DATA " << triagPoints_.size() << "\n";
     }

     // first plot everything in a scalar field
     for (int s = 0 ; s < nrScalarField_ ; s++)
     {
    	 // plot each scalar field as a scalar field
    	 myfile << "SCALARS " << triagSurfSpaceNames_[s] << " float \n";
    	 myfile << "LOOKUP_TABLE default \n";
    	 for (int ii = 0 ; ii < triagSurfSpaceValues_[s].size() ; ii++)
    	 {
    		 myfile << triagSurfSpaceValues_[s][ii] << " \n";
    	 }
         myfile << "\n";
     }

     // since this is in 3D we could pair 2X to plot these 3 pairs as vectors
     // then plot 3 by the the scalar fields as vectors
     if ( nrScalarField_ > 2 )
     {
         for (int s = 0 ; s < nrScalarField_ ; s = s + 3)
         {
        	 myfile << "VECTORS " << triagSurfSpaceNames_[s] << triagSurfSpaceNames_[s+1] << triagSurfSpaceNames_[s+2] << "_vect float \n";
        	 for (int ii = 0 ; ii < triagSurfSpaceValues_[s].size() ; ii++)
        	 {
        		 myfile << triagSurfSpaceValues_[s][ii] <<  " " <<  triagSurfSpaceValues_[s+1][ii] << " " <<  triagSurfSpaceValues_[s+2][ii] << "\n";
        	 }
             myfile << "\n";
         }
     }

     myfile.close();
}

void TriangleSurf3D::addEvaluationPointValues(const Mesh& mesh ,
		int maxCellLID , int nQuad ,
		const double* coeffPtr ,
		const Array<Point>& quadPts) const{
	// - only when we want to record the information only then go into the routine
	// - look for the points which are in the same cell
	// - for the found points and for the chosen polygon space add the mean value of the nearest point

	//SUNDANCE_MSG3( 5 , "addEvaluationPointValues() isRecording_= " << isRecording_ << " recordedScalarField_="<< recordedScalarField_ );
	if (isRecording_)
	{ // only if the flag is set, then move to the recording
		Array<int> pointsToLookAt(0);
		int nrPointsLookAt = 0;
		int nrQuadPoint = quadPts.size();
		int verb = 0;

		SUNDANCE_MSG3( verb , "addEvaluationPointValues() maxCellLID= " << maxCellLID << " nQuad="<< nQuad );

		// if the polygon has no Mesh then throw an error
		TEUCHOS_TEST_FOR_EXCEPTION(  hasMesh_ == false , std::runtime_error,
			      " TriangleSurf3D::addEvaluationPointValues , Polygon must have a valid mesh");

		// look for the impacted points from the polygon
		// look for the impacted points from the polygon
		for (int ii = 0 ; ii < pointMaxCellLID_.size() ; ii++){
			// test if the polygon point "ii" is inside the "maxCellLID" cell
			if ( maxCellLID == pointMaxCellLID_[ii] ){
				pointsToLookAt.push_back( ii );
				nrPointsLookAt = nrPointsLookAt + 1;
			}
		}

		// transform all the points into real coordinates
		// then look for the two nearest one, and take the average one
		// we could also take the nearest one, but the average one
		Array<Point> realCoordPoints(quadPts.size());
		Array<Point> cellsPoint(8);
		int tmp;
		cellsPoint[0] = mesh_->nodePosition( mesh_->facetLID( mesh_->spatialDim() ,maxCellLID,0,0,tmp) );
		cellsPoint[7] = mesh_->nodePosition( mesh_->facetLID( mesh_->spatialDim() ,maxCellLID,0,7,tmp) );
		for (int p = 0 ; p < nrQuadPoint ; p++){
			// calculate the real coordinated
			realCoordPoints[p] = cellsPoint[0] + Point( (cellsPoint[7][0] - cellsPoint[0][0])*quadPts[p][0]  ,
					(cellsPoint[7][1] - cellsPoint[0][1])*quadPts[p][1] , (cellsPoint[7][2] - cellsPoint[0][2])*quadPts[p][2] );
			//SUNDANCE_MSG3( verb , "transformed reference coords= " << quadPts[p] << " to real coords="<< realCoordPoints[p] );
		}

		// for each polygon point inside this cell
		for ( int polyP = 0 ; polyP < nrPointsLookAt ; polyP++)
		{
			// for this polygon point look for the 2 nearest point
			int firstI = -1 , secondI = -1;
			double dist = 1e+100 , dist_tmp =0.0;
			// look for the nearest point
			for( int p = 0 ; p < nrQuadPoint ; p++){
				dist_tmp = (realCoordPoints[p] - triagPoints_[pointsToLookAt[polyP]])*(realCoordPoints[p] - triagPoints_[pointsToLookAt[polyP]]);
				if ( dist_tmp < dist ) {
					dist = dist_tmp;
					firstI = p;
				}
			}
			// look for the second nearest point
			dist = 1e+100;
			for( int p = 0 ; p < nrQuadPoint ; p++){
				dist_tmp = (realCoordPoints[p] - triagPoints_[pointsToLookAt[polyP]])*(realCoordPoints[p] - triagPoints_[pointsToLookAt[polyP]]);
				if ( (dist_tmp < dist) && ( p != firstI ) ) {
					dist = dist_tmp;
					secondI = p;
				}
			}
			//having the two nearest point we add the average value
			SUNDANCE_MSG3( verb , "pointsToLookAt["<<polyP<<"]=" << pointsToLookAt[polyP] << " , v1="
					<< coeffPtr[firstI] << " , v2="<< coeffPtr[secondI]);
			SUNDANCE_MSG3( verb , " Set Point =" << triagPoints_[pointsToLookAt[polyP]] << " , p0="
					<< realCoordPoints[firstI] << " , p1="<< realCoordPoints[secondI]);
			// here we just set the value since one point should be set only once
			triagSurfSpaceValues_[recordedScalarField_][pointsToLookAt[polyP]] = 0.5*(coeffPtr[firstI] + coeffPtr[secondI]);
		}
	}
}


void TriangleSurf3D::setSpaceValues(const FunctionalEvaluatorBase& scalarFunctional , int fieldIndex ){

	int verb = 0;
	SUNDANCE_MSG3( verb , "setSpaceValues() fieldIndex= " << fieldIndex );

	// signaling that we are starting recording values
	isRecording_ = true;
	recordedScalarField_ = fieldIndex;
	SUNDANCE_MSG3( verb , " set values to zero " );
	// set all values to zero because now the recorded values will be added to the vector elements
	for (int ii = 0 ; ii < triagSurfSpaceValues_[fieldIndex].size() ; ii++ ){
		triagSurfSpaceValues_[fieldIndex][ii] = 0.0;
	}

	SUNDANCE_MSG3( verb , " BEFORE evaluate the FunctionalEvaluator object curveID:" );
	// --- evaluate the scalar function ---
	// this function has to be a curve integral of a scalar function
	// and this polygon has to be the argument for the curve integral
	scalarFunctional.evaluate();

	SUNDANCE_MSG3( verb , " AFTER evaluate the FunctionalEvaluator object  curveID:" );

	// in parallal case make a vector reduction
	if ( this->mesh_ != 0 ){
		  if (mesh_->comm().getNProc() > 1)
		  {
			  SUNDANCE_MSG3( verb , " make array vector reduction for parallel case " );
			  // create one temporary vector
			  Array<double> reduceVector( triagSurfSpaceValues_[fieldIndex].size() , 0.0 );
			  // call the MPI function to reduce the vector
			  SUNDANCE_MSG3( verb , " call the reduce function " );
			  mesh_->comm().allReduce((void*) &(triagSurfSpaceValues_[fieldIndex][0]), (void*) &(reduceVector[0]),
          triagSurfSpaceValues_[fieldIndex].size() ,MPIDataType::doubleType(), MPIOp::sumOp());
			  // copy the reduced vector back
			  SUNDANCE_MSG3( verb , " copy the reduced values back " );
			  for (int ii = 0 ; ii < triagSurfSpaceValues_[fieldIndex].size() ; ii++ ){
				  triagSurfSpaceValues_[fieldIndex][ii] = reduceVector[ii];
			  }
		  }
	}

	// if the twin polygon exists then copy also there the value
	if (twinTriangleSurf_ != 0){
		SUNDANCE_MSG3( verb , " update twin polygon " );
		// copy the recorded values to the twin polygon
		for (int ii = 0 ; ii < triagSurfSpaceValues_[fieldIndex].size() ; ii++ ){
			twinTriangleSurf_->triagSurfSpaceValues_[fieldIndex][ii] = this->triagSurfSpaceValues_[fieldIndex][ii];
		}
	}

	// we set the flags so that we stop recording
	recordedScalarField_ = -1;
	isRecording_ = false;
}

void TriangleSurf3D::eval0(const double* vars, double* f , int scalarFieldIndex ) const {

	  Point inCoord(vars[0],vars[1],vars[2]);
	  int verb = 0;

	  // first look for the triangle which is the closest
	  Point p0 , p1, p2 , n , intP , tmpP1 , tmpP0 = inCoord , P = inCoord , PP = inCoord;
	  double distTmp = 1e+100 , dist = 1e+100 , tseg = 0.0;
	  int spaceTreeCellID = -1 , triagIndex = -1;
	  bool isIntersected , insideTriangle;
	  const Array<int>& triangleIDs = getSpaceTree( inCoord , spaceTreeCellID);
	  for (int t = 0; t < triangleIDs.size() ; t++){
			p0 = triagPoints_[triagIndexes_[3*triangleIDs[t]]];
			p1 = triagPoints_[triagIndexes_[3*triangleIDs[t]+1]];
			p2 = triagPoints_[triagIndexes_[3*triangleIDs[t]+2]];
			n = cross( p1-p0 , p2 - p0 );
			n = ::sqrt(n*n) * n;
			tmpP1 = tmpP0 + n;
			//SUNDANCE_MSG3( verb , "p0=" << p0 << " , p1=" << p1 << " , p2=" << p2 << " , n=" << n);
			// get the intersection point
			intP = triangleLineIntersect( p0 , p1 , p2 , tmpP0 , tmpP1 , isIntersected , tseg);
			// since the directions normal there should be always an intersection point
			insideTriangle = pointInTriangle( p0 , p1, p2 , intP );
			if (!insideTriangle) {
				// alternatively one could just take the average point and then the difference or other less computing intensive variants
				P = (1/3)*(p0+p1+p2);
				// we take the projected point and take the distance from this point
				distTmp = ::sqrt( (P-tmpP0)*(P-tmpP0) );
			}
			else
			{
				distTmp = ::sqrt( (intP-tmpP0)*(intP-tmpP0) );
				P = intP;
			}

	       	//SUNDANCE_MSG3( verb , "curveEquation() distTmp=" << distTmp );
	       	// we store the absolute minimal distance
	       	if (::fabs(dist) > ::fabs(distTmp)){
	       		//SUNDANCE_MSG3( verb , "curveEquation() store distance " << distTmp << " prev dist: " << dist);
	       		dist = distTmp;
	       		triagIndex = t;
	       		PP = P;
	       	}
	  }
	  // triagIndex has the nearest triangle
      p0 = triagPoints_[triagIndexes_[3*triagIndex]];
	  p1 = triagPoints_[triagIndexes_[3*triagIndex+1]];
	  p2 = triagPoints_[triagIndexes_[3*triagIndex+2]];
	  //SUNDANCE_MSG3( verb , " triagIndex = " << triagIndex << " , p0=" << p0 << " , p1=" << p1 << " , p2=" << p2 << " , n=" << n);
	  double d0 = ::sqrt((PP-p0)*(PP-p0));
	  double d1 = ::sqrt((PP-p1)*(PP-p1));
	  double d2 = ::sqrt((PP-p2)*(PP-p2));
	  double sumd = d0+d1+d2;

	  // TODO: check this
	  f[0] = 0.5*((d1+d2)/sumd) * triagSurfSpaceValues_[scalarFieldIndex][triagIndexes_[3*triagIndex]];
	  f[0] += 0.5*((d0+d2)/sumd) * triagSurfSpaceValues_[scalarFieldIndex][triagIndexes_[3*triagIndex+1]];
	  f[0] += 0.5*((d0+d1)/sumd) * triagSurfSpaceValues_[scalarFieldIndex][triagIndexes_[3*triagIndex+2]];
	  SUNDANCE_MSG3( verb , "eval0() p= " << inCoord );
	  //SUNDANCE_MSG3( verb , "v0 = " << triagSurfSpaceValues_[scalarFieldIndex][triagIndexes_[3*triagIndex]] <<
		//	             " , v1 =" << triagSurfSpaceValues_[scalarFieldIndex][triagIndexes_[3*triagIndex+1]] <<
		//	             " , v2 =" << triagSurfSpaceValues_[scalarFieldIndex][triagIndexes_[3*triagIndex+2]] );
	  SUNDANCE_MSG3( verb , "eval0() ret = " << f[0] << " d0=" << d0 << " d1=" << d1 << " d2=" << d2);
}

TriangleSurf3D* TriangleSurf3D::createTwinTriangleSurf(double shiftX , double shiftY , double shiftZ ,
		 double scaleX , double scaleY , double scaleZ ){

	Array<Point> points(triagPoints_.size());

	// make the shift and the scaling
	for (int j = 0 ; j < triagPoints_.size() ; j++ ){
		points[j] = Point( shiftX + scaleX*triagPoints_[j][0] , shiftY + scaleY*triagPoints_[j][1] , shiftZ + scaleZ*triagPoints_[j][2]);
	}

	// create the polygon
	TriangleSurf3D* surf = new TriangleSurf3D( points , this->triagIndexes_ ,this->_alpha1 , this->_alpha2 , flipDomains_ );

	// add all the scalar fields in the same order
	for (int j = 0 ; j < nrScalarField_ ; j++){
		surf->addNewScalarField( triagSurfSpaceNames_[j], 0.0 );
	}

	// store the twin polygon mutually
	twinTriangleSurf_ = surf;
	surf->twinTriangleSurf_ = this;

	return surf;
}

TriangleSurf3D* TriangleSurf3D::importGTSSurface( const std::string& filename , double a1 , double a2 , Point innerP) {

	std::string str_tmp;
	std::ifstream myfile;
	std::vector<std::string> elems;
	//int verb = 0;

	char *pEnd , *pEnd1;
	double d1 , d2, d3;
	int    i1 , i2, i3;
	int nr_points = 0 , nr_lines = 0 , nr_triags = 0;
	Array<Point> triagPoints;
	Array<int> lineIndex;
	Array<int> triagIndex;

	myfile.open( filename.c_str() , std::ios::in );

	// the first line which continas the number of points lines and triangles
	getline( myfile , str_tmp );
	nr_points = round(strtod (str_tmp.c_str() , &pEnd));
	nr_lines = round(strtod (pEnd , &pEnd1));
	nr_triags = round(strtod (pEnd1 , NULL));

	// each line of the file should contain a point
	triagPoints.resize(nr_points);
    for (int p = 0 ; p < nr_points ; p++)
    {
    	// read the line
    	getline( myfile , str_tmp );
    	d1 = strtod (str_tmp.c_str() , &pEnd);
    	d2 = strtod (pEnd , &pEnd1);
    	d3 = strtod (pEnd1 , NULL);
		// we add the point to the polygon
		Point tmpP(d1,d2,d3);
		triagPoints[p] = tmpP;
		//SUNDANCE_MSG3( verb , " point[" << p << "] = " << triagPoints[p]);
    }
    lineIndex.resize(2*nr_lines);
    for (int p = 0 ; p < nr_lines ; p++)
    {
    	// read the line
    	getline( myfile , str_tmp );
    	i1 = round(strtod (str_tmp.c_str() , &pEnd));
    	i2 = round(strtod (pEnd , &pEnd1));
    	lineIndex[2*p] = i1-1;  lineIndex[2*p+1] = i2-1;
    }
    triagIndex.resize(3*nr_triags);
    for (int p = 0 ; p < nr_triags ; p++)
    {
    	// read the line
    	getline( myfile , str_tmp );
    	i1 = round(strtod (str_tmp.c_str() , &pEnd)) - 1 ;
    	i2 = round(strtod (pEnd , &pEnd1)) - 1;
    	i3 = round(strtod (pEnd1 , NULL)) - 1;
    	//SUNDANCE_MSG3( verb , " i1 = " << i1 << " i2 = " << i2 << " i3 = " << i3);
    	// there are the index from the line, find the corresponding
    	if ( (lineIndex[2*i1] == lineIndex[2*i2]) || (lineIndex[2*i1] == lineIndex[2*i2+1]) ) {
           triagIndex[3*p] = lineIndex[2*i1];  triagIndex[3*p+1] = lineIndex[2*i1+1];
    	   if ( lineIndex[2*i2] != lineIndex[2*i1]){
    		  triagIndex[3*p+2] = triagIndex[3*p+1];  triagIndex[3*p+1] = lineIndex[2*i2];
    	   }
    	   else{
    		  triagIndex[3*p+2] = triagIndex[3*p+1];  triagIndex[3*p+1] = lineIndex[2*i2+1];
    	   }
    	} else { // this means that the second point from the first line is in the second line
            triagIndex[3*p] = lineIndex[2*i1+1];  triagIndex[3*p+1] = lineIndex[2*i1];
     	   if ( lineIndex[2*i2] != lineIndex[2*i1+1]){
     		  triagIndex[3*p+2] = triagIndex[3*p+1];  triagIndex[3*p+1] = lineIndex[2*i2];
     	   }
     	   else{
     		  triagIndex[3*p+2] = triagIndex[3*p+1];  triagIndex[3*p+1] = lineIndex[2*i2+1];
     	   }
    	}

    	// this is only for
    	Point p1 = triagPoints[triagIndex[3*p]];
    	Point p2 = triagPoints[triagIndex[3*p+1]];
    	Point p3 = triagPoints[triagIndex[3*p+2]];
    	Point n = cross( p2 - p1 , p3 - p1 );
    	Point m = (1.0/3.0)*(p1+p2+p3);
    	double dOrig = (innerP-m)*(innerP-m);
    	double dNew = (innerP-(m+1e-3*n))*(innerP-(m+1e-3*n));

    	if ( dOrig < dNew){
    		// this means the normal vector points outwards, so we have to flip
    		int tmp_swap = triagIndex[3*p+1];
    		triagIndex[3*p+1] = triagIndex[3*p+2];
    		triagIndex[3*p+2] = tmp_swap;
    	}
		//SUNDANCE_MSG3( verb , " tirag indexes [" << triagIndex[3*p] << " , " << triagIndex[3*p+1] << " , " << triagIndex[3*p+2] << "]");
    }
	// close the file
	myfile.close();
    // create the object and return
	return (new TriangleSurf3D( triagPoints , triagIndex , a1 , a2 ));
}

TriangleSurf3D* TriangleSurf3D::importSTLSurface( const std::string& filename , double a1 , double a2 , Point innerP){
	std::string str_tmp;
	std::ifstream myfile;
	int verb = 0;

	char *pEnd , *pEnd1;
	double d1 , d2 , d3;
	int index[3];
	int index_tmp = 0 , tmpI , nrTriag = 0;
	Array<Point> triagPoints(0);
	Array<int> triagIndex(0);

	SUNDANCE_MSG3( verb , "TriangleSurf3D::importSTLSurface");

	myfile.open( filename.c_str() , std::ios::in );

	getline( myfile , str_tmp );
	getline( myfile , str_tmp );
	getline( myfile , str_tmp );

	// each line of the file should contain a point
	while (!myfile.eof())
	{
		for ( int d = 0 ; d < 3 ; d++){
			getline( myfile , str_tmp );
			int tm = 0;
			while (str_tmp[tm] != 'x') {tm++;}
	    	d1 = strtod(&(str_tmp.c_str()[tm+1]) , &pEnd);
	    	d2 = strtod(pEnd , &pEnd1);
	    	d3 = strtod(pEnd1 , NULL);
	    	Point poin(d1,d2,d3);
	    	SUNDANCE_MSG3( verb , "Point read in " << d1 << " , " << d2 << " , " << d3);
	    	tmpI = -1;
	    	for( int p = 0 ; p < index_tmp ; p++){
	    		if ( ((triagPoints[p] - poin)*(triagPoints[p] - poin)) < 1e-15){
	    			tmpI = p;
	    		}
	    	}
	    	// look for this point index
	    	if ( tmpI > -1){
	    		index[d] = tmpI;
	    	} else {
	    		triagPoints.resize(index_tmp+1);
		    	triagPoints[index_tmp] = poin;
		    	SUNDANCE_MSG3( verb , "New Point:" << triagPoints[index_tmp] );
		    	index[d] = index_tmp;
		    	index_tmp++;
	    	}
		}
		SUNDANCE_MSG3( verb , "TriangleSurf3D::importSTLSurface i0:" << index[0] << " , i1:" << index[1] << " , i2:" << index[2]);
		SUNDANCE_MSG3( verb , "TriangleSurf3D::importSTLSurface Points t0:" << triagPoints[index[0]] << " , t1:" << triagPoints[index[1]]
		                       << " , t2:" << triagPoints[index[2]]);
		triagIndex.resize(3*(nrTriag+1));
		triagIndex[3*nrTriag] = index[0];
		triagIndex[3*nrTriag+1] = index[1];
		triagIndex[3*nrTriag+2] = index[2];
		nrTriag++;

		int p = nrTriag -1;
    	Point p1 = triagPoints[triagIndex[3*p]];
    	Point p2 = triagPoints[triagIndex[3*p+1]];
    	Point p3 = triagPoints[triagIndex[3*p+2]];
    	Point n = cross( p2 - p1 , p3 - p1 );
    	Point m = (1.0/3.0)*(p1+p2+p3);
    	double dOrig = (innerP-m)*(innerP-m);
    	double dNew = (innerP-(m+1e-3*n))*(innerP-(m+1e-3*n));

    	if ( dOrig > dNew){
    		// this means the normal vector points inwards, so we have to flip
    		int tmp_swap = triagIndex[3*p+1];
    		triagIndex[3*p+1] = triagIndex[3*p+2];
    		triagIndex[3*p+2] = tmp_swap;
    	}

		getline( myfile , str_tmp );
		getline( myfile , str_tmp );
		getline( myfile , str_tmp );
		getline( myfile , str_tmp );
	}

	// close the file
	myfile.close();
	return (new TriangleSurf3D( triagPoints , triagIndex , a1 , a2 ));
}
