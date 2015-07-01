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

#include "SundancePolygon2D.hpp"
#include "SundancePoint.hpp"
#include "SundanceMesh.hpp"
#include "SundanceDefs.hpp"
#include "PlayaMPIComm.hpp"

#include <iostream>
#include <fstream>

using namespace Sundance;
using Playa::MPIOp;
using Playa::MPIDataType;

int Polygon2D::intersectionEdge_ = -1;

Polygon2D::Polygon2D(const Mesh& mesh , const Array<Point>& points , double a1, double a2, bool closedPolygon ,bool flipD ) :
	CurveBase(1, a2, a1 , flipD), hasMesh_(true), closedPolygon_(closedPolygon) , mesh_(&(mesh)) ,
	polygonSpaceValues_() , polygonSpaceNames_() , nrScalarField_(0) , isRecording_(false),twinPolygon_(0) // here we twist the alpha parameters for the polygon convention
{
	int verb = 0;
	// just store the input points
	polyPoints_.resize(points.size());
	SUNDANCE_MSG3( verb , "Polygon2D() Ctor nrPoints=" << points.size() );
	maxX_ = minX_ = points[0][0];
	maxY_ = minY_ = points[0][1];
	for (int j = 0 ; j < points.size() ; j++ ){
		polyPoints_[j] = points[j];
		SUNDANCE_MSG3( verb , " point[" << j << "] = " << polyPoints_[j]);
		maxX_ = (maxX_ > points[j][0]) ? maxX_ : points[j][0];
		minX_ = (minX_ < points[j][0]) ? minX_ : points[j][0];
		maxY_ = (maxY_ > points[j][1]) ? maxY_ : points[j][1];
		minY_ = (minY_ < points[j][1]) ? minY_ : points[j][1];
	}
    // get the maxCellsLID for each point
	computeMaxCellLIDs();
}

Polygon2D::Polygon2D(const Array<Point>& points , double a1, double a2, bool closedPolygon , bool flipD ) :
	CurveBase(1, a2, a1, flipD), hasMesh_(false), closedPolygon_(closedPolygon), mesh_(0),
	polygonSpaceValues_() , polygonSpaceNames_() , nrScalarField_(0) , isRecording_(false) ,twinPolygon_(0) // here we twist the alpha parameters for the polygon convention
{
	int verb = 0;
	// just store the input points
	polyPoints_.resize(points.size());
	SUNDANCE_MSG3( verb , "Polygon2D() Ctor nrPoints=" << points.size() );
	maxX_ = minX_ = points[0][0];
	maxY_ = minY_ = points[0][1];
	for (int j = 0 ; j < points.size() ; j++ ){
		polyPoints_[j] = points[j];
		SUNDANCE_MSG3( verb , " point[" << j << "] = " << polyPoints_[j]);
		maxX_ = (maxX_ > points[j][0]) ? maxX_ : points[j][0];
		minX_ = (minX_ < points[j][0]) ? minX_ : points[j][0];
		maxY_ = (maxY_ > points[j][1]) ? maxY_ : points[j][1];
		minY_ = (minY_ < points[j][1]) ? minY_ : points[j][1];
	}
    // get the maxCellsLID for each point
	computeMaxCellLIDs();
}

Polygon2D::Polygon2D(const Mesh& mesh , const std::string& filename , double a1, double a2, bool closedPolygon ,bool flipD ) :
CurveBase(1, a2, a1, flipD), hasMesh_(true), closedPolygon_(closedPolygon) , mesh_(&(mesh)),
polygonSpaceValues_() , polygonSpaceNames_() , nrScalarField_(0) , isRecording_(false), twinPolygon_(0) // here we twist the alpha parameters for the polygon convention
{
	std::string str_tmp;
	std::ifstream myfile;
	std::vector<std::string> elems;
	int verb = 6;

	char *pEnd;
	double d1 , d2;
	int index_tmp = 0;

	SUNDANCE_MSG3( verb , "Polygon2D() read from file ");

	myfile.open( filename.c_str() , std::ios::in );
	polyPoints_.resize(index_tmp);

	getline( myfile , str_tmp );
	d1 = strtod (str_tmp.c_str() , &pEnd);
	d2 = strtod (pEnd , NULL);
	// each line of the file should contain a point
	while (!myfile.eof())
	{
		// we add the point to the polygon
		index_tmp = polyPoints_.size();
		polyPoints_.resize(index_tmp+1);
		Point tmpP(d1,d2);
		polyPoints_[index_tmp] = tmpP;
		SUNDANCE_MSG3( verb , " point[" << index_tmp << "] = " << polyPoints_[index_tmp]);
		getline( myfile , str_tmp );
		//SUNDANCE_MSG3( verb , " next line = " << str_tmp << " line size = " << str_tmp.size());
		// if the line is less then 2 then break
		if (str_tmp.size() < 2 ) break;
		d1 = strtod (str_tmp.c_str() , &pEnd);
		d2 = strtod (pEnd , NULL);

		// look for the bounding box
		if (index_tmp == 0){
			maxX_ = minX_ = tmpP[0];
			maxY_ = minY_ = tmpP[1];
		}else{
			maxX_ = (maxX_ > tmpP[0]) ? maxX_ : tmpP[0];
			minX_ = (minX_ < tmpP[0]) ? minX_ : tmpP[0];
			maxY_ = (maxY_ > tmpP[1]) ? maxY_ : tmpP[1];
			minY_ = (minY_ < tmpP[1]) ? minY_ : tmpP[1];
		}
    }

	// close the file
	myfile.close();
	// compute the maxCells in which are the points
	computeMaxCellLIDs();
}

Polygon2D::Polygon2D( const std::string& filename , double a1, double a2, bool closedPolygon ,bool flipD ) :
CurveBase(1, a2, a1, flipD), hasMesh_(false), closedPolygon_(closedPolygon), mesh_(0),
polygonSpaceValues_() , polygonSpaceNames_() , nrScalarField_(0) , isRecording_(false), twinPolygon_(0)// here we twist the alpha parameters for the polygon convention
{
	std::string str_tmp;
	std::ifstream myfile;
	std::vector<std::string> elems;
	int verb = 6;

	char *pEnd;
	double d1 , d2;
	int index_tmp = 0;

	SUNDANCE_MSG3( verb , "Polygon2D() read from file ");

	myfile.open( filename.c_str() , std::ios::in );
	polyPoints_.resize(index_tmp);

	getline( myfile , str_tmp );
	d1 = strtod (str_tmp.c_str() , &pEnd);
	d2 = strtod (pEnd , NULL);
	// each line of the file should contain a point
	while (!myfile.eof())
	{
		// we add the point to the polygon
		index_tmp = polyPoints_.size();
		polyPoints_.resize(index_tmp+1);
		Point tmpP(d1,d2);
		polyPoints_[index_tmp] = tmpP;
		SUNDANCE_MSG3( verb , " point[" << index_tmp << "] = " << polyPoints_[index_tmp]);
		getline( myfile , str_tmp );
		//SUNDANCE_MSG3( verb , " next line = " << str_tmp << " line size = " << str_tmp.size());
		// if the line is less then 2 then break
		if (str_tmp.size() < 2 ) break;
		d1 = strtod (str_tmp.c_str() , &pEnd);
		d2 = strtod (pEnd , NULL);

		// look for the bounding box
		if (index_tmp == 0){
			maxX_ = minX_ = tmpP[0];
			maxY_ = minY_ = tmpP[1];
		}else{
			maxX_ = (maxX_ > tmpP[0]) ? maxX_ : tmpP[0];
			minX_ = (minX_ < tmpP[0]) ? minX_ : tmpP[0];
			maxY_ = (maxY_ > tmpP[1]) ? maxY_ : tmpP[1];
			minY_ = (minY_ < tmpP[1]) ? minY_ : tmpP[1];
		}
    }

	// close the file
	myfile.close();
	// compute the maxCells in which are the points
	computeMaxCellLIDs();
}

void Polygon2D::setMesh(const Mesh& mesh){

	mesh_ = &(mesh);
	hasMesh_ = true;
	// set the maxLID for each point of the polygon
	computeMaxCellLIDs();
}

void Polygon2D::computeMaxCellLIDs(){

	// test if we have a mesh, otherwise
	if (hasMesh_ == false) return;

	int meshDim = mesh_->spatialDim();
	int nrLID = mesh_->numCells(meshDim);
	//int verb = 6;

	//SUNDANCE_MSG3( verb , " Polygon2D::computeMaxCellLIDs " << polyPoints_.size() );

	// initially set everything to -1, so in parallel case we'll know which do not belong to this processor
	pointsMaxCellLID_.resize(polyPoints_.size());
	for (int j = 0 ; j < polyPoints_.size() ; j++){ pointsMaxCellLID_[j] = -1; }

	for (int cellLID = 0 ; cellLID < nrLID ; cellLID++){
		// set the LID first for all points to -1
		// run through each cell, and in each cell look for each point which is contained in that cell
		//
		switch (mesh_->cellType(meshDim)){
			case QuadCell:{
				// in the case of Quadcells we always have structured cells (Quadcells)
				int tmp;
				Point p0 = mesh_->nodePosition( mesh_->facetLID(meshDim,cellLID,0,0,tmp) );
				//Point p1 = mesh_.nodePosition( mesh_.facetLID(meshDim,cellLID,0,1,tmp) );
				//Point p2 = mesh_.nodePosition( mesh_.facetLID(meshDim,cellLID,0,2,tmp) );
				Point p3 = mesh_->nodePosition( mesh_->facetLID(meshDim,cellLID,0,3,tmp) );
				//SUNDANCE_MSG3( verb , " Polygon2D::computeMaxCellLIDs cellLID =" << cellLID << " ,p0:" << p0 << " ,p3:" << p3 );
				for (int j = 0 ; j < polyPoints_.size() ; j++)
				{
	               //SUNDANCE_MSG3( verb , " test j = " << j );
                   Point& tmp = polyPoints_[j];
                   if ( (tmp[0] >= p0[0]) &&  (tmp[0] <= p3[0]) &&
                		(tmp[1] >= p0[1]) &&  (tmp[1] <= p3[1]))
                   {
                	   // store the maxCellLID
                	   //SUNDANCE_MSG3( verb , " Polygon2D::computeMaxCellLIDs j=" << j << " maxCellLID=" << polyPoints_.size() );
                	   pointsMaxCellLID_[j] = cellLID;
                   }
				}
				break;
			}
			case TriangleCell:{
				// todo:
				break;
			}
			default: { /* throw error */}
		}
	} // - from the for loop
}

// update the state of the curve of the control point were changed
void Polygon2D::update() {
	maxX_ = minX_ = polyPoints_[0][0];
	maxY_ = minY_ = polyPoints_[0][1];
	// get the maximum and minimum values
	for (int j = 0 ; j < polyPoints_.size() ; j++ ){
		maxX_ = (maxX_ > polyPoints_[j][0]) ? maxX_ : polyPoints_[j][0];
		minX_ = (minX_ < polyPoints_[j][0]) ? minX_ : polyPoints_[j][0];
		maxY_ = (maxY_ > polyPoints_[j][1]) ? maxY_ : polyPoints_[j][1];
		minY_ = (minY_ < polyPoints_[j][1]) ? minY_ : polyPoints_[j][1];
	}
	// computes the cellLIDs
	computeMaxCellLIDs();
}

Expr Polygon2D::getParams() const
{
	//todo: return the points of the polygon, only later for the optimization
	// if necessary at all ...
	return Expr(List(0.0));
}


bool Polygon2D::shortDistanceCalculation(const Point& p, double &res) const{
	double px = minX_;
	double ox = maxX_ - minX_;
	double py = minY_;
	double oy = maxY_ - minY_;

	px = px + 0.5*ox;
	py = py + 0.5*oy;
	double r = (1.0 + 1e-9)*(::sqrt( (ox/2)*(ox/2) + (oy/2)*(oy/2) ));

	double d = ( (p[0]-px)*(p[0]-px) + (p[1]-py)*(p[1]-py) - r * r);

	if (d > 1e-14){
		res = d;
		return true;
	}else{
		return false;
	}
}

double Polygon2D::curveEquation_intern(const Point& evalPoint) const
{
	//int verb = 0;
	double dist = 1e+100;
	double dist_tmp = 0.0;
	double sign_tmp = 0.0;
	Point p0 , p1;
	double eps = 1e-14;


	// make a shorter version of distance calculation
	if ( shortDistanceCalculation(evalPoint, dist) == true ) { return dist; }


	// for each line segment
    for (int ii = 0 ; ii < polyPoints_.size() ; ii++)
    {
    	if (ii < polyPoints_.size() - 1)
    	{
        	p0 = polyPoints_[ii];
        	p1 = polyPoints_[ii+1];
    	}
    	else
    	{   // the last line from size-1 to 0
    		if (!closedPolygon_) continue; // only if the polygon is closed
        	p0 = polyPoints_[ii];
        	p1 = polyPoints_[0];
    	}
    	Point intP(0.0,0.0);

    	// this is the formula from the equation V*U = 0 and det|P-P0,P-P1| = 0 (from paper)
    	double a[4] = { p1[0]-p0[0] , p1[1]-p0[1] ,
    			        p0[1]-p1[1] , p1[0]-p0[0] };
    	double b[2] = { evalPoint[0]*(p1[0] - p0[0]) + evalPoint[1]*(p1[1] - p0[1]) ,
    			        p0[1]*(p1[0] - p0[0]) - p0[0]*(p1[1] - p0[1]) };

    	// simple Gauss elimination
    	if (::fabs(a[0]) < eps ){
        	double fakt = -a[0]/a[2];
        	intP[1] = (b[0] + fakt*b[1])/(a[1] + fakt*a[3]);
        	intP[0] = (b[1] - a[3]*intP[1]) / a[2];
    	}
    	else{
        	double fakt = -a[2]/a[0];
        	intP[1] = (b[1] + fakt*b[0])/(a[3] + fakt*a[1]);
        	intP[0] = (b[0] - a[1]*intP[1]) / a[0];
    	}
    	// now we have the intersection point
    	//SUNDANCE_MSG3( verb , " intP= " << intP );

    	if  (  (    ((intP[0] >= p0[0]-eps) && (intP[0] <= p1[0]+eps))
    		     || ((intP[0] <= p0[0]+eps) && (intP[0] >= p1[0]-eps)) )
    		&& (    ((intP[1] >= p0[1]-eps) && (intP[1] <= p1[1]+eps))
        		 || ((intP[1] <= p0[1]+eps) && (intP[1] >= p1[1]-eps)) )
			)
    	{
    		dist_tmp = 0.0;
    	}
    	else
    	{
    		double d1 = ::sqrt((intP-p0)*(intP-p0));
    		double d2 = ::sqrt((intP-p1)*(intP-p1));
    		dist_tmp = ( d1 < d2 )? d1 : d2;
    	}
   		//SUNDANCE_MSG3( verb , "curveEquation() evalPoint=" << evalPoint << " intP=" << intP <<
   		//		       " is between points" << p0 << " and " << p1 << " dist_tmp=" << dist_tmp );
       	dist_tmp = ::sqrt(dist_tmp*dist_tmp + (intP-evalPoint)*(intP-evalPoint) );
       	Point v1 = p1 - p0;
       	Point v2 = evalPoint - p0;
       	// v1 X v2
       	sign_tmp = v1[0]*v2[1] - v1[1]*v2[0];
       	// if the vector product is positive then the point is inside
       	if ( sign_tmp > 0.0){
      		dist_tmp = -dist_tmp;
       	}
       	//SUNDANCE_MSG3( verb , "curveEquation() dist_tmp=" << dist_tmp );
       	// we store the absolute minimal distance
       	if (::fabs(dist) > ::fabs(dist_tmp)){
       		//SUNDANCE_MSG3( verb , "curveEquation() store distance " << dist_tmp << " prev dist: " << dist);
       		dist = dist_tmp;
       	}

    }

    //SUNDANCE_MSG3( verb , "curveEquation() valPoint=" << evalPoint << " return distance = " << dist);
	return dist;
}


bool Polygon2D::shortIntersectCalculation(const Point& st, const Point& end) const{
	double px = minX_;
	double ox = maxX_ - minX_;
	double py = minY_;
	double oy = maxY_ - minY_;

	px = px + 0.5*ox;
	py = py + 0.5*oy;
	double r = (1.0 + 1e-9)*(::sqrt( (ox/2)*(ox/2) + (oy/2)*(oy/2) ));

	double d1 = ( (st[0]-px)*(st[0]-px) + (st[1]-py)*(st[1]-py) - r * r);
	double d2 = ( (end[0]-px)*(end[0]-px) + (end[1]-py)*(end[1]-py) - r * r);

	// if one is inside then we can exit here
	if ( (d1 < -1e-14) || ((d2 < -1e-14)) ) { return false; }

	// direction vector; line == start + t * vec
	Point vec = end - st;

	// center point of circle
	Point center( px, py );

	// intersection of circle and line leads to
	// norm(start + t * vec - center) == r
	// i.e. following coefficients of quadratic equation
	double a = vec * vec;
	double b = 2 * ((st - center) * vec);
	double c = ((st - center) * (st - center)) - r * r;

	double det = b * b - 4 * a * c;
	if (det < 0.0) { return false; } // no solution, both are outside and no intersection with the circle
	det = ::sqrt(det);

	// numerically more favorable than usual usage of solution formula
	int sign = b < 0.0 ? -1 : 1;
	double q = -0.5 * (b + sign * det);

	// first solution
	double t = c / q;

	// Only consider intersection points between 'start' and 'end'
	if (t >= 0.0 && t <= 1.0) return true;

	// second solution (only considered if different from first)
	if (det < 1.e-14) return false;
	t = q / a;
	if (t >= 0.0 && t <= 1.0) return true;

	return false;
}


void Polygon2D::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
		Array<Point>& result) const
{
	int verb = 0;
	nrPoints = 0;
	result.resize(nrPoints);
	Point p0 , p1;
	double eps = 1e-9;

	// the edge index which will be intersecting as last
	intersectionEdge_ = -1;


	SUNDANCE_MSG3( verb , "returnIntersectPoints() start= " << start << " end="<< end );
	// make a shorter version of distance calculation
	if ( shortIntersectCalculation(start, end) == true ) { return; }

    for (int ii = 0 ; ii < polyPoints_.size() ; ii++)
    {
    	if (ii < polyPoints_.size() - 1)
    	{
        	p0 = polyPoints_[ii];
        	p1 = polyPoints_[ii+1];
    	}
    	else
    	{   // the last line from size-1 to 0
    		if (!closedPolygon_) continue; // only if the polygon is closed
        	p0 = polyPoints_[ii];
        	p1 = polyPoints_[0];
    	}
    	Point intP(0.0,0.0);
    	// this is the formula where twice the determinant should be zero
    	double a[4] = { p0[1]-p1[1] , p1[0]-p0[0] ,
    			        start[1]-end[1] , end[0]-start[0]};
    	double b[2] = { p0[1]*(p1[0]-p0[0]) - p0[0]*(p1[1]-p0[1]),
    			        start[1]*(end[0] - start[0]) - start[0]*(end[1] - start[1]) };
    	if ( ::fabs( a[0]*a[3]-a[1]*a[2]) < 1e-12){
    		// matrix is singular det|A|~0, no solution
    		continue;
    	}

    	// simple Gauss elimination
       	if (::fabs(a[0]) < eps ){
           	double fakt = -a[0]/a[2];
           	intP[1] = (b[0] + fakt*b[1])/(a[1] + fakt*a[3]);
           	intP[0] = (b[1] - a[3]*intP[1]) / a[2];
        }
        else
        {
           	double fakt = -a[2]/a[0];
           	intP[1] = (b[1] + fakt*b[0])/(a[3] + fakt*a[1]);
           	intP[0] = (b[0] - a[1]*intP[1]) / a[0];
        }
/*
       	bool between1 = (    ((intP[0] >= p0[0]-eps) && (intP[0] <= p1[0]+eps))
   		     || ((intP[0] <= p0[0]+eps) && (intP[0] >= p1[0]-eps)) )
   		&& (    ((intP[1] >= p0[1]-eps) && (intP[1] <= p1[1]+eps))
       		 || ((intP[1] <= p0[1]+eps) && (intP[1] >= p1[1]-eps)) ) ;
       	bool between2 = (    ((intP[0] >= start[0]-eps) && (intP[0] <= end[0]+eps))
      		     || ((intP[0] <= start[0]+eps) && (intP[0] >= end[0]-eps)) )
      		&& (    ((intP[1] >= start[1]-eps) && (intP[1] <= end[1]+eps))
          		 || ((intP[1] <= start[1]+eps) && (intP[1] >= end[1]-eps)) );
       	bool between21 = ((intP[0] >= start[0]-eps) && (intP[0] <= end[0]+eps)) || ((intP[0] <= start[0]+eps) && (intP[0] >= end[0]-eps));
        bool between22 = ((intP[1] >= start[1]-eps) && (intP[1] <= end[1]+eps)) || ((intP[1] <= start[1]+eps) && (intP[1] >= end[1]-eps));*/
        // now we have the intersection point
        //SUNDANCE_MSG3( verb , "returnIntersectPoints() for line ii= " << ii << " , intP= " << intP << " , p0= " << p0 << " , p1= " << p1
        //		<< " I1:" << between1 << " I2:" << between2 );
        //SUNDANCE_MSG3( verb , " I21:" << between21 << " I22:" << between22);

       	if  (      (    ((intP[0] >= p0[0]-eps) && (intP[0] <= p1[0]+eps))
        		     || ((intP[0] <= p0[0]+eps) && (intP[0] >= p1[0]-eps)) )
        		&& (    ((intP[1] >= p0[1]-eps) && (intP[1] <= p1[1]+eps))
            		 || ((intP[1] <= p0[1]+eps) && (intP[1] >= p1[1]-eps)) )
            	&& (    ((intP[0] >= start[0]-eps) && (intP[0] <= end[0]+eps))
           		     || ((intP[0] <= start[0]+eps) && (intP[0] >= end[0]-eps)) )
           		&& (    ((intP[1] >= start[1]-eps) && (intP[1] <= end[1]+eps))
               		 || ((intP[1] <= start[1]+eps) && (intP[1] >= end[1]-eps)) )
    		)
        {
       		    // the the intersection point to the results since it is between the two end points
        		//SUNDANCE_MSG3( verb , "returnIntersectPoints() found intP=" << intP << " is between points");
        		result.resize(nrPoints+1);
        		result[nrPoints] = intP;
               	nrPoints++;
               	// store the line index which was intersected at last
               	intersectionEdge_ = ii;
        }
    }// from the for loop
    //SUNDANCE_MSG3( verb , "returnIntersectPoints() start= " << start << " , end="<< end << " , nrPoints=" << nrPoints);
}


void Polygon2D::returnIntersect(const Point& start, const Point& end,
		int& nrPoints, Array<double>& result) const
{
	Array<Point> t;
	returnIntersectPoints(start, end, nrPoints, t);

	result.resize(nrPoints);

	// Return coordinates instead of t values
	for (int i = 0; i < nrPoints; i++)
	{
		Point tmp( end[0]-t[i][0]/(end[0]-start[0]) , end[1]-t[i][1]/(end[1]-start[1]) );
		result[i] =  sqrt( tmp*tmp );
	}
}


void Polygon2D::getCellsPolygonesPoint( int maxCellLID , Array<Point>& points) const {
	int siz = 0;
	points.resize(siz);

	for(int j = 0 ; j < pointsMaxCellLID_.size() ; j++ ){
		// test if this point is in the maxCell
		if ( maxCellLID == pointsMaxCellLID_[j] ){
			points.resize(siz+1);
			// add the point to the returned array
			points[siz] = polyPoints_[j];
			siz = siz + 1;
		}
	}
}


void Polygon2D::writeToVTK(const std::string& filename) const
{
     std::ofstream myfile;
     myfile.open(filename.c_str());
     myfile << "# vtk DataFile Version 2.0 \n";
     myfile << "Generated by Sundance::Polygon2D Author: Janos Benk \n";
     myfile << "ASCII\n";
     myfile << "\n";
     myfile << "DATASET UNSTRUCTURED_GRID\n";
     myfile << "POINTS "<< polyPoints_.size() << " float \n";
     myfile << "\n";
     for (int ii = 0 ; ii < polyPoints_.size() ; ii++){
    	 myfile << polyPoints_[ii][0] << " " << polyPoints_[ii][1] << " 0\n";
     }
     myfile << "\n \n";
     myfile << "CELLS " << polyPoints_.size()<< " " << polyPoints_.size()*3 << "\n";
     myfile << "\n";
     for (int ii = 0 ; ii < polyPoints_.size()-1 ; ii++){
    	 myfile << "2 "<< ii << " " << ii +1 << "\n";
     }
     myfile << "2 "<< polyPoints_.size()-1 << " 0\n";
     myfile << "\n \n";
     myfile << "CELL_TYPES " << polyPoints_.size() << "\n";
     myfile << "\n";
     for (int ii = 0 ; ii < polyPoints_.size() ; ii++){
    	 myfile << "3\n";
     }
     myfile << "\n";

     if ( nrScalarField_ > 0) {
    	 myfile << "POINT_DATA " << polyPoints_.size() << "\n";
     }

     // first plot everything in a scalar field
     for (int s = 0 ; s < nrScalarField_ ; s++)
     {
    	 // plot each scalar field as a scalar field
    	 myfile << "SCALARS " << polygonSpaceNames_[s] << " float \n";
    	 myfile << "LOOKUP_TABLE default \n";
    	 for (int ii = 0 ; ii < polygonSpaceValues_[s].size() ; ii++){
    		 myfile << polygonSpaceValues_[s][ii] << " \n";
    	 }
         myfile << "\n";
     }

     // since this is in 2D we could pair 2X to plot these 2 pairs as vectors
     // then plot 2 by the the scalar fields as vectors
     if ( nrScalarField_ > 1 ) {
         for (int s = 0 ; s < nrScalarField_ ; s = s + 2)
         {
        	 myfile << "VECTORS " << polygonSpaceNames_[s] << polygonSpaceNames_[s+1] << "_vect float \n";
        	 for (int ii = 0 ; ii < polygonSpaceValues_[s].size() ; ii++){
        		 myfile << polygonSpaceValues_[s][ii] <<  " " <<  polygonSpaceValues_[s+1][ii] << " 0.0\n";
        	 }
             myfile << "\n";
         }
     }

     myfile.close();
}

void Polygon2D::addEvaluationPointValues(const Mesh& mesh ,
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
			      " Polygon2D::addEvaluationPointValues , Polygon must have a valid mesh");

		// look for the impacted points from the polygon
		for (int ii = 0 ; ii < pointsMaxCellLID_.size() ; ii++){
			// test if the polygon point "ii" is inside the "maxCellLID" cell
			if ( maxCellLID == pointsMaxCellLID_[ii] ){
				pointsToLookAt.push_back( ii );
				nrPointsLookAt = nrPointsLookAt + 1;
			}
		}

		// transform all the points into real coordinates
		// then look for the two nearest one, and take the average one
		// we could also take the nearest one, but the average one
		Array<Point> realCoordPoints(quadPts.size());
		Array<Point> cellsPoint(4);
		int tmp;
		cellsPoint[0] = mesh_->nodePosition( mesh_->facetLID( mesh_->spatialDim() ,maxCellLID,0,0,tmp) );
		//cellsPoint[1] = mesh_->nodePosition( mesh_->facetLID( mesh_->spatialDim() ,maxCellLID,0,1,tmp) );
		//cellsPoint[2] = mesh_->nodePosition( mesh_->facetLID( mesh_->spatialDim() ,maxCellLID,0,2,tmp) );
		cellsPoint[3] = mesh_->nodePosition( mesh_->facetLID( mesh_->spatialDim() ,maxCellLID,0,3,tmp) );
		for (int p = 0 ; p < nrQuadPoint ; p++){
			// calculate the real coordinated
			realCoordPoints[p] = cellsPoint[0] + Point( (cellsPoint[3][0] - cellsPoint[0][0])*quadPts[p][0]  ,
					(cellsPoint[3][1] - cellsPoint[0][1])*quadPts[p][1] );
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
				dist_tmp = (realCoordPoints[p] - polyPoints_[pointsToLookAt[polyP]])*(realCoordPoints[p] - polyPoints_[pointsToLookAt[polyP]]);
				if ( dist_tmp < dist ) {
					dist = dist_tmp;
					firstI = p;
				}
			}
			// look for the second nearest point
			dist = 1e+100;
			for( int p = 0 ; p < nrQuadPoint ; p++){
				dist_tmp = (realCoordPoints[p] - polyPoints_[pointsToLookAt[polyP]])*(realCoordPoints[p] - polyPoints_[pointsToLookAt[polyP]]);
				if ( (dist_tmp < dist) && ( p != firstI ) ) {
					dist = dist_tmp;
					secondI = p;
				}
			}
			//having the two nearest point we add the average value
			SUNDANCE_MSG3( verb , "pointsToLookAt["<<polyP<<"]=" << pointsToLookAt[polyP] << " , v1="
					<< coeffPtr[firstI] << " , v2="<< coeffPtr[secondI]);
			SUNDANCE_MSG3( verb , " Set Point =" << polyPoints_[pointsToLookAt[polyP]] << " , p0="
					<< realCoordPoints[firstI] << " , p1="<< realCoordPoints[secondI]);
			// here we just set the value since one point should be set only once
			polygonSpaceValues_[recordedScalarField_][pointsToLookAt[polyP]] = 0.5*(coeffPtr[firstI] + coeffPtr[secondI]);
		}
	}
}


void Polygon2D::setSpaceValues(const FunctionalEvaluatorBase& scalarFunctional , int fieldIndex ){

	int verb = 0;
	SUNDANCE_MSG3( verb , "setSpaceValues() fieldIndex= " << fieldIndex );

	// signaling that we are starting recording values
	isRecording_ = true;
	recordedScalarField_ = fieldIndex;
	SUNDANCE_MSG3( verb , " set values to zero " );
	// set all values to zero because now the recorded values will be added to the vector elements
	for (int ii = 0 ; ii < polygonSpaceValues_[fieldIndex].size() ; ii++ ){
		polygonSpaceValues_[fieldIndex][ii] = 0.0;
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
			  Array<double> reduceVector( polygonSpaceValues_[fieldIndex].size() , 0.0 );
			  // call the MPI function to reduce the vector
			  SUNDANCE_MSG3( verb , " call the reduce function " );
			  mesh_->comm().allReduce((void*) &(polygonSpaceValues_[fieldIndex][0]), (void*) &(reduceVector[0]),
			    polygonSpaceValues_[fieldIndex].size() ,MPIDataType::doubleType(), MPIOp::sumOp());
			  // copy the reduced vector back
			  SUNDANCE_MSG3( verb , " copy the reduced values back " );
			  for (int ii = 0 ; ii < polygonSpaceValues_[fieldIndex].size() ; ii++ ){
				  polygonSpaceValues_[fieldIndex][ii] = reduceVector[ii];
			  }
		  }
	}

	// if the twin polygon exists then copy also there the value
	if (twinPolygon_ != 0){
		SUNDANCE_MSG3( verb , " update twin polygon " );
		// copy the recorded values to the twin polygon
		for (int ii = 0 ; ii < polygonSpaceValues_[fieldIndex].size() ; ii++ ){
			twinPolygon_->polygonSpaceValues_[fieldIndex][ii] = this->polygonSpaceValues_[fieldIndex][ii];
		}
	}

	// we set the flags so that we stop recording
	recordedScalarField_ = -1;
	isRecording_ = false;
}

void Polygon2D::eval0(const double* vars, double* f , int scalarFieldIndex ) const {
	  double x = vars[0];
	  double y = vars[1];
	  Point inCoord(x,y);
	  double eps = 1e-14;
	  int verb = 0;

	  SUNDANCE_MSG3( verb , "eval0() x= " << x << " y=" << y);
	  // just look for the two nearest point, then look which neighbor point is the nearest
	  // calculate the projection of that point to the line
	  int nrPolygonPoints = polyPoints_.size() , pointIndex = -1;
	  double dist = 1e+100 , dist_tmp;
	  for (int ii = 0 ; ii < nrPolygonPoints ; ii++){
		  dist_tmp = (inCoord - polyPoints_[ii])*(inCoord - polyPoints_[ii]);
		  //SUNDANCE_MSG4( verb ,"polyPoints_["<<ii<<"]="<<polyPoints_[ii]<< "  ,  dist["<<ii<<"]= " << dist_tmp);
		  if ( dist_tmp < dist ) {
			  dist = dist_tmp;
			  pointIndex = ii;
		  }
	  }

	  // get the two neighbor point candidate
	  int pointL = -1 , pointR = -1;
	  if ( pointIndex == nrPolygonPoints-1 ) {
		  if ( closedPolygon_ ){
			  pointL = pointIndex-1 , pointR = 0;
		  } else {
			  pointL = pointIndex-1 , pointR = pointIndex-1;
		  }
	  }
	  else {
		  if ( pointIndex == 0 ) {
			  if ( closedPolygon_ ){
				  pointL = nrPolygonPoints-1 , pointR = pointIndex+1;
			  } else {
				  pointL = pointIndex+1 , pointR = pointIndex+1;
			  }
		  }
		  else {
			  pointL = pointIndex-1 , pointR = pointIndex+1;
		  }
	  }

	  double distL = (inCoord - polyPoints_[pointL])*(inCoord - polyPoints_[pointL]);
	  double distR = (inCoord - polyPoints_[pointR])*(inCoord - polyPoints_[pointR]);

	  int pointIndexN = ( distL > distR) ? pointR : pointL;

	  // now calculate the projection of the point inCoord on the line [pointIndex , pointIndexN]
	  // and then limit the contribution of the two points between [0,1] (make a saturation to 0 or 1)

	  Point p0 = polyPoints_[pointIndex];
	  Point p1 = polyPoints_[pointIndexN];
	  SUNDANCE_MSG3( verb , "pI0 = " << pointIndex << " pI1=" << pointIndexN );
	  SUNDANCE_MSG3( verb , "p0 = " << p0 << " p1=" << p1);
	  double oneDimScale = 0.0;
	  Point intP(0.0,0.0);
	  // to calculate the projection of a point on a line we set up a system with two equations (for more detail see the PDF)
	  // we solve this system and we get the intersection point (if the system is not singular)
	  double a[4] = { p1[0]-p0[0] , p1[1]-p0[1] ,
				      p0[1]-p1[1] , p1[0]-p0[0] };
	  double b[2] = { - inCoord[0]*(p1[0]-p0[0]) - inCoord[1]*(p1[1]-p0[1]),
					  - p0[1]*(p1[0]-p0[0]) + p0[0]*(p1[1]-p0[1]) };
	  if ( ::fabs( a[0]*a[3]-a[1]*a[2]) < eps){
		  // matrix is singular det|A|~0, no solution, so the point must be on the line
		  if ( ::fabs(p1[0] - p0[0]) > eps )
			  { oneDimScale = (p1[0] - inCoord[0])/(p1[0] - p0[0]); }
		  else
			  { oneDimScale = (p1[1] - inCoord[1])/(p1[1] - p0[1]); }
	  }
	  else
	  {   // simple Gauss elimination
		  if (::fabs(a[0]) < eps ){
         	double fakt = -a[0]/a[2];
         	intP[1] = (b[0] + fakt*b[1])/(a[1] + fakt*a[3]);
         	intP[0] = (b[1] - a[3]*intP[1]) / a[2];
         	intP = (-1)*intP;
		  }
		  else
		  {
         	double fakt = -a[2]/a[0];
         	intP[1] = (b[1] + fakt*b[0])/(a[3] + fakt*a[1]);
         	intP[0] = (b[0] - a[1]*intP[1]) / a[0];
         	intP = (-1)*intP;
		  }
		  if ( ::fabs(p1[0] - p0[0]) > eps )
			  { oneDimScale = (p1[0] - intP[0])/(p1[0] - p0[0]); }
		  else
			  { oneDimScale = (p1[1] - intP[1])/(p1[1] - p0[1]); }
	  }
	  // limit this factor to [0,1] and then calculate the value and return it back in f[0]
	  SUNDANCE_MSG3( verb , "intP = " << intP << " oneDimScale=" << oneDimScale);
	  oneDimScale = ( oneDimScale > 1.0 ) ? 1.0 : oneDimScale;
	  oneDimScale = ( oneDimScale < 0.0 ) ? 0.0 : oneDimScale;
	  f[0] = oneDimScale*polygonSpaceValues_[scalarFieldIndex][pointIndex] +
			  (1.0-oneDimScale)*polygonSpaceValues_[scalarFieldIndex][pointIndexN];
	  SUNDANCE_MSG3( verb , "eval0() ret = " << f[0] << " oneDimScale=" << oneDimScale);
}

Polygon2D* Polygon2D::createTwinPolygon(double shiftX , double shiftY , double scaleX , double scaleY) {

	Array<Point> points(polyPoints_.size());

	// make the shift and the scaling
	for (int j = 0 ; j < polyPoints_.size() ; j++ ){
		points[j] = Point( shiftX + scaleX*polyPoints_[j][0] , shiftY + scaleY*polyPoints_[j][1] );
	}

	// create the polygon
	Polygon2D* polyg = new Polygon2D( points , this->_alpha1 , this->_alpha2 , closedPolygon_ , flipDomains_ );

	// add all the scalar fields in the same order
	for (int j = 0 ; j < nrScalarField_ ; j++){
		polyg->addNewScalarField( polygonSpaceNames_[j], 0.0 );
	}

	// store the twin polygon mutually
	twinPolygon_ = polyg;
	polyg->twinPolygon_ = this;

	return polyg;
}

RCP<CurveBase> Polygon2D::unite(ParametrizedCurve& c1 , ParametrizedCurve& c2)
{
   const CurveBase* pb1 = c1.ptr().get();
   const CurveBase* pb2 = c2.ptr().get();

   const Polygon2D* pol1 = dynamic_cast<const Polygon2D*>(pb1);
   const Polygon2D* pol2 = dynamic_cast<const Polygon2D*>(pb2);

   if ( (pol1 == 0) || (pol2 == 0)){
	   TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error, "Polygon2D::unite one of the inputs is not a polygon 2D ");
   }

   // the array where the resulting points will be stored
   Array<Point> allPoints(0);

   int nrP1 = pol1->polyPoints_.size();
   int nrP2 = pol2->polyPoints_.size();
   int p1Index = 0;
   int p2Index = 0;
   int pNewIndex = 0;
   int polyIndex = 0;
   bool doLoop = true;
   int verb = 0;

   SUNDANCE_MSG3( verb , "unite() starting ... " );
   // 1) start with the frst polygon, iterate as long one point will be outside
   while (doLoop)
   {
	   SUNDANCE_MSG3( verb , "unite() testing p1Index=" << p1Index );
	  if ((pol2->curveEquation(pol1->polyPoints_[p1Index+1]) > 0) || (p1Index >= nrP1 - 1))
	  {
		  doLoop = false;
		  // when already the first point is inside then add to the new polygon the first point
		  if (pol2->curveEquation(pol1->polyPoints_[p1Index]) > 0)
		  {
			   allPoints.resize(pNewIndex+1);
			   allPoints[ pNewIndex ] = pol1->polyPoints_[p1Index];
			   pNewIndex = pNewIndex + 1;
		  }
	  }
	  else
	  {
		  p1Index = p1Index + 1;
	  }
   }
   SUNDANCE_MSG3( verb , "unite() for starting found p1Index=" << p1Index << " pNewIndex=" << pNewIndex );

   // 2) iterate as long it is outside the other, with the current polygon
   // 3) if we have switch , then determine intersection point, and continue with the other polygon
   // 4) stop if we run out of points in the case of one polygon
   while ( (p1Index < nrP1) && (p2Index < nrP2))
   {
	   SUNDANCE_MSG3( verb , "unite() next point polyIndex = " << polyIndex);
	   if (polyIndex == 0)
	   { // first polygon
		   if ( (pol2->curveEquation(pol1->polyPoints_[p1Index]) > 0) &&
				(pol2->curveEquation(pol1->polyPoints_[p1Index+1]) > 0) ) {
			   // the next point is also inside and
			   allPoints.resize(pNewIndex+1);
			   allPoints[ pNewIndex ] = pol1->polyPoints_[p1Index];
			   SUNDANCE_MSG3( verb , "unite() step forward with polygon1 " << p1Index << " pNewIndex="
					    << pNewIndex << " allPoints[ pNewIndex ]=" << allPoints[ pNewIndex ]);
			   pNewIndex = pNewIndex + 1;
			   p1Index = p1Index + 1;
		   } else {
			   // calculate intersection point
			   Array<Point> intesectPoints;
			   int nrIntersectPoints;
			   pol2->returnIntersectPoints( pol1->polyPoints_[p1Index] , pol1->polyPoints_[p1Index+1] , nrIntersectPoints , intesectPoints );
			   SUNDANCE_MSG3( verb , "unite() calculate intersection of p2 with p1Index = " << p1Index << " nrP=" << nrIntersectPoints
					   << " , " << pol1->polyPoints_[p1Index] << " , " << pol1->polyPoints_[p1Index+1]
					  << " , " << intersectionEdge_ << " p2I=" << p2Index);
			   if (nrIntersectPoints > 0)
			   {
				   if (intersectionEdge_ < p2Index) {
					   p2Index = nrP2; // exit the whole while loop
					   break;
				   }
				   // add the first point and switch if we are not at the beginning
				   if (pNewIndex > 1){
					   polyIndex = 1;
					   allPoints.resize(pNewIndex+1);
					   allPoints[ pNewIndex ] = pol1->polyPoints_[p1Index];
					   pNewIndex = pNewIndex + 1;
					   p2Index = intersectionEdge_ + 1;
				   }
				   else{
					   p1Index = p1Index + 1;
				   }
				   allPoints.resize(pNewIndex+1);
				   allPoints[ pNewIndex ] = intesectPoints[nrIntersectPoints-1];
				   SUNDANCE_MSG3( verb , "unite() switch from p1 to p2 p2Index=" << p2Index << " pNewIndex=" <<
						   pNewIndex << " allPoints[ pNewIndex ]=" << allPoints[ pNewIndex ] );
				   pNewIndex = pNewIndex + 1;
				   // switch the polygon , only if this is the first point
			   }
			   else
			   {   // error
				   TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error , "Polygon2D::unite nrIntersectPoints == 0");
			   }
		   }
	   }
	   else
	   { // second polygon outside
		   if ( (pol1->curveEquation(pol2->polyPoints_[p2Index]) > 0) &&
				(pol1->curveEquation(pol2->polyPoints_[p2Index+1]) > 0) ) {
			   // the next point is also inside and
			   allPoints.resize(pNewIndex+1);
			   allPoints[ pNewIndex ] = pol2->polyPoints_[p2Index];
			   SUNDANCE_MSG3( verb , "unite() step forward with polygon2 " << p2Index << " pNewIndex="
					    << pNewIndex << " allPoints[ pNewIndex ]=" << allPoints[ pNewIndex ]);
			   pNewIndex = pNewIndex + 1;
			   p2Index = p2Index + 1;
		   } else {
			   // calculate intersection point
			   Array<Point> intesectPoints;
			   int nrIntersectPoints;
			   pol2->returnIntersectPoints( pol2->polyPoints_[p2Index] , pol2->polyPoints_[p2Index+1] , nrIntersectPoints , intesectPoints );
			   SUNDANCE_MSG3( verb , "unite() calculate intersection of p1 with p2Index = " << p2Index << " nrP=" << nrIntersectPoints
					   << " , " << pol2->polyPoints_[p2Index] << " , " << pol2->polyPoints_[p2Index+1]
					   << " , " << intersectionEdge_ << " p1I=" << p1Index);
			   // todo: for circle some reason no intersection point is found
			   if ( (nrIntersectPoints > 0) || (intersectionEdge_ < p1Index) )
			   {
				   if (intersectionEdge_ < p1Index) {
					   p1Index = nrP1; // exit the whole while loop
					   break;
				   }
				   p1Index = intersectionEdge_ + 1;
				   allPoints.resize(pNewIndex+1);
				   allPoints[ pNewIndex ] = pol2->polyPoints_[p2Index];
				   pNewIndex = pNewIndex + 1;
				   allPoints.resize(pNewIndex+1);
				   allPoints[ pNewIndex ] = intesectPoints[nrIntersectPoints-1];
				   SUNDANCE_MSG3( verb , "unite() switch from p2 to p1 p1Index=" << p1Index << " pNewIndex=" <<
						   pNewIndex << " allPoints[ pNewIndex ]=" << allPoints[ pNewIndex ] );
				   pNewIndex = pNewIndex + 1;
				   // switch the polygon
				   polyIndex = 0;
			   }
			   else
			   {   // error
				   TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error , "Polygon2D::unite nrINtersectPoints == 0");
			   }
		   }
	   }
   }

   // now create one polygon
   if (pol1->hasMesh_)
	   return rcp( new Polygon2D(*(pol1->mesh_) , allPoints , pol1->_alpha1 , pol2->_alpha2 ) );
   else
	   return rcp( new Polygon2D( allPoints , pol1->_alpha1 , pol2->_alpha2 ) );
}
