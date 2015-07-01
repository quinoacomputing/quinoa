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
 * SundanceTriangleSurf3D.hpp
 *
 *  Created on: Oct 13, 2011
 *      Author: benk
 */


#ifndef SUNDANCETRIANGLESURF3D_HPP_
#define SUNDANCETRIANGLESURF3D_HPP_

#include "SundanceCurveBase.hpp"
#include "SundanceParametrizedCurve.hpp"

namespace Sundance
{

/**  */
class TriangleSurf3D: public Sundance::CurveBase
{
public:

	/** Ctor with the points which form the polygon (the direction of the points matters for IN & OUT )*/
	TriangleSurf3D(const Mesh& mesh , const Array<Point>& points ,
			const Array<int>& triags , double a1 , double a2, bool flipD = false);

	/** Ctor with the points which form the polygon (the direction of the points matters for IN & OUT )*/
	TriangleSurf3D(const Array<Point>& points ,
			       const Array<int>& triags , double a1 , double a2, bool flipD = false);

	/** Ctor to read in the polygon from file from file */
	TriangleSurf3D(const Mesh& mesh , const std::string& filename , double a1 , double a2, bool flipD = false);

	/** Ctor to read in the polygon from file from file */
	TriangleSurf3D(const std::string& filename , double a1 , double a2, bool flipD = false);

	/** */
	virtual ~TriangleSurf3D() {;}

	/** @return Expr The parameters of the curve which uniquely defines the curve*/
	virtual Expr getParams() const;

	/** return the points of the polygon */
	virtual Array<Point>& getControlPoints() { return triagPoints_; }

	/** update the state of the curve of the control point were changed */
	virtual void update();

	/**
	 * This function is important for nonstructural mesh integration.<br>
	 * The function returns the intersection points with a given line (only 2D at the moment)<br>
	 * The line is defined with a starting point and an end point<br>
	 * The resulting points should be between the two points
	 * @param start, the start point of the line
	 * @param end, the end point of the line
	 * @param nrPoints , number of resulted (intersected) points
	 * @param result , the resulted points of intersections */
	virtual void returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
			Array<Point>& result) const;

	/**
	 * As above, but, instead of coordinates, returns intersection values t in [0,1]
	 * and the line is defined by "start + t * (end-start)"
	 */
	virtual void returnIntersect(const Point& start, const Point& end, int& nrPoints,
			Array<double>& result) const;

	/** Return a ref count pointer to self */
	virtual RCP<CurveBase> getRcp()
	{
		return rcp(this);
	}

	/** This curve is a real curve */
	virtual bool isCurveValid() const
	{
		return true;
	}

	/** sets the mesh for the polygon , empty implementation*/
	virtual void setMesh(const Mesh& mesh);

	/** Writes the geometry into a VTK file for visualization purposes
	 * @param filename */
	virtual void writeToVTK(const std::string& filename) const;

	/** see super class */
	virtual int addNewScalarField(std::string fieldName , double initialValue){
		triagSurfSpaceValues_.resize(nrScalarField_+1);
		triagSurfSpaceNames_.resize(nrScalarField_+1);
		triagSurfSpaceValues_[nrScalarField_].resize(triagPoints_.size() , initialValue);
		triagSurfSpaceNames_[nrScalarField_] = fieldName;
		nrScalarField_ = nrScalarField_ + 1;
		return nrScalarField_ - 1;
	}

	/** returns one array of doubles which contain the values of the scalar field (for each point, nodal basis)
	 * @param scalarFieldIndex [IN] the index
	 * @return array with the scalar field values */
	virtual Array<double>& getScalarFieldValues(int scalarFieldIndex) { return triagSurfSpaceValues_[scalarFieldIndex]; }

	/** function to set the values of one specified scalar field <br>
	 * it is important that the Functional will be defined on a curve, otherwise the value will stay zero
	 * @param scalarFunctional [IN] , the scalar functions*/
	virtual void setSpaceValues(const FunctionalEvaluatorBase& scalarFunctional , int fieldIndex );

	/** function to record the function values along one interface. <br>
	 * This function will be called from the integration routine.*/
	virtual void addEvaluationPointValues(const Mesh& mesh ,
			int maxCellLID , int nQuad ,
			const double* coeffPtr ,
			const Array<Point>& quadPts) const ;

	/** function which will be called for curve Integral evaluation
	 * @param vars [IN] coordinates
	 * @param f [OUT] the output of the expression
	 * @param scalarFieldIndex [IN] scalar field index */
	virtual void eval0(const double* vars, double* f , int scalarFieldIndex ) const;

	/** create the twin polygon and then also scale and shift the polygon */
	TriangleSurf3D* createTwinTriangleSurf(double shiftX , double shiftY , double shiftZ ,
			 double scaleX , double scaleY , double scaleZ );

	/** returns the pointer to the twin polygon */
	TriangleSurf3D* getTwinTriangleSurf() const { return twinTriangleSurf_; }

	/** sets the pointer to the twin polygon */
	void setTwinTriangleSurf(TriangleSurf3D* surf) const { twinTriangleSurf_ = surf; }

	/** reads in a geometry from a GTS file format */
    static TriangleSurf3D* importGTSSurface( const std::string& filename , double a1 , double a2 , Point innerP);

	/** reads in a geometry from a STL file format */
    static TriangleSurf3D* importSTLSurface( const std::string& filename , double a1 , double a2 , Point innerP);

protected:
	/**
	 * This function should be implemented
	 * @param evalPoint point where the polygon equation is evaluated <br>
	 * @return double the value of the curve equation at the evaluation point  */
	virtual double curveEquation_intern(const Point& evalPoint) const;
private:

	/** looks for each point in which cell it is contained (later for evaluation purposes)*/
	void computeMaxCellLIDs();

	/** */
	bool shortDistanceCalculation(const Point& p, double &res) const;

	/** */
	bool shortIntersectCalculation(const Point& st, const Point& end) const;

	/** The function returns the */
	const Array<int>& getSpaceTree(const Point& p , int& spaceTreeCellID) const
	{
		// todo: later here should be a real space tree
		spaceTreeCellID = 0;
		return triagIDs_;
	}

	/** The function returns the */
	void setupSpaceTree() { /* TODO: implement this*/ }

	/** returns true if a point is inside the given triangle*/
	inline bool pointInTriangle(const Point& A, const Point& B, const Point& C , const Point& P) const {
		double eps = 1e-10;
	    //Compute vectors
		Point v0 = C - A;
		Point v1 = B - A;
		Point v2 = P - A;
	    //Compute dot products
	    double dot00 = v0 * v0;
	    double dot01 = v0 * v1;
	    double dot02 = v0 * v2;
	    double dot11 = v1 * v1;
	    double dot12 = v1 * v2;
	    // Compute barycentric coordinates
	    double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	    if ((u >= 0-eps ) && (v >= 0-eps) && (u + v <= 1+eps ))
	    	return true;
	    else
	    	return false;
	}

	/** returns the one intersection points by one triangle and one line */
	inline Point triangleLineIntersect(const Point& A, const Point& B, const Point& C ,
			             const Point& L0 , const Point& L1 , bool& isIntersected , double& t) const {
		/*	x = zeros(3,1); sI = -1.0; u = P1 - P0; w = P0 - A;
			n = cross(B-A,C-A); D = dot(n,u); N = -dot(n,w);
			if (abs(D)< 1e-10) % no solution
			else sI = N/D; x = P0 + sI*u; end; */
			 Point x(0.0,0.0,0.0);
			 t = -1.0; isIntersected = true;
			 Point u = L1-L0;
			 Point w = L0-A;
			 Point n = cross(C-A,B-A);
			 double D = n*u;
		     double N = -n*w;
		 	if (::fabs(D) < 1e-10){
		 		isIntersected = false;
		    } else {
		 		t = N/D;
		 		x = L0 + t*u;
		 	}
		 	return x;
	}

	/** Projects one point to one line */
	inline void projectPointToLine(const Point& A, const Point& B, const Point& P , Point& Pr) const {
		/*	D = dot((b-a),(b-a)); N = - dot( (b-a) , (a-p)); ProjectPoint = zeros(3,1);
			if (abs(D)>1e-10) t = N/D; t = min( max(0,t) , 1.0); ProjectPoint = a + t*(b-a); end */
			double D = (B-A)*(B-A);
			double N = -(B-A)*(A-P);
			if (::fabs(D) > 1e-10){
					double t = N/D;
					t = ::fmin(1,::fmax(t,0));
					Pr = A + t*(B-A);
			}
	}

	/** indicates weather we have a mesh or not */
	bool hasMesh_;

	/** the mesh will be needed for point location, since the geometry is
	 * non-conform with respect to the mesh */
	const Mesh* mesh_;

	/** The points in the specified order which form */
	Array<Point> triagPoints_;

	/** for each point we store in which cell it is located so by setting the values at surface points
	 *  it will be done more quickly */
	Array<int> pointMaxCellLID_;

	/** The index of the points in one triangle */
	Array<int> triagIndexes_;

	/** The index of the points in one triangle */
	Array<int> triagIDs_;

	/** coordinates to store the bounding box*/
	double minX_;
	double maxX_;
	double minY_;
	double maxY_;
	double minZ_;
	double maxZ_;

	/** values created on the polygon */
	mutable Array< Array<double> > triagSurfSpaceValues_;

	/** values created on the polygon */
	mutable Array< std::string > triagSurfSpaceNames_;

	/** number of scalar fields*/
	int nrScalarField_;

	/** variable needed for setting the values on the polygon*/
	mutable bool isRecording_;

	/** store which field are we recording */
	mutable int recordedScalarField_;

// -------------- twin polygon for Lagrange - Euler framework interface ----------
	/** pointer to the polygon to which this is coupled (like a twin)*/
	mutable TriangleSurf3D* twinTriangleSurf_;

};

}
#endif /* SUNDANCETRIANGLESURF3D_HPP_ */
