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

#ifndef SUNDANCEPOLYGON2D_H_
#define SUNDANCEPOLYGON2D_H_

#include "SundanceCurveBase.hpp"
#include "SundanceParametrizedCurve.hpp"

namespace Sundance
{

 //class FunctionalEvaluator;

/**  */
class Polygon2D: public Sundance::CurveBase
{
public:

	/** Ctor with the points which form the polygon (the direction of the points matters for IN & OUT )*/
	Polygon2D(const Mesh& mesh , const Array<Point>& points , double a1 , double a2, bool closedPolygon = true ,bool flipD = false);

	/** Ctor with the points which form the polygon (the direction of the points matters for IN & OUT )*/
	Polygon2D(const Array<Point>& points , double a1 , double a2, bool closedPolygon = true , bool flipD = false);

	/** Ctor to read in the polygon from file from file */
	Polygon2D(const Mesh& mesh , const std::string& filename , double a1 , double a2, bool closedPolygon = true ,bool flipD = false);

	/** Ctor to read in the polygon from file from file */
	Polygon2D(const std::string& filename , double a1 , double a2, bool closedPolygon = true ,bool flipD = false);

	/** */
	virtual ~Polygon2D() {;}

	/** @return Expr The parameters of the curve which uniquely defines the curve*/
	virtual Expr getParams() const;

	/** return the points of the polygon */
	virtual Array<Point>& getControlPoints() { return polyPoints_; }

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

	/** returns the points which are inside one maxCell */
	void getCellsPolygonesPoint( int maxCellLID , Array<Point>& points) const ;

	/** see super class */
	virtual int addNewScalarField(std::string fieldName , double initialValue){
		polygonSpaceValues_.resize(nrScalarField_+1);
		polygonSpaceNames_.resize(nrScalarField_+1);
		polygonSpaceValues_[nrScalarField_].resize(polyPoints_.size() , initialValue);
		polygonSpaceNames_[nrScalarField_] = fieldName;
		nrScalarField_ = nrScalarField_ + 1;
		return nrScalarField_ - 1;
	}

	/** returns one array of doubles which contain the values of the scalar field (for each point, nodal basis)
	 * @param scalarFieldIndex [IN] the index
	 * @return array with the scalar field values */
	virtual Array<double>& getScalarFieldValues(int scalarFieldIndex) { return polygonSpaceValues_[scalarFieldIndex]; }

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

	/** create the twin polygon and then also scale and shift the polygon
	 * @param shiftX
	 * @param shiftY
	 * @param scaleX
	 * @param scaleY */
	Polygon2D* createTwinPolygon(double shiftX , double shiftY , double scaleX , double scaleY);

	/** generate the polygon which is the unification (not intersection) of two polygons */
    static RCP<CurveBase> unite(ParametrizedCurve& c1 , ParametrizedCurve& c2);

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

	/** indicates weather we have a mesh or not */
	bool hasMesh_;

	/** this flag shows if the last point is connected to the first one*/
	bool closedPolygon_;

	/** the mesh will be needed for point location, since the geometry is
	 * non-conform with respect to the mesh */
	const Mesh* mesh_;

	/** The points in the specified order which form */
	Array<Point> polyPoints_;

	/** each points should know in which cell is evaluable (DiscreteSpace)*/
	Array<int> pointsMaxCellLID_;

	/** coordinates to store the bounding box*/
	double minX_;
	double maxX_;
	double minY_;
	double maxY_;

// ----------- values on the polygon ----------------
	/** values created on the polygon */
	mutable Array< Array<double> > polygonSpaceValues_;

	/** values created on the polygon */
	mutable Array< std::string > polygonSpaceNames_;

	/** number of scalar fields*/
	int nrScalarField_;

	/** variable needed for setting the values on the polygon*/
	mutable bool isRecording_;

	/** store which field are we recording */
	mutable int recordedScalarField_;

// -------------- twin polygon for Lagrange - Euler framework interface ----------
	/** pointer to the polygon to which this is coupled (like a twin)*/
	mutable Polygon2D* twinPolygon_;

// ------------- tmp variable for polygon union ----------
	/** the index of the line which is at last the intersection point,
	 * used only in the union of two polygons */
	static int intersectionEdge_;

};

}
#endif /* SUNDANCEPOLYGON2D_H_ */
