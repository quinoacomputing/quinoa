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

#ifndef SUNDANCEPARAMETRIZEDCURVE_H_
#define SUNDANCEPARAMETRIZEDCURVE_H_

#include "PlayaHandle.hpp"
#include "SundanceDummyParametrizedCurve.hpp"
#include "SundanceCurveBase.hpp"

/** This specifies the different filter modes
 *  If we have a closed curve (see Circle) than we inside that curve the "curveEquation"
 *  function should be negative <br>
 *  On the outside the curveEquation function should return a positive value. */
enum CurveCellFilterMode {Outside_Curve = 1 ,
	                      Inside_Curve = -1,
	                      On_Curve = 0 };

namespace Sundance
{

/**
 * This class a "handle" to the parameterized curves, which might be used
 * in Adaptive Cell Integration.<br> e.g : the curve specifies
 * which is my computational domain and which is not.<br>
 * But this might be used for other purposes as well.
 *
 */
class ParametrizedCurve : public Playa::Handle<CurveBase> {
public:

	/* */
	/*HANDLE_CTORS(ParametrizedCurve, CurveBase);*/

    /** Construct an empty ParametrizedCurve type object */
	ParametrizedCurve();

    /** Construct from a raw pointer to a ParametrizedCurve type subtype */
	ParametrizedCurve(Playa::Handleable<CurveBase>* rawPtr);

    /** Construct from a smart pointer to a ParametrizedCurve type subtype */
	ParametrizedCurve(const RCP<CurveBase>& smartPtr);

	virtual ~ParametrizedCurve();

	/** Returns the parameters of the curve*/
	Expr getParams() const {
		return ptr()->getParams();
	}

	/** return the control points of the parameterized curve */
	Array<Point>& getControlPoints() {
		return ptr()->getControlPoints();
	}

	/** update the state of the curve of the control point were changed */
	void update() {
		ptr()->update();
	}

	/** function which will be called for curve/surface Integral evaluation */
	void eval0(const double* vars, double* f , int scalarFieldIndex ) const {
		ptr()->eval0(vars, f , scalarFieldIndex );
	}

	/** adds new scalar field to the interface*/
	int addNewScalarField(std::string fieldName , double initialValue){
		return ptr()->addNewScalarField(fieldName , initialValue);
	}

	/** returns one array of doubles which contain the values of the scalar field (for each point, nodal basis)*/
	Array<double>& getScalarFieldValues(int scalarFieldIndex) {
		return ptr()->getScalarFieldValues(scalarFieldIndex);
	}

	/** function to set the values of one specified scalar field <br>
	 * it is important that the Functional will be defined on a curve, otherwise the value will stay zero */
	void setSpaceValues(const FunctionalEvaluatorBase& scalarFunctional , int fieldIndex ){
		ptr()->setSpaceValues(scalarFunctional , fieldIndex );
	}

	/*Function to record the function values along one curve. <br>
	 *This function will be called from the integration routine.*/
	void addEvaluationPointValues(const Mesh& mesh ,
			int maxCellLID , int nQuad ,
			const double* coeffPtr ,
			const Array<Point>& quadPts) const {
		ptr()->addEvaluationPointValues( mesh , maxCellLID , nQuad , coeffPtr , quadPts);
	}

	/** List integration parameters for the FCM method*/
	void getIntegrationParams(Array<double>& alphas) const {
		return ptr()->getIntegrationParams(alphas);
	}

	/** return the integration coefficient for the inner */
	double getAlpha1() const { return ptr()->getAlpha1(); }

	/** return the integration coefficient for the outer */
	double getAlpha2() const { return ptr()->getAlpha2(); }

	/** return the dimension of the curve , in which it is defined */
	int getCurveDim() const { return ptr()->getCurveDim(); }

	/** Returns the integration parameter for the FCM method*/
	double integrationParameter(const Point& evaluationPoint) const {
		return ptr()->integrationParameter(evaluationPoint);
	}

	/** The curve(2D, 3D curve, 3D surface) equation */
	double curveEquation(const Point& evaluationPoint) const {
		return ptr()->curveEquation(evaluationPoint);
	}

	/** flip the domains from the ficticious to the real */
	void flipDomains() const {
		ptr()->flipDomains();
	}

	/** Returns the points of intersection of the curve with the line defined by the points <br>
	 *  Only intersection between start and end point will be stated */
	void returnIntersectPoints(const Point& startEdgePoint, const Point& endEdgePoint,
	                      int& nrPoints ,Array<Point>& result) const {
		ptr()->returnIntersectPoints( startEdgePoint , endEdgePoint , nrPoints , result);
	}

	/**
	 *  Same as above, but intersection parameters t will be returned where
	 *  the line is defined by start+t*(end-start). Consequently all t's will be
	 *  in the range [0,1]
	 */
	void returnIntersect(const Point& startEdgePoint, const Point& endEdgePoint,
		                      int& nrPoints ,Array<double>& result) const {
			ptr()->returnIntersect( startEdgePoint , endEdgePoint , nrPoints , result);
		}

	/** In the case of simple geometries the geometry it can be transformed to a polygon, which
	 * reflects the original geometry
	 * @param mesh
	 * @param resolution , the global resolution */
	const RCP<CurveBase> getPolygon(const Mesh& mesh , double resolution) const {
		return ptr()->getPolygon( mesh , resolution);
	}

	/** Some parameterized curve types need also the mesh
	 * @param mesh */
	void setMesh(const Mesh& mesh) {
		ptr()->setMesh( mesh );
	}

	/** Writes the geometry into a VTK file for visualization purposes
	 * @param filename */
	void writeToVTK(const std::string& filename) const {
		ptr()->writeToVTK(filename);
	}

	/** Shows if the curve is a valid curve*/
	inline bool isCurveValid() const { return ptr()->isCurveValid(); }

	/** function which shows if some integral schould be calculated along the curve, <br>
	 * True means we have to calculate a curve/surf integral */
	inline bool isCurveIntegral() const { return ptr()->isCurveIntegral(); }

	/** Static class returns the dummy class when this is needed */
	static const ParametrizedCurve& returnDummyCurve() {return handle_; }

	/** The operator to determine the uniqueness of one curve
	 *  compare to another ParametrizedCurve, used for placement in STL containers*/
	bool operator<(const ParametrizedCurve& other) const { return (myID_ < other.myID_)?true:false;}

	/** return the ID of the curve */
	int myID() const { return myID_; }

private:

	/** The ID to identify the ParamCurve uniquely*/
	int myID_;

	/* The dummy curve which does not do anything */
	static const ParametrizedCurve handle_;

	/** used to give global IDs*/
	static int globalId_;
};

}
#endif /* SUNDANCEPARAMETRIZEDCURVE_H_ */
