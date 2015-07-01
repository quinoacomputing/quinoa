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

#ifndef SUNDANCECURVEBASE_H_
#define SUNDANCECURVEBASE_H_

#include "SundanceExpr.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"
#include "PlayaHandleable.hpp"
#include "SundanceFunctionalEvaluatorBase.hpp"

namespace Sundance
{

     // forward declaration of the mesh class
     class Mesh;

/**
 * The base class for parameterized curves<br>
 * This class is mean to used in the Adaptive Cell Integration. <br>
 * This class optionally can passed to one integration and the curves defines
 * where the integral is defined, and where it is not defined. <br>
 * The class defines the basic interface methods, which are needed in the
 * Adaptive Cell Integration.
 */

class CurveBase : public Playa::Handleable<CurveBase>
{
public:

	CurveBase(){;}

	/** Base constructor
	 * @param dim , the dimension of the curve, where it is defined
	 * @param alpha1 , the integration coefficient for inside the domain
	 * @param alpha2 , the integration coefficient for outside the domain
	 * @param flipD, to flip the domains for this curve*/
	CurveBase(int dim , double alpha1 , double alpha2 , bool flipD = false):_dimension(dim),
	_alpha1(alpha1),_alpha2(alpha2) , flipDomains_(1.0){
		if (flipD) flipDomains_ = -1.0;
	}

	virtual ~CurveBase(){;}

	/** This function should be used later in the optimization process
	 * TODO: This interface should be extended later
	 * @return Expr The parameters of the curve which uniquely defines the curve*/
	virtual Expr getParams() const = 0 ;

	/** update the state of the curve of the control point were changed*/
	virtual void update() {;}

	/** function to return the control points of a given curve */
	virtual Array<Point>& getControlPoints() {
		TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error, " getControlPoints() method is not overwritten ");
	}

	/** function which will be called for curve/surface Integral evaluation
	 * @param vars [IN] coordinates
	 * @param f [OUT] the output of the expression
	 * @param scalarFieldIndex [IN] scalar field index */
	virtual void eval0(const double* vars, double* f , int scalarFieldIndex ) const {
		TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error, " eval0() method is not overwritten ");
	}

	/** adds new scalar field to the interface
	 * @param fieldName [IN] the new name of the field
	 * @param initialValue [IN] initial value of the scalar field
	 * @return the scalar field index */
	virtual int addNewScalarField(std::string fieldName , double initialValue){
		TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error, " addNewScalarField() method is not overwritten ");
	}

	/** returns one array of doubles which contain the values of the scalar field (for each point, nodal basis)
	 * @param scalarFieldIndex [IN] the index
	 * @return array with the scalar field values */
	virtual Array<double>& getScalarFieldValues(int scalarFieldIndex) {
		TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error, " getScalarFieldValues() method is not overwritten ");
	}

	/** function to set the values of one specified scalar field <br>
	 * it is important that the Functional will be defined on a curve, otherwise the value will stay zero
	 * @param scalarFunctional [IN] , the scalar functions*/
	virtual void setSpaceValues(const FunctionalEvaluatorBase& scalarFunctional , int fieldIndex ){
		TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error, " setSpaceValues() method is not overwritten ");
	}

	/** function to record the function values along one curve. <br>
	 * This function will be called from the integration routine.
	 * @param mesh [IN] the mesh (only for consistency test if the curve mesh is the same as the integral mesh)
	 * @param maxCellLID [IN] the max dim cell LID on which is currently integrated
	 * @param nQuad [IN] number of quadrature points (and implicitly evaluation point)
	 * @param coeffPtr [IN] the values of the expression at the quadrature points
	 * @param quadPts [IN] the evaluation points in the reference coordinates */
	virtual void addEvaluationPointValues(const Mesh& mesh ,
			int maxCellLID , int nQuad ,
			const double* coeffPtr ,
			const Array<Point>& quadPts) const { /*std::cout << " dummy addEvaluationPointValues" << std::endl*/;}

	/**
	 * Return the integration parameters
	 * @param alphas, integration parameters, first one corresponds to Equation > 0
	 */
	void getIntegrationParams(Array<double>& alphas) const {
		alphas.resize(2);
		alphas[0] = _alpha1;
		alphas[1] = _alpha2;
	}

	/** return the integration coefficient for the inner */
	double getAlpha1() const { return _alpha1; }

	/** return the integration coefficient for the outer */
	double getAlpha2() const { return _alpha2; }

	/** return the dimension of the curve , in which it is defined */
	int getCurveDim() const { return _dimension; }

	/**
	 * @param evaluationPoint the point where we want the alpha integration parameter <br>
	 * The parameter alpha is used for the Adaptive Cell Integration
	 * @return double the alpha parameter for the integration */
	double integrationParameter(const Point& evaluationPoint) const{// here we return the integration parameter
		// this function can be overwritten
		if (curveEquation(evaluationPoint) > 0 )
			return _alpha1;
		else
			return _alpha2;
	}

	/** If we want to flip the fictitious with the real computational domain then call this method */
	void flipDomains() const {
      if ( flipDomains_ > 0){
    	  flipDomains_ = -1.0;
      } else {
    	  flipDomains_ = 1.0;
      }
	}

	/**
	 * Returns the value of the parameterized curve at one given location
	 * @param evaluationPoint the point where we want the alpha integration parameter <br>
	 * @return double the value of the curve equation at the evaluation point  */
	double curveEquation(const Point& evaluationPoint) const {
		return flipDomains_ * curveEquation_intern(evaluationPoint);
	}

protected:

	/**
	 * This function should be implemented
	 * @param evaluationPoint the point where we want the alpha integration parameter <br>
	 * @return double the value of the curve equation at the evaluation point  */
	virtual double curveEquation_intern(const Point& evaluationPoint) const = 0 ;

public:

	/**
	 * This function is important for nonstructural mesh integration.<br>
	 * The function returns the intersection points with a given line (in 2D and in 3D)<br>
	 * The line is defined with a starting point and an end point
	 * @param startEdgePoint , the start point of the line
	 * @param endEdgePoint , the end point of the line
	 * @param nrPoints , number of resulted (intersected) points
	 * @param result , the resulted points of intersections */
	virtual void returnIntersectPoints(const Point& startEdgePoint, const Point& endEdgePoint,
			                      int& nrPoints ,Array<Point>& result) const = 0 ;

	/**
	 * As above, but, instead of coordinates, returns intersection values t in [0,1]
	 * and the line is defined by "start + t * (end-start)"
	 */
	virtual void returnIntersect(const Point& startEdgePoint, const Point& endEdgePoint,
				                      int& nrPoints ,Array<double>& result) const = 0 ;

	/** In some constructors we need to pass a Parametrized curve even if there is any
	 * so we need sometimes to pass a dummy curve, and this function shows which curve is dummy and which is not */
	virtual bool isCurveValid() const = 0;

	/** function which shows if some integral schould be calculated along the curve <br>
	 * default this returns false , only one wrapper class should overwrite this and return false */
	virtual bool isCurveIntegral() const { return false; }

	/** In the case of simple geometries the geometry it can be transformed to a polygon, which
	 * reflects the original geometry
	 * @param mesh
	 * @param resolution , the global resolution */
	virtual const RCP<CurveBase> getPolygon(const Mesh& mesh , double resolution) const {
		TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error, " getPolygon() method is not overwritten ");
		return rcp((CurveBase*)0);
	}

	/** e.g. the polygon needs the mesh for additional informations
	 * empty implementation in the base class
	 * @param mesh */
	virtual void setMesh(const Mesh& mesh) {;}

	/** is case of wrapper objects this should return the real object behind the wrapper
	 * e.g. such wrapper is the ParametrizedCurveIntegral class */
	virtual const CurveBase* getUnderlyingCurveObj() const { return this; }

	/** Writes the geometry into a VTK file for visualization purposes
	 * @param filename */
	virtual void writeToVTK(const std::string& filename) const {
		TEUCHOS_TEST_FOR_EXCEPTION( true , std::runtime_error, " writeToVTK() method is not overwritten ");
	}

protected:

	/** The dimension where the curve is defined */
	int _dimension;

	/** The integration parameter when the curve equation is positive*/
	double _alpha1;

	/** The integration parameter when the curve equation is negative*/
	double _alpha2;

	/** by default this value is 1.0 , if we want to flip the domains then -1.0*/
	mutable double flipDomains_;
};

} // end from namespace
#endif /* SUNDANCECURVEBASE_H_ */
