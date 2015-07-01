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

#ifndef SUNDANCECIRCLE_H_
#define SUNDANCECIRCLE_H_

#include "SundanceCurveBase.hpp"

namespace Sundance
{

/** This class models a circle in 2D.<br>
 * This is one of the simplest curves */
class Circle: public Sundance::CurveBase
{
public:

	/** The constructor of the circle in 2D
	 * @param centerx , the X coordinate of the center
	 * @param centery, the Y coordinate of the center
	 * @param radius
	 * @param a1 alpha1 coefficient for FCM
	 * @param a2 alpha2 coefficient for FCM */
	Circle(double centerx, double centery, double radius, double a1, double a2, bool flipD = false);

	virtual ~Circle();

	/** @return Expr The parameters of the curve which uniquely defines the curve*/
	virtual Expr getParams() const;

protected:
	/**
	 * This function should be implemented
	 * @param evalPoint point where the circle equation is evaluated <br>
	 * @return double the value of the curve equation at the evaluation point  */
	virtual double curveEquation_intern(const Point& evalPoint) const;
public:

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

	/** In the case of simple geometries the geometry it can be transformed to a polygon, which
	 * reflects the original geometry
	 * @param mesh
	 * @param resolution , the global resolution */
	virtual const RCP<CurveBase> getPolygon(const Mesh& mesh ,double resolution) const ;

private:

	/** */
	double _centerx;

	/** */
	double _centery;

	/** */
	double _radius;

};

}
#endif /* SUNDANCECIRCLE_H_ */
