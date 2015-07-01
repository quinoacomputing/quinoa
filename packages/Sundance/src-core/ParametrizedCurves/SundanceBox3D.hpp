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

#ifndef SUNDANCEBOX3D_H_
#define SUNDANCEBOX3D_H_

#include "SundanceCurveBase.hpp"

namespace Sundance
{

/** This class models a rectangle (Box) in 3D.<br>
 * This is one of the simplest curves */
class Box3D: public Sundance::CurveBase
{
public:

	/** The constructor of a 3D
	 * @param px the lowest X coordinate of the box
	 * @param py the lowest Y coordinate of the box
	 * @param pz the lowest Z coordinate of the box
	 * @param ox offset of the box in X direction
	 * @param oy offset of the box in Y direction
	 * @param oz offset of the box in Z direction
	 * @param a1 alpha1 coefficient for FCM
	 * @param a2 alpha2 coefficient for FCM */
	Box3D(double px, double py, double pz,
		  double ox, double oy, double oz, double a1, double a2, bool flipD = false);

	virtual ~Box3D();

	/** @return Expr The parameters of the curve which uniquely defines the curve*/
	virtual Expr getParams() const;

protected:
	/**
	 * This function should be implemented
	 * @param evalPoint point where the box equation is evaluated <br>
	 * @return double the value of the curve equation at the evaluation point  */
	virtual double curveEquation_intern(const Point& evalPoint) const;
public:

	/**
	 * This function is important for nonstructural mesh integration.<br>
	 * The function returns the intersection points with a given line <br>
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

private:

	/** */
	double px_;

	/** */
	double py_;

	/** */
	double pz_;

	/** */
	double ox_;

	/** */
	double oy_;

	/** */
	double oz_;

};

}
#endif /* SUNDANCEBOX3D_H_ */
