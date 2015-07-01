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

#include "SundanceLine2D.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"
#include "SundancePolygon2D.hpp"

using namespace Sundance;

Line2D::Line2D(double slope, double b, double a1, double a2, bool flipD ) :
	CurveBase(1, a1, a2, flipD), slope_(slope), b_(b)
{
}

Line2D::~Line2D()
{
}

Expr Line2D::getParams() const
{
	// return the parameters of the 2D circle
	return Expr(List(slope_, b_ ));
}

double Line2D::curveEquation_intern(const Point& evalPoint) const
{
	TEUCHOS_TEST_FOR_EXCEPTION(evalPoint.dim() != 2, std::runtime_error,
			"Line2D::curveEquation() evaluation point dimension must be 2");

	// y = slope *x + b
	return (slope_ * evalPoint[0] + b_ - evalPoint[1]);
}

void Line2D::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
		Array<Point>& result) const
{
	Array<double> t;
	returnIntersect(start, end, nrPoints, t);

	result.resize(2);

	// Return coordinates instead of t values
	for (int i = 0; i < nrPoints; i++)
	{
		result[i] = start + t[i] * (end - start);
	}
}

void Line2D::returnIntersect(const Point& start, const Point& end, int& nrPoints, Array<
		double>& result) const
{
	 // we calculate two different lines which have the same slope , and take the free term
     double b1 = start[1] - slope_* start[0];
     double b2 = end[1] - slope_* end[0];

     nrPoints = 0;

     // test if b_ is between b1 and b2
     if ( ((b_ >= b1) && (b_ <= b2)) || ((b_ <= b1) && (b_ >= b2)) ){
    	 result.resize(1);
    	 result[0] = fabs((b1-b_)/(b1-b2));
    	 nrPoints = 1;
     }
}

const RCP<CurveBase> Line2D::getPolygon(const Mesh& mesh , double resolution) const {

	int verb = 0;

	/* y = slope*x + b */

	// dummy implementation
	// have always 100 points, and x goes from 0 to 1
	int nrPoints = 100;
	SUNDANCE_MSG3( verb , " Line2D::getPolygon nrPoints = " << nrPoints );
	Array<Point> points(nrPoints);


	for (int pI = 0 ; pI < nrPoints ; pI++){
		Point p( (-5.0 +10.0*(double)(pI)/(double)(nrPoints-1)) , slope_*(-5.0 +10.0*(double)(pI)/(double)(nrPoints-1)) + b_);
		SUNDANCE_MSG3( verb , " Line2D::getPolygon add pint " << pI << " p=" << p );
		points[pI] = p;
	}

	// return the polygon
	return rcp(new Polygon2D( mesh , points , _alpha1 , _alpha2 , false, (flipDomains_<0)));
}

