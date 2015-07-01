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

#include "SundanceCircle.hpp"
#include "SundancePoint.hpp"
#include "SundancePolygon2D.hpp"
#include "SundanceDefs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;

Circle::Circle(double centerx, double centery, double radius, double a1,
		double a2, bool flipD ) :
	CurveBase(1, a1, a2, flipD), _centerx(centerx), _centery(centery), _radius(radius)
{
}

Circle::~Circle()
{
}

Expr Circle::getParams() const
{
	// return the parameters of the 2D circle
	return Expr(List(_centerx, _centery, _radius));
}

double Circle::curveEquation_intern(const Point& evalPoint) const
{
	TEUCHOS_TEST_FOR_EXCEPTION(evalPoint.dim() != 2, std::runtime_error,
			"Circle::curveEquation() evaluation point dimension must be 2");

	Point center(_centerx, _centery);

	// the circle equation is (x-cx)^2 + (y-cy)^2 - r^2 = 0
	return (((evalPoint - center) * (evalPoint - center)) - _radius * _radius);
}

void Circle::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
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

void Circle::returnIntersect(const Point& start, const Point& end, int& nrPoints, Array<
		double>& result) const
{
	// direction vector; line == start + t * vec
	Point vec = end - start;

	// center point of circle
	Point center(_centerx, _centery);

	nrPoints = 0;
	result.resize(2);

	// intersection of circle and line leads to
	// norm(start + t * vec - center) == r
	// i.e. following coefficients of quadratic equation
	double a = vec * vec;
	double b = 2 * ((start - center) * vec);
	double c = ((start - center) * (start - center)) - _radius * _radius;

	double det = b * b - 4 * a * c;
	if (det < 0.0)
		return; // no solution
	det = sqrt(det);

	// numerically more favorable than usual usage of solution formula
	int sign = b < 0.0 ? -1 : 1;
	double q = -0.5 * (b + sign * det);

	// first solution
	double t = c / q;

	// Only consider intersection points between 'start' and 'end'
	if (t >= 0.0 && t <= 1.0)
		result[nrPoints++] = t;

	// second solution (only considered if different from first)
	if (det < 1.e-14)
		return;
	t = q / a;
	if (t >= 0.0 && t <= 1.0)
	{
		// Sorting: first solution is nearest to 'start'
		if (t < result[0])
		{
			result[nrPoints++] = result[0];
			result[0] = t;
		}
		else
		{
			result[nrPoints++] = t;
		}
	}
}

const RCP<CurveBase> Circle::getPolygon(const Mesh& mesh , double resolution) const {

	int verb = 0;
	// 2*pi*r/h will give the angle
	double average_angle = resolution/(2.0*3.14159265358979*_radius);

	int nrPoints = ::ceil (2.0*3.14159265358979/average_angle);
	double stepAngle = (2.0*3.14159265358979/(double)nrPoints);

	SUNDANCE_MSG3( verb , " Circle::getPolygon average_angle=" << average_angle << " nrPoints = " << nrPoints << " stepAngle=" << stepAngle);
	Array<Point> points(nrPoints);
	for (int pI = 0 ; pI < nrPoints ; pI++){
		Point p(_centerx + _radius*::cos(stepAngle*(double)pI), _centery + _radius*::sin(stepAngle*(double)pI));
		SUNDANCE_MSG3( verb , " Circle::getPolygon add pint " << pI << " p=" << p );
		points[pI] = p;
	}

	// return the polygon
	return rcp(new Polygon2D( mesh , points , _alpha1 , _alpha2 , true , (flipDomains_ < 0) ));
}
