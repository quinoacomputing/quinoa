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

#include "SundancePlane3D.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"

using namespace Sundance;

Plane3D::Plane3D(double a, double b, double c, double a1, double a2, bool flipD ) :
	CurveBase(2, a1, a2, flipD), a_(a) , b_(b) , c_(c)
{
}

Plane3D::~Plane3D()
{
}

Expr Plane3D::getParams() const
{
	// return the parameters of the 3D plane z = a*x + b*y + c
	return Expr(List( a_ , b_ , c_ ));
}

double Plane3D::curveEquation_intern(const Point& evalPoint) const
{
	TEUCHOS_TEST_FOR_EXCEPTION(evalPoint.dim() != 3, std::runtime_error,
			"Plane3D::curveEquation() evaluation point dimension must be 3");

	// z = a*x + b*y + c
	return ( a_*evalPoint[0] + b_*evalPoint[1] + c_ -  evalPoint[2]);
}

void Plane3D::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
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

void Plane3D::returnIntersect(const Point& start, const Point& end, int& nrPoints, Array<
		double>& result) const
{
	 // we calculate two different planes which have the same slope , and take the free term
	// z = a*x + b*y + c
    double c1 = start[2] - a_* start[0] - b_* start[1];
    double c2 = end[2] - a_* end[0] - b_* end[1];

    nrPoints = 0;

    // test if c_ is between c1 and c2
    if ( ((c_ >= c1) && (c_ <= c2)) || ((c_ <= c1) && (c_ >= c2)) ){
    	result.resize(1);
    	result[0] = fabs((c1-c_)/(c1-c2));
    	nrPoints = 1;
    }
}

