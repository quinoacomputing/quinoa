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

#include "SundanceBox3D.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"

using namespace Sundance;

Box3D::Box3D(double px, double py, double pz,
		     double ox, double oy, double oz, double a1, double a2, bool flipD ) :
	CurveBase(2, a1, a2,flipD), px_(px), py_(py), pz_(pz), ox_(ox), oy_(oy), oz_(oz)
{
}

Box3D::~Box3D()
{
}

Expr Box3D::getParams() const
{
	// return the parameters of the box
	return Expr(List(px_, py_, py_, ox_ , oy_, oz_));
}

double Box3D::curveEquation_intern(const Point& evalPoint) const
{
	int verb = 0;
	TEUCHOS_TEST_FOR_EXCEPTION(evalPoint.dim() != 3, std::runtime_error,
			"Box3D::curveEquation() evaluation point dimension must be 3");

	// calculate the distance compared to the middle point
	double distX =  fabs(px_ + 0.5*ox_ - evalPoint[0]) - 0.5*ox_;
	double distY =  fabs(py_ + 0.5*oy_ - evalPoint[1]) - 0.5*oy_;
	double distZ =  fabs(pz_ + 0.5*oz_ - evalPoint[2]) - 0.5*oz_;
	distX = (distX > distY) ? distX : distY ;
	distX = (distX > distZ) ? distX : distZ ;

	SUNDANCE_OUT(verb > 3, " Box3D::curveEquation for:" << evalPoint << " is: " << distX);

	return distX;
}

void Box3D::returnIntersect(const Point& start, const Point& end, int& nrPoints,
		Array<double>& result) const
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

void Box3D::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
		Array<Point>& result) const
{
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
				"Box3D::returnIntersectPoints() not implemented yet");
	// todo: implement this
}

