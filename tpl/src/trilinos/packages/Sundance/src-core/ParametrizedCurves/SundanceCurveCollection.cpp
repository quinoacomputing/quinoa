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

#include "SundanceCurveCollection.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"

using namespace Sundance;

CurveCollection::CurveCollection(double a1, double a2,int dim) :
	CurveBase(dim, a1, a2, false) , nrCurves_(0)
{
}

CurveCollection::~CurveCollection()
{
}

Expr CurveCollection::getParams() const
{
	// no parameters
	return Expr(List(0.0));
}

double CurveCollection::curveEquation_intern(const Point& evalPoint) const
{
	int verb = 0;

	// calculate the distance compared to the middle point
	double distX = 1e+100;

	SUNDANCE_OUT(verb > 3, " CurveCollection::curveEquation before:" << evalPoint << " is: " << distX);
	if ( nrCurves_ > 0) distX = curves_[0]->curveEquation(evalPoint);
	for (int c = 1 ; c < nrCurves_; c++){
		double tmp = curves_[c]->curveEquation(evalPoint);
		distX = (tmp > distX) ? distX : tmp;
	}

	SUNDANCE_OUT(verb > 3, " CurveCollection::curveEquation for:" << evalPoint << " is: " << distX);

	return distX;
}

void CurveCollection::returnIntersect(const Point& start, const Point& end, int& nrPoints,
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

void CurveCollection::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
		Array<Point>& result) const
{
    int verb = 0;

    result.resize(0); nrPoints = 0;
	for (int c = 0 ; c < nrCurves_; c++){
		Array<Point> pnts;
		int nrP;
		curves_[c]->returnIntersectPoints( start, end, nrP, pnts);
		for (int p = 0 ; p < nrP ; p++){
			result.append( pnts[p] );
			nrPoints++;
		}
	}

	SUNDANCE_OUT(verb > 3, " CurveCollection::returnIntersectPoints END , nrPoints:" << nrPoints
			<< " , start:" << start << " , end:" << end);
}
