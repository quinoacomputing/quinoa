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

#include "SundanceFeketeBrickQuadrature.hpp"
#include "SundanceFeketeQuadQuadrature.hpp"
#include "SundanceFeketeQuadrature.hpp"
#include "SundanceFeketeTriangleQuadrature.hpp"
#include "SundanceGaussLobatto1D.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include <stack>

using namespace Sundance;
using namespace Teuchos;

/* declare LAPACK subroutines */
extern "C"
{
/* matrix-vector multiplication */
void dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda,
		double *x, int *incx, double *beta, double *y, int *incy);
}

/* constructor */
FeketeQuadrature::FeketeQuadrature(int order) :
	QuadratureFamilyBase(order)
{
	_hasBasisCoeffs = false;
}

XMLObject FeketeQuadrature::toXML() const
{
	XMLObject rtn("FeketeQuadrature");
	rtn.addAttribute("order", Teuchos::toString(order()));
	return rtn;
}

void FeketeQuadrature::getLineRule(Array<Point>& quadPoints,
		Array<double>& quadWeights) const
{
	int p = order() + 3;
	p = p + (p % 2);
	int n = p / 2;

	quadPoints.resize(n);
	quadWeights.resize(n);

	GaussLobatto1D q1(n, 0.0, 1.0);

	for (int i = 0; i < n; i++)
	{
		quadWeights[i] = q1.weights()[i];
		quadPoints[i] = Point(q1.nodes()[i]);
	}
}

void FeketeQuadrature::getTriangleRule(Array<Point>& quadPoints,
		Array<double>& quadWeights) const
{
	Array<double> x;
	Array<double> y;
	Array<double> w;

	FeketeTriangleQuadrature::getPoints(order(), w, x, y);
	quadPoints.resize(w.length());
	quadWeights.resize(w.length());
	for (int i = 0; i < w.length(); i++)
	{
		quadWeights[i] = 0.5 * w[i];
		quadPoints[i] = Point(x[i], y[i]);
	}
}

void FeketeQuadrature::getQuadRule(Array<Point>& quadPoints,
		Array<double>& quadWeights) const
{
	Array<double> x;
	Array<double> y;
	Array<double> w;

	FeketeQuadQuadrature::getPoints(order(), w, x, y);
	quadPoints.resize(w.length());
	quadWeights.resize(w.length());
	for (int i = 0; i < w.length(); i++)
	{
		quadWeights[i] = w[i];
		quadPoints[i] = Point(x[i], y[i]);
	}
}

void FeketeQuadrature::getBrickRule(Array<Point>& quadPoints,
		Array<double>& quadWeights) const
{
	Array<double> x;
	Array<double> y;
	Array<double> z;
	Array<double> w;

	FeketeBrickQuadrature::getPoints(order(), w, x, y, z);
	quadPoints.resize(w.length());
	quadWeights.resize(w.length());
	for (int i = 0; i < w.length(); i++)
	{
		quadWeights[i] = w[i];
		quadPoints[i] = Point(x[i], y[i], z[i]);
	}
}

void FeketeQuadrature::getAdaptedWeights(const CellType& cellType, int cellDim,
		int cellLID, int facetIndex, const Mesh& mesh, const ParametrizedCurve&
		globalCurve, Array<Point>& quadPoints, Array<double>& quadWeights,
		bool& weightsChanged) const
{
	//int verb = 4;
	Tabs tabs;

	// Cache lookup
	/*if (mesh.IsSpecialWeightValid() && mesh.hasSpecialWeight(cellDim, cellLID))
	{
		weightsChanged = true;
		mesh.getSpecialWeight(cellDim, cellLID, quadWeights);
		SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::getAdaptedWeights Found cached weights for cell LID " << cellLID)
		return;
	}*/

	// Maximal cell type
	CellType maxCellType = mesh.cellType(mesh.spatialDim());

	switch (maxCellType)
	{
	// We have a triangle based mesh
	case TriangleCell:
	{
		switch (cellType)
		{
		case TriangleCell:
		{
			getAdaptedTriangleWeights(cellLID, mesh, globalCurve, quadPoints,
					quadWeights, weightsChanged);
			break;
		}
		default:
#ifndef TRILINOS_7
			SUNDANCE_ERROR ("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for submaximal cell " << maxCellType << " in a triangle mesh")
			;
#else
			SUNDANCE_ERROR7("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for submaximal cell " << maxCellType << " in a triangle mesh");
#endif
		}
		break;
	}
	case QuadCell:
	{
		switch(cellType)
		{
		case QuadCell:
		{
			getAdaptedQuadWeights(cellLID, mesh, globalCurve, quadPoints,
				quadWeights, weightsChanged);
			break;
		}
		default:
#ifndef TRILINOS_7
		SUNDANCE_ERROR ("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for submaximal cell " << maxCellType << " in a triangle mesh")
		;
#else
		SUNDANCE_ERROR7("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for submaximal cell " << maxCellType << " in a triangle mesh");
#endif
		}
		break;
	}
	default:
#ifndef TRILINOS_7
		SUNDANCE_ERROR("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for cell type " << cellType)
		;
#else
		SUNDANCE_ERROR7("getAdaptedWeights rule for ACI(Adaptive Cell Integration) not available for cell type " << cellType);
#endif
	}
}

void FeketeQuadrature::getAdaptedTriangleWeights(int cellLID, const Mesh& mesh,
		const ParametrizedCurve& globalCurve, Array<Point>& quadPoints, Array<
				double>& quadWeights, bool& weightsChanged) const
{
	int verb = 0;
	Tabs tabs;

	// ToDo: What tolerance?
	double tol = 1.e-8;
	weightsChanged = true;

	// For pushForward() we need the cellLID in an array
	Array<int> cellLIDs(1);
	cellLIDs[0] = cellLID;

	// Intersection points (parameters) with the edges of triangle
	Array<double> ip01, ip02, ip12;
	int nr_ip01, nr_ip02, nr_ip12;

	Array<double> originalWeights = quadWeights;
	ParametrizedCurve noCurve = new DummyParametrizedCurve();
	Array<double> integrals(originalWeights.size());
	for (int i = 0; i < integrals.size(); i++)
		integrals[i] = 0.0;

	// Create stack and put the initial triangle on stack
	std::stack<Point> s;
	std::stack<double> sa;
	s.push(Point(0.0, 0.0));
	s.push(Point(1.0, 0.0));
	s.push(Point(0.0, 1.0));
	sa.push(0.5);
	Array<Point> nodes;
	Array<Point> nodesref(3);

	// Number of investigated triangles
	unsigned int nTris = 0;

	while (s.size() >= 3)
	{
		// Increment number of triangles investigated
		++nTris;

		SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::getAdaptedTriangleWeights Current number of triangles on stack: " << sa.size());

		// Integration results
		bool refine = false;
		Array<double> results1(originalWeights.size());
		Array<double> results2(originalWeights.size());

		// Fetch the coordinates of current triangle from stack
		double area = sa.top();
		sa.pop();
		for (int i = 2; i >= 0; i--)
		{
			nodesref[i] = s.top();
			s.pop();
		}
		mesh.pushForward(2, cellLIDs, nodesref, nodes);

		// Intersection points (parameters) with the edges of triangle
		globalCurve.returnIntersect(nodes[0], nodes[1], nr_ip01, ip01);
		globalCurve.returnIntersect(nodes[0], nodes[2], nr_ip02, ip02);
		globalCurve.returnIntersect(nodes[1], nodes[2], nr_ip12, ip12);

		//  Determine kind of intersection
		Point x, xref;
		Point vec1, vec1ref;
		Point vec2, vec2ref;
		bool xInOmega = false;

		// 'No intersection' is also handled here
		if ((nr_ip01 + nr_ip02 + nr_ip12 == 0) || (nr_ip01 == 0 && nr_ip02 == 0
				&& nr_ip12 == 1) || (nr_ip01 == 1 && nr_ip02 == 1 && nr_ip12
				== 0))
		{
			x = nodes[0];
			vec1 = nodes[2] - nodes[0];
			vec2 = nodes[1] - nodes[0];
			xref = nodesref[0];
			vec1ref = nodesref[2] - nodesref[0];
			vec2ref = nodesref[1] - nodesref[0];
		}
		else if ((nr_ip01 == 0 && nr_ip02 == 1 && nr_ip12 == 0) || (nr_ip12
				== 1 && nr_ip01 == 1 && nr_ip02 == 0))
		{
			x = nodes[1];
			vec1 = nodes[0] - nodes[1];
			vec2 = nodes[2] - nodes[1];
			xref = nodesref[1];
			vec1ref = nodesref[0] - nodesref[1];
			vec2ref = nodesref[2] - nodesref[1];
		}
		else if ((nr_ip01 == 1 && nr_ip02 == 0 && nr_ip12 == 0) || (nr_ip02
				== 1 && nr_ip12 == 1 && nr_ip01 == 0))
		{
			x = nodes[2];
			vec1 = nodes[1] - nodes[2];
			vec2 = nodes[0] - nodes[2];
			xref = nodesref[2];
			vec1ref = nodesref[1] - nodesref[2];
			vec2ref = nodesref[0] - nodesref[2];
		}
		else
		{
			refine = true;
		}

		// If we found a case we can deal with, calculate the integrals
		if (!refine)
		{
			if (globalCurve.curveEquation(x) < 0)
				xInOmega = true;
			integrateRegion(TriangleCell, 2, order() + 1, x, xref, vec1, vec2,
					vec1ref, vec2ref, globalCurve, results1);
			if (results1.size() == 0)
			{
				SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::getAdaptedTriangleWeights Unknown error in integration. Refine.")
				refine = true;
			}
			else
			{
				integrateRegion(TriangleCell, 2, order() + 2, x, xref, vec1,
						vec2, vec1ref, vec2ref, globalCurve, results2);

				// Check results
				for (int i = 0; i < results1.size(); i++)
				{
					if (::fabs(results2[i] - results1[i]) > sqrt(area * tol))
					{
						SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::getAdaptedTriangleWeights Insufficient accuracy. Refine.")
						refine = true;
						break;
					}
				}
			}
		}

		// Refine (either because of lacking accuracy or strange intersections)
		if (refine)
		{
			// Divide triangle into 4 smaller ones (introduce new nodes at the
			// center of each edge)
			s.push(nodesref[0]);
			s.push((nodesref[0] + nodesref[1]) / 2);
			s.push((nodesref[0] + nodesref[2]) / 2);
			sa.push(area / 4);
			s.push((nodesref[0] + nodesref[1]) / 2);
			s.push(nodesref[1]);
			s.push((nodesref[1] + nodesref[2]) / 2);
			sa.push(area / 4);
			s.push((nodesref[0] + nodesref[2]) / 2);
			s.push((nodesref[1] + nodesref[2]) / 2);
			s.push(nodesref[2]);
			sa.push(area / 4);
			s.push((nodesref[0] + nodesref[2]) / 2);
			s.push((nodesref[0] + nodesref[1]) / 2);
			s.push((nodesref[1] + nodesref[2]) / 2);
			sa.push(area / 4);
			continue;
		}

		// We found an accurate solution to that part of the triangle!
		if (xInOmega)
		{
			for (int i = 0; i < integrals.size(); i++)
				integrals[i] += results2[i];
		}
		else
		{
			integrateRegion(TriangleCell, 2, order() + 1, x, xref, vec1, vec2,
					vec1ref, vec2ref, noCurve, results1);
			for (int i = 0; i < integrals.size(); i++)
				integrals[i] += (results1[i] - results2[i]);
		}
	} // Stack depleted here! We're done.

	SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::getAdaptedTriangleWeights Calculations in cell " << cellLID << " completed. " << nTris << " triangle(s) used.");

	Array<double> alphas;
	globalCurve.getIntegrationParams(alphas);
	for (int i = 0; i < quadWeights.size(); i++)
		quadWeights[i] = alphas[1] * integrals[i] + alphas[0]
				* (originalWeights[i] - integrals[i]);
	mesh.setSpecialWeight(2, cellLID, quadWeights);
}

void FeketeQuadrature::getAdaptedQuadWeights(int cellLID, const Mesh& mesh,
		const ParametrizedCurve& globalCurve, Array<Point>& quadPoints, Array<
				double>& quadWeights, bool& weightsChanged) const
{
	// Verbosity
	int verb = 0;
	Tabs tabs;

	// ToDo: What tolerance?
	double tol = 1.e-8;

	// Initialization
	weightsChanged = true;
	Array<double> originalWeights = quadWeights;

	ParametrizedCurve noCurve = new DummyParametrizedCurve();
	Array<double> integrals(originalWeights.size());
	for (int i = 0; i < integrals.size(); i++)
		integrals[i] = 0.0;

	// For pushForward() we need the cellLID in an array
	Array<int> cellLIDs(1);
	cellLIDs[0] = cellLID;
	SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::getAdaptedQuadWeights: Current cell LID " << cellLID)

	//        IV
	//   2 __________ 3
	//    |          |
	//    |          |
	// II |          | III
	//    |__________|
	//   0            1
	//         I

	Array<int> nIntEdge(4);
	Array<Array<double> > IntEdgePars(4);
	int edgeIndex[4][2] = { {0,1} , {0,2} , {1,3} , {2,3} };

	// Create stack and put the initial quad
	std::stack<Point> s;
	std::stack<double> sa;
	s.push(Point(0.0, 0.0));
	s.push(Point(1.0, 0.0));
	s.push(Point(0.0, 1.0));
	s.push(Point(1.0, 1.0));
	sa.push(1.);
	Array<Point> nodes;
	Array<Point> nodesref(4);

	// Number of investigated quads
	unsigned int nQuads = 0;

	while (s.size() >= 4)
	{
		// Increment number of quads investigated
		++nQuads;
		SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::getAdaptedQuadWeights: Current number of quads on stack: " << sa.size())

		// Integration results
		bool refine = false;
		bool xInOmega = false;
		Array<double> results1(originalWeights.size());
		Array<double> results2(originalWeights.size());

		// Fetch the coordinates of current quad from stack
		double area = sa.top();
		sa.pop();
		for (int i = 3; i >= 0; i--)
		{
			nodesref[i] = s.top();
			s.pop();
		}
		mesh.pushForward(2, cellLIDs, nodesref, nodes);

		// Intersection points (parameters) with the edges of triangle
		for (int jj = 0 ; jj < 4 ; jj++ ) {
			globalCurve.returnIntersect(nodes[edgeIndex[jj][0]], nodes[edgeIndex[jj][1]], nIntEdge[jj], IntEdgePars[jj] );
	/*		if ( nIntEdge[jj] == 0 ) {
				x = nodes[edgeIndex[jj][0]];
				vec1 =
				vec2 =
			}
	*/	}

		//  Determine kind of intersection
		Point x, xref;
		Point vec1, vec1ref;
		Point vec2, vec2ref;

		// Baseline is not intersected
		if ( nIntEdge[0] == 0 )
		{
			x = nodes[0];
			vec1 = nodes[1] - nodes[0];
			vec2 = nodes[2] - nodes[0];
			xref = nodesref[0];
			vec1ref = nodesref[1] - nodesref[0];
			vec2ref = nodesref[2] - nodesref[0];
		}
		else if ( nIntEdge[1] == 0 )
		{
			x = nodes[2];
			vec1 = nodes[0] - nodes[2];
			vec2 = nodes[3] - nodes[2];
			xref = nodesref[2];
			vec1ref = nodesref[0] - nodesref[2];
			vec2ref = nodesref[3] - nodesref[2];
		}
		else if ( nIntEdge[2] == 0 )
		{
			x = nodes[1];
			vec1 = nodes[3] - nodes[1];
			vec2 = nodes[0] - nodes[1];
			xref = nodesref[1];
			vec1ref = nodesref[3] - nodesref[1];
			vec2ref = nodesref[0] - nodesref[1];
		}
		else
		{
			refine = true;
		}

		// If we found a case we can deal with, calculate the integrals
		if (!refine)
		{
			if (globalCurve.curveEquation(x) < 0)
				xInOmega = true;
			integrateRegion(QuadCell, 2, order() + 1, x, xref, vec1, vec2,
					vec1ref, vec2ref, globalCurve, results1);
			if (results1.size() == 0)
			{
				SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::getAdaptedQuadWeights Unknown error in integration. Refine.")
				refine = true;
			}
			else
			{
				integrateRegion(QuadCell, 2, order() + 2, x, xref, vec1,
						vec2, vec1ref, vec2ref, globalCurve, results2);

				// Check results
				double reltol = ::sqrt(area * tol);
				for (int i = 0; i < results1.size(); i++)
				{
					if (::fabs(results2[i] - results1[i]) > reltol)
					{
						SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::getAdaptedQuadWeights Insufficient accuracy. Refine.")
						refine = true;
						break;
					}
				}
			}
		}

		// Refine (either because of lacking accuracy or strange intersections)
		if (refine)
		{
			// Divide quad into 4 smaller ones (introduce new nodes at the
			// center of each edge)
			s.push(nodesref[0]);
			s.push((nodesref[0] + nodesref[1]) / 2);
			s.push((nodesref[0] + nodesref[2]) / 2);
			s.push((nodesref[0] + nodesref[3]) / 2);
			sa.push(area / 4);
			s.push((nodesref[0] + nodesref[1]) / 2);
			s.push(nodesref[1]);
			s.push((nodesref[1] + nodesref[2]) / 2);
			s.push((nodesref[1] + nodesref[3]) / 2);
			sa.push(area / 4);
			s.push((nodesref[0] + nodesref[2]) / 2);
			s.push((nodesref[1] + nodesref[2]) / 2);
			s.push(nodesref[2]);
			s.push((nodesref[2] + nodesref[3]) / 2);
			sa.push(area / 4);
			s.push((nodesref[0] + nodesref[3]) / 2);
			s.push((nodesref[1] + nodesref[3]) / 2);
			s.push((nodesref[2] + nodesref[3]) / 2);
			s.push(nodesref[3]);
			sa.push(area / 4);
			continue;
		}

		// We found an accurate solution to that part of the triangle!
		if (xInOmega)
		{
			for (int i = 0; i < integrals.size(); i++)
				integrals[i] += results2[i];
		}
		else
		{
			integrateRegion(QuadCell, 2, order() + 1, x, xref, vec1, vec2,
					vec1ref, vec2ref, noCurve, results1);
			for (int i = 0; i < integrals.size(); i++)
				integrals[i] += (results1[i] - results2[i]);
		}
	} // Stack depleted here! We're done.

	SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::getAdaptedQuadWeights Calculations in cell " << cellLID << " completed. " << nQuads << " quad(s) used.");

	Array<double> alphas;
	globalCurve.getIntegrationParams(alphas);
	for (int i = 0; i < quadWeights.size(); i++)
		quadWeights[i] = alphas[1] * integrals[i] + alphas[0]
				* (originalWeights[i] - integrals[i]);
	if (verb > 1 ) {
	    writeTable(Out::os(), tabs, quadWeights, quadWeights.size());
	}
	mesh.setSpecialWeight(2, cellLID, quadWeights);
}

void FeketeQuadrature::integrateRegion(const CellType& cellType,
		const int cellDim, const int innerOrder, const Point& x,
		const Point& xref, const Point& vec1, const Point& vec2,
		const Point& vec1ref, const Point& vec2ref,
		const ParametrizedCurve& curve, Array<double>& integrals) const
{
	int verb = 0;
	Tabs tabs;

	switch (cellType)
	{
	case TriangleCell:
	{
		// Prepare output array
		integrals.resize((order() + 1) * (order() + 2) / 2);
		for (int i = 0; i < integrals.size(); i++)
			integrals[i] = 0.0;

		// Get 'inner' integration points (along a line)
		int nrInnerQuadPts = innerOrder + 3;
		nrInnerQuadPts = (nrInnerQuadPts + nrInnerQuadPts % 2) / 2;
		GaussLobatto1D innerQuad(nrInnerQuadPts, 0.0, 1.0);

		Array<double> innerQuadPts(nrInnerQuadPts);
		Array<double> innerQuadWgts(nrInnerQuadPts);
		innerQuadPts = innerQuad.nodes();
		innerQuadWgts = innerQuad.weights();

		// Calculate intersection points of 'integration ray' with curve
		// Here one needs physical coordinates because of the curve
		for (int edgePos = 0; edgePos < nrInnerQuadPts; edgePos++)
		{
			// Direction of current integration
			Point ray = innerQuadPts[edgePos] * vec1 + (1.0
					- innerQuadPts[edgePos]) * vec2;
			Point rayref = innerQuadPts[edgePos] * vec1ref + (1.0
					- innerQuadPts[edgePos]) * vec2ref;

			int nrIntersect = 0;
			Array<double> t;
			Point rayEnd = x + ray;
			curve.returnIntersect(x, rayEnd, nrIntersect, t);

			double mu = 1.0; // w/o an intersection we integrate the whole line
			if (nrIntersect == 1)
			{
				mu = t[0];
			}
			else if (nrIntersect >= 2)
			{
				// More than one intersection -> We have to refine!
				integrals.resize(0);
				return;
			}

			// Array for values of basis functions at inner integration points
			Array<double> innerIntegrals(integrals.size());
			for (int i = 0; i < innerIntegrals.size(); i++)
				innerIntegrals[i] = 0.0;

			// We follow one of the rays (in reference coordinates)
			for (int rayPos = 0; rayPos < nrInnerQuadPts; rayPos++)
			{
				Point q = xref + mu * innerQuadPts[rayPos] * rayref;

				// Evaluate all basis functions at quadrature point
				Array<double> funcVals;
				evaluateAllBasisFunctions(TriangleCell, q, funcVals);

				// Do quadrature for all basis functions
				for (int i = 0; i < funcVals.size(); i++)
					innerIntegrals[i] += funcVals[i] * innerQuadWgts[rayPos]
							* innerQuadPts[rayPos];
			}
			for (int j = 0; j < integrals.size(); j++)
				integrals[j] += innerQuadWgts[edgePos] * mu * mu
						* innerIntegrals[j];
		}

		// We calculated on parts of the triangle -> Scale to right size
		double det = ::fabs(vec1ref[0] * vec2ref[1] - vec2ref[0] * vec1ref[1]);
		for (int j = 0; j < integrals.size(); j++)
			integrals[j] *= det;
	} break;

	case QuadCell:
	{
		// Prepare output array
		integrals.resize((order() + 1) * (order() + 1));
		for (int i = 0; i < integrals.size(); i++)
			integrals[i] = 0.0;
		SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::integrateRegion: Initializing quadrature on a QuadCell")

		// Get 'inner' integration points (along a line)
		int nrInnerQuadPts = innerOrder + 3;
		nrInnerQuadPts = (nrInnerQuadPts + nrInnerQuadPts % 2) / 2;
		GaussLobatto1D innerQuad(nrInnerQuadPts, 0.0, 1.0);

		Array<double> innerQuadPts(nrInnerQuadPts);
		Array<double> innerQuadWgts(nrInnerQuadPts);
		innerQuadPts = innerQuad.nodes();
		innerQuadWgts = innerQuad.weights();

		// We assume: no intersection along first axis
		for (int xPos = 0; xPos < nrInnerQuadPts; xPos++)
		{
			// Direction of current integration
			Point x_cur = x + innerQuadPts[xPos] * vec1;
			Point x_curref = xref + innerQuadPts[xPos] * vec1ref;

			SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::integrateRegion: Integrating at xPos " << xPos << " Ref point " << x_curref)

			// Intersection along second axis?
			int nrIntersect = 0;
			Array<double> t;
			Point x_curEnd = x_cur + vec2;
			curve.returnIntersect(x_cur, x_curEnd, nrIntersect, t);
			SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::integrateRegion: Nr. of intersections along y: " << t.size())

			double mu = 1.0; // w/o an intersection we integrate the whole line
			if (nrIntersect == 1)
			{
				mu = t[0];
				SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::integrateRegion: Intersection parameter is " << mu)
			}
			else if (nrIntersect >= 2)
			{
				// More than one intersection -> We have to refine!
				integrals.resize(0);
				SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::integrateRegion: More than one intersection -> Abort!")

				return;
			}

			// Array for values of basis functions at inner integration points
			Array<double> innerIntegrals(integrals.size());
			for (int jj = 0; jj < innerIntegrals.size(); jj++)
				innerIntegrals[jj] = 0.0;

			for (int yPos = 0; yPos < nrInnerQuadPts; yPos++)
			{
				Point y_curref = x_curref + mu * innerQuadPts[yPos] * vec2ref;
				SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::integrateRegion: Integrating at yPos " << yPos << " Ref point " << y_curref)

				// Evaluate all basis functions at quadrature point
				Array<double> funcVals;
				evaluateAllBasisFunctions(QuadCell, y_curref, funcVals);

				// Do quadrature for all basis functions
				for (int jj = 0; jj < funcVals.size(); jj++)
					innerIntegrals[jj] += funcVals[jj] * innerQuadWgts[yPos] * mu;
			}
			for (int ii = 0; ii < integrals.size(); ii++)
				integrals[ii] += innerQuadWgts[xPos] * innerIntegrals[ii];
		}

		// Scale to right size
		double det = ::fabs(vec1ref[0] * vec2ref[1] - vec2ref[0] * vec1ref[1]);
		for (int i = 0; i < integrals.size(); i++)
			integrals[i] *= det;
		SUNDANCE_MSG3(verb, tabs << "FeketeQuadrature::integrateRegion: Integration of region finished. Integral values " << integrals)
	} break;

	default:
#ifndef TRILINOS_7
		SUNDANCE_ERROR("FeketeQuadrature::integrateRegion() not available for cell type " << cellType)
		;
#else
		SUNDANCE_ERROR7("FeketeQuadrature::integrateRegion() not available for cell type " << cellType);
#endif
	}

}

void FeketeQuadrature::evaluateAllBasisFunctions(const CellType cellType,
		const Point& q, Array<double>& result) const
{

	switch (cellType)
	{
	case TriangleCell:
	{
		// Coefficients in _basisCoeffs together with PKD polynomials form
		// a Lagrange basis at Fekete points
		if (!_hasBasisCoeffs)
		{
			FeketeTriangleQuadrature::computeBasisCoeffs(order(), _basisCoeffs);
			_hasBasisCoeffs = true;
		}

    // KL added explicit cast of sqrt argument to double to allow
    // compilation on windows.
		// One has nFeketePts basis functions
		int nFeketePts = sqrt((double) _basisCoeffs.size());
		result.resize(nFeketePts);

		// Determine values of all PKD polynomials at evaluation point q
		Array<double> pkdVals(nFeketePts);
		FeketeTriangleQuadrature::evalPKDpolynomials(order(), q[0], q[1],
				&(pkdVals[0]));

		// For each Lagrange basis function, sum up weighted PKD polynomials
		// (The weights are found in _basisCoeffs)
		char trans = 'N';
		double alpha = 1.;
		double beta = 0.;
		int inc = 1;
		::dgemv_(&trans, &nFeketePts, &nFeketePts, &alpha, &(_basisCoeffs[0]),
				&nFeketePts, &(pkdVals[0]), &inc, &beta, &(result[0]), &inc);
		break;
	}
	case QuadCell:
	{
		// Coefficients in _basisCoeffs form a Lagrange basis at Fekete points
		if (!_hasBasisCoeffs)
		{
			FeketeQuadQuadrature::computeBasisCoeffs(order(), _basisCoeffs);
			_hasBasisCoeffs = true;
		}

    // KL added explicit cast of sqrt argument to double to allow
    // compilation on windows.
		// One has nFeketePts basis functions
		int nFeketePts = sqrt((double) _basisCoeffs.size());
		result.resize(nFeketePts);

		// Determine values of all polynomials at evaluation point q
		Array<double> polVals(nFeketePts);
		FeketeQuadQuadrature::evalPolynomials(nFeketePts, q[0], q[1],
				&(polVals[0]));

		// For each Lagrange basis function, sum up weighted polynomials
		// (The weights are found in _basisCoeffs)
		char trans = 'N';
		double alpha = 1.;
		double beta = 0.;
		int inc = 1;
		::dgemv_(&trans, &nFeketePts, &nFeketePts, &alpha, &(_basisCoeffs[0]),
				&nFeketePts, &(polVals[0]), &inc, &beta, &(result[0]), &inc);
		break;
	}
	default:
#ifndef TRILINOS_7
		SUNDANCE_ERROR("FeketeQuadrature::evaluateAllBasisFunctions() not available for cell type " << cellType)
		;
#else
		SUNDANCE_ERROR7("FeketeQuadrature::evaluateAllBasisFunctions() not available for cell type " << cellType);
#endif
	}
}
