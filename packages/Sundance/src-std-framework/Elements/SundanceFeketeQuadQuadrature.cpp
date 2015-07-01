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


#include "SundanceFeketeQuadQuadrature.hpp"
#include "SundanceOut.hpp"
#include "SundancePoint.hpp"
#include "SundanceGaussLobatto1D.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;
using namespace Teuchos;

extern "C"
{

/* LAPACK factorization */
void dgetrf_(const int* M, const int* N, double* A, const int* lda,
		const int* iPiv, int* info);

/* LAPACK inversion of factorized matrix */
void dgetri_(const int* n, double* a, const int* lda, const int* iPiv,
		double* work, const int* lwork, int* info);
}

void FeketeQuadQuadrature::getPoints(int order, Array<double>& wgt, Array<
		double>& x, Array<double>& y)
{
	int p = order + 3;
	p = p + (p % 2);
	int nNodes = p / 2;
	GaussLobatto1D rule(nNodes, 0.0, 1.0);
	Array<double> s = rule.nodes();
	Array<double> t = s;
	Array<double> w = rule.weights();
	int n = rule.nPoints();

	wgt.resize(n * n);
	x.resize(n * n);
	y.resize(n * n);

	int k = 0;
	for (int i = 0; i < n; i++)
	{
		double p = s[i];
		for (int j = 0; j < n; j++, k++)
		{
			double q = t[j]; //similar to the p value, caz we have quad
			x[k] = p;
			y[k] = q;
			wgt[k] = w[i] * w[j];
		}
	}
}

void FeketeQuadQuadrature::computeBasisCoeffs(const int order, Array<double>& basisCoeffs)
{
	// Get Fekete points of chosen order
	Array<Point> feketePts;

	Array<double> x;
	Array<double> y;
	Array<double> w;
	getPoints(order, w, x, y);
	int nFeketePts = w.length();
	feketePts.resize(nFeketePts);
	for (int i = 0; i < nFeketePts; i++)
		feketePts[i] = Point(x[i], y[i]);

	// We construct a Lagrange basis at feketePts (there are nFeketePts basis functions)
	basisCoeffs.resize(nFeketePts * nFeketePts);

	// Let's compute the coefficients:
	// Build Vandermonde matrix at Fekete points
	for (int n = 0; n < nFeketePts; n++)
	{
		// Set pointer to beginning of n-th row and determine values of
		// polynomials at n-th Fekete point
		double* start = &(basisCoeffs[n * nFeketePts]);
		evalPolynomials(nFeketePts, feketePts[n][0], feketePts[n][1], start);
	}

	// Invert Vandermonde matrix to obtain basis coefficients:
	// LAPACK error flag, array for switched rows, work array
	int lapack_err = 0;
	Array<int> pivot;
	pivot.resize(nFeketePts);
	Array<double> work;
	work.resize(1);
	int lwork = -1;

	// LU factorization
	::dgetrf_(&nFeketePts, &nFeketePts, &(basisCoeffs[0]), &nFeketePts,
			&(pivot[0]), &lapack_err);

	TEUCHOS_TEST_FOR_EXCEPTION(
			lapack_err != 0,
			std::runtime_error,
			"FeketeQuadQuadrature::computeBasisCoeffs(): factorization of generalized Vandermonde matrix failed");

	// Determine work array size and invert factorized matrix
	::dgetri_(&nFeketePts, &(basisCoeffs[0]), &nFeketePts, &(pivot[0]),
			&(work[0]), &lwork, &lapack_err);
	lwork = (int) work[0];
	work.resize(lwork);
	::dgetri_(&nFeketePts, &(basisCoeffs[0]), &nFeketePts, &(pivot[0]),
			&(work[0]), &lwork, &lapack_err);

	TEUCHOS_TEST_FOR_EXCEPTION(
			lapack_err != 0,
			std::runtime_error,
			"FeketeQuadQuadrature::computeBasisCoeffs(): inversion of generalized Vandermonde matrix failed");
}

void FeketeQuadQuadrature::evalPolynomials(int nPts, double x,
		double y, double* resultPtr)
{
	// Calculate all (normalized and shifted) Legendre polynomials < order
	int i = 0;
	int order = (int) sqrt((double) nPts);

	double xLegendre  = 1.;
	double xLegendre1 = 0.;

	for (int k = 0; k < order; k++)
	{
		double yLegendre  = 1.;
		double yLegendre1 = 0.;
		double normxLegendre = sqrt( 2. * k + 1 ) * xLegendre;

		// 'i' implements a bijection polynomial index (k,l) onto basis index (i)
		for (int l = 0; l < order ; l++)
		{
			double normyLegendre = sqrt( 2. * l + 1 ) * yLegendre;

			// We will write order*order values
			resultPtr[i++] = normxLegendre * normyLegendre;

			double yLegendre2 = yLegendre1;
			yLegendre1 = yLegendre;
			yLegendre = ( (2 * l + 1) * (2 * y - 1) * yLegendre1 - l * yLegendre2 ) / ( l + 1 );
		}

		double xLegendre2 = xLegendre1;
		xLegendre1 = xLegendre;
		xLegendre = ( (2 * k + 1) * (2 * x - 1) * xLegendre1 - k * xLegendre2 ) / ( k + 1 );
	}
}

bool FeketeQuadQuadrature::test(int p)
{
	Array<double> w;
	Array<double> x;
	Array<double> y;

	getPoints(p, w, x, y);
	bool pass = true;
	for (int a = 0; a <= p; a++)
	{
		int bMax = p - a;
		for (int b = 0; b <= bMax; b++)
		{
			double sum = 0.0;
			for (int q = 0; q < w.length(); q++)
			{
				sum += w[q] * pow(x[q], (double) a) * pow(y[q], (double) b);
			}
			double err = fabs(sum - exact(a, b));
			bool localPass = err < 1.0e-14;
			pass = pass && localPass;
			if (!localPass)
			{
				fprintf(stderr,
						"order=%d m (%d, %d) q=%22.15g exact=%22.15g\n", p, a,
						b, sum, exact(a, b));
				std::cerr << "error = " << err << std::endl;
			}
		}
	}
	return pass;
}

double FeketeQuadQuadrature::exact(int a, int b)
{
	return 1.0 / (a + 1) / (b + 1);
}

