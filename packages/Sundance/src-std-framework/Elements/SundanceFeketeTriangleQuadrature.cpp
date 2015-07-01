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


#include "SundanceFeketeTriangleQuadrature.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundancePoint.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;
using namespace Teuchos;

extern "C"
{
/* matrix-vector multiplication */
void dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda,
		double *x, int *incx, double *beta, double *y, int *incy);

/* LAPACK factorization */
void dgetrf_(const int* M, const int* N, double* A, const int* lda,
		const int* iPiv, int* info);

/* LAPACK inversion of factorized matrix */
void dgetri_(const int* n, double* a, const int* lda, const int* iPiv,
		double* work, const int* lwork, int* info);
}

/**
 *  Reference:
 *	T. Warburton, An explicit construction of interpolation nodes on the simplex
 *  J. Eng. Math. (2006) 56, pp. 247-262
 */
void FeketeTriangleQuadrature::getPoints(int order, Array<double>& wgt, Array<
		double>& x, Array<double>& y)
{
	Array<double> w;
	Array<int> multiplicity;
	Array<Array<double> > q;

	if (order == 1)
	{
		// Delete it? One gets 2nd order for the same price...
		multiplicity = tuple(3);
		q.resize(1);
		w = tuple(1.0 / 3.0);
		q[0] = tuple(1.0, 0.0, 0.0);
	}
	else if (order == 2)
	{
		// Corners have weights == 0, but we list them anyway for adaptive cell integration
		multiplicity = tuple(3, 3);
		q.resize(2);
		w = tuple(0.0, 1.0 / 3.0);
		q[0] = tuple(1.0, 0.0, 0.0);
		q[1] = tuple(0.0, 0.5, 0.5);
	}
	else if (order == 3)
	{
		multiplicity = tuple(1, 3, 6);
		q.resize(3);
		w = tuple(0.45, 1.0 / 60.0, 1.0 / 12.0);
		q[0] = tuple(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
		q[1] = tuple(1.0, 0.0, 0.0);
		q[2] = tuple(0.000000000000000, 0.276393202250021, 0.723606797749979);
	}
	else if (order == 4)
	{
		multiplicity = tuple(3, 3, 3, 6);
		q.resize(4);
		w = tuple(-0.002443433685378, 0.040684938686721, 0.200360939516545,
				0.047365444407723);
		q[0] = tuple(1.0, 0.0, 0.0);
		q[1] = tuple(0.0, 0.5, 0.5);
		q[2] = tuple(0.551583507555305, 0.224208246222347, 0.224208246222347);
		q[3] = tuple(0.000000000000000, 0.172673164646012, 0.827326835353988);
	}
	else if (order == 5)
	{
		multiplicity = tuple(3, 3, 3, 6, 6);
		q.resize(5);
		w = tuple(0.005695992854379, 0.125094845896272, 0.116460019007601,
				0.012270695896559, 0.030770541890982);
		q[0] = tuple(1.0, 0.0, 0.0);
		q[1] = tuple(0.684472514501908, 0.157763742749046, 0.157763742749046);
		q[2] = tuple(0.171245477332074, 0.414377261333963, 0.414377261333963);
		q[3] = tuple(0.000000000000000, 0.117472338035268, 0.882527661964732);
		q[4] = tuple(0.000000000000000, 0.357384241759678, 0.642615758240322);
	}
	else if (order == 6)
	{
		multiplicity = tuple(1, 3, 3, 3, 6, 6, 6);
		q.resize(7);
		w = tuple(0.078877036160935, -0.001814501354491, 0.021979435329947,
				0.055816419204784, 0.013312944349338, 0.013197985920130,
				0.089018887113590);
		q[0] = tuple(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
		q[1] = tuple(1.0, 0.0, 0.0);
		q[2] = tuple(0.0, 0.5, 0.5);
		q[3] = tuple(0.769520861253251, 0.115239569373374, 0.115239569373375);
		q[4] = tuple(0.000000000000000, 0.084888051860717, 0.915111948139283);
		q[5] = tuple(0.000000000000000, 0.265575603264643, 0.734424396735357);
		q[6] = tuple(0.127944523240232, 0.317617605129902, 0.554437871629866);
	}
	else if (order == 9)
	{
		multiplicity = tuple(1, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6);
		q.resize(12);
		w = tuple(0.054830037550851, 0.001827879484114, 0.019001036234018,
				0.036420799722133, 0.030554508665561, -0.000151089387317,
				0.005590374483628, 0.004964725734335, 0.006760863782205,
				0.020422905832759, 0.035099032364350, 0.040939402211986);
		q[0] = tuple(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
		q[1] = tuple(1.0, 0.0, 0.0);
		q[2] = tuple(0.888511356254118, 0.055744321872941, 0.055744321872941);
		q[3] = tuple(0.643273196352075, 0.178363401823962, 0.178363401823962);
		q[4] = tuple(0.067157296821421, 0.466421351589290, 0.466421351589290);
		q[5] = tuple(0.000000000000000, 0.040233045916771, 0.959766954083229);
		q[6] = tuple(0.000000000000000, 0.130613067447248, 0.869386932552752);
		q[7] = tuple(0.000000000000000, 0.261037525094778, 0.738962474905222);
		q[8] = tuple(0.000000000000000, 0.417360521166807, 0.582639478833193);
		q[9] = tuple(0.063288094133481, 0.162084689378845, 0.774627216487674);
		q[10] = tuple(0.066372381175690, 0.304692583561054, 0.628935035263256);
		q[11] = tuple(0.185067879630320, 0.326702931348863, 0.488229189020817);
	}
	/* Points and weights according to Fekete approach by
	 * Taylor, Wingate, Vincent 2000
	 * Unfortunately not sufficient digits for type 'double'!
	 else if (order == 6)
	 {
	 multiplicity = tuple(1, 3, 3, 3, 6, 6, 6);
	 q.resize(7);
	 w = tuple(0.2178563571, 0.1104193374, 0.0358939762, 0.0004021278,
	 0.1771348660, 0.0272344079, 0.0192969460);
	 q[0] = tuple(0.3333333333, 0.3333333333, 0.3333333334);
	 q[1] = tuple(0.7873290632, 0.1063354684, 0.1063354684);
	 q[2] = tuple(0.0000000000, 0.5000000000, 0.5000000000);
	 q[3] = tuple(1.0000000000, 0.0000000000, 0.0000000000);
	 q[4] = tuple(0.1171809171, 0.3162697959, 0.5665492870);
	 q[5] = tuple(0.0000000000, 0.2655651402, 0.7344348598);
	 q[6] = tuple(0.0000000000, 0.0848854223, 0.9151145777);
	 }
	 else if (order == 9)
	 {
	 multiplicity = tuple(1, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6);
	 q.resize(12);
	 w = tuple(0.1096011288, 0.0767491008, 0.0646677819, 0.0276211659,
	 0.0013925011, 0.0933486453, 0.0619010169, 0.0437466450,
	 0.0114553907, 0.0093115568, 0.0078421987, 0.0022457501);
	 q[0] = tuple(0.3333333333, 0.3333333333, 0.3333333334);
	 q[1] = tuple(0.6591363598, 0.1704318201, 0.1704318201);
	 q[2] = tuple(0.0600824712, 0.4699587644, 0.4699587644);
	 q[3] = tuple(0.9021308608, 0.0489345696, 0.0489345696);
	 q[4] = tuple(1.0000000000, 0.0000000000, 0.0000000000);
	 q[5] = tuple(0.1784337588, 0.3252434900, 0.4963227512);
	 q[6] = tuple(0.0588564879, 0.3010242110, 0.6401193011);
	 q[7] = tuple(0.0551758079, 0.1543901944, 0.7904339977);
	 q[8] = tuple(0.0000000000, 0.4173602935, 0.5826397065);
	 q[9] = tuple(0.0000000000, 0.2610371960, 0.7389628040);
	 q[10] = tuple(0.0000000000, 0.1306129092, 0.8693870908);
	 q[11] = tuple(0.0000000000, 0.0402330070, 0.9597669930);
	 }*/

	else
	{
#ifndef TRILINOS_7
		SUNDANCE_ERROR("symmetric Fekete quadrature rule order "
				<< order <<
				" for triangles not available");
#else
		SUNDANCE_ERROR7("symmetric Fekete quadrature rule order "
				<< order <<
				" for triangles not available");
#endif
	}

	for (int i = 0; i < q.length(); i++)
	{
		Array<Array<double> > qPerm;
		permute(multiplicity[i], q[i], qPerm);
		for (int j = 0; j < multiplicity[i]; j++)
		{
			x.append(qPerm[j][0]);
			y.append(qPerm[j][1]);
			wgt.append(w[i]);
		}
	}

}

bool FeketeTriangleQuadrature::supportsOrder(int order)
{
	if ((order >= 1 && order <= 6) || order == 9)
		return true;
	return false;
}

/**
 * Here we calculate coefficients for Proriol-Koornwinder-Dubiner polynomials
 * so that they form a Lagrange basis at given (Fekete quadrature) points in the
 * triangle
 */
void FeketeTriangleQuadrature::computeBasisCoeffs(const int order, Array<double>& basisCoeffs)
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

	// We construct a Lagrange basis at feketePts (there are nFeketePts of them)
	// Each basis polynomial itself is given by a linear combination of
	// nFeketePts PKD polynomials and their coefficients
	basisCoeffs.resize(nFeketePts * nFeketePts);

	// Let's compute the coefficients:
	// Build Vandermonde matrix of PKD basis at Fekete points
	for (int n = 0; n < nFeketePts; n++)
	{
		// Set pointer to beginning of n-th row and determine values of
		// PKD polynomials at n-th Fekete point
		double* start = &(basisCoeffs[n * nFeketePts]);
		evalPKDpolynomials(order, feketePts[n][0], feketePts[n][1], start);
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
			"FeketeTriangleQuadrature::computeBasisCoeffs(): factorization of generalized Vandermonde matrix failed");

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
			"FeketeTriangleQuadrature::computeBasisCoeffs(): inversion of generalized Vandermonde matrix failed");
}

/**
 * Evaluates all basis functions of a Proriol-Koornwinder-Dubiner basis
 * up to the given order at (x,y) in reference (barycentric) coordinates of a triangle;
 * Missing third coordinate z = 1-x-y
 */
void FeketeTriangleQuadrature::evalPKDpolynomials(int order, double x,
		double y, double* resultPtr)
{
	int i = 0;
	double Legendre = 1.0;
	double Legendre1 = 0.0;
	for (int k = 0; k <= order; k++)
	{
		double Jacobi = 1.0;
		double Jacobi1 = 0.0;

		// 'i' implements a bijection polynomial index (k,l) onto basis index (i)
		for (int l = 0; l <= order - k; l++)
		{
			// Overall we will write (order+1)*(order+2)/2 values
			resultPtr[i++] = Legendre * pow(x + y, k) * Jacobi;

			// Update Jacobi polynomial
			double Jacobi2 = Jacobi1;
			Jacobi1 = Jacobi;
			int c1 = 2 * (l + 1) * (l + 2 * k + 2) * (2 * l + 2 * k + 1);
			int c2 = (2 * l + 2 * k + 2) * (2 * k + 1) * (2 * k + 1);
			int c3 = (2 * l + 2 * k + 1) * (2 * l + 2 * k + 2) * (2 * l + 2 * k
					+ 3);
			int c4 = 2 * (l + 2 * k + 1) * l * (2 * l + 2 * k + 3);
			Jacobi = ((c2 + c3 * (1.0 - 2.0 * x - 2.0 * y)) * Jacobi1 - c4
					* Jacobi2) / c1;
		}

		// Update Legendre polynomial
		double Legendre2 = Legendre1;
		Legendre1 = Legendre;
		double xy = x + y;
		if (::fabs(xy) > 0.0)
		{
			xy = (2.0 * k + 1.0) * ((x - y) / xy) * Legendre1;
		}
		Legendre = (xy - k * Legendre2) / (k + 1);
	}
}

void FeketeTriangleQuadrature::permute(int m, const Array<double>& q, Array<
		Array<double> >& qPerm)
{
	qPerm.resize(m);
	if (m == 1)
	{
		qPerm[0] = q;
	}
	else if (m == 3)
	{
		qPerm[0] = tuple(q[0], q[1], q[2]);
		qPerm[1] = tuple(q[1], q[0], q[2]);
		qPerm[2] = tuple(q[2], q[1], q[0]);
	}
	else if (m == 6)
	{
		qPerm[0] = tuple(q[0], q[1], q[2]);
		qPerm[1] = tuple(q[0], q[2], q[1]);
		qPerm[2] = tuple(q[1], q[0], q[2]);
		qPerm[3] = tuple(q[1], q[2], q[0]);
		qPerm[4] = tuple(q[2], q[1], q[0]);
		qPerm[5] = tuple(q[2], q[0], q[1]);
	}
	else
	{
#ifndef TRILINOS_7
		SUNDANCE_ERROR("invalid multiplicity "
				<< m <<
				" in FeketeTriangleQuadrature::permute()");
#else
		SUNDANCE_ERROR7("invalid multiplicity "
				<< m <<
				" in FeketeTriangleQuadrature::permute()");
#endif
	}
}

bool FeketeTriangleQuadrature::test(int p)
{
	Array<double> w;
	Array<double> x;
	Array<double> y;

	getPoints(p, w, x, y);
	bool pass = true;

	for (int a = 0; a <= p; a++)
	{
		for (int b = 0; b < p - a; b++)
		{
			int cMax = p - a - b;
			for (int c = 0; c <= cMax; c++)
			{
				double sum = 0.0;
				for (int q = 0; q < w.length(); q++)
				{
					sum += 0.5 * w[q] * pow(x[q], (double) a) * pow(y[q],
							(double) b) * pow(1.0 - x[q] - y[q], (double) c);
				}
				double err = fabs(sum - exact(a, b, c));
				bool localPass = err < 1.0e-14;
				pass = pass && localPass;
				if (!localPass)
				{
					fprintf(
							stderr,
							"order=%d m (%d, %d, %d) q=%22.15g exact=%22.15g\n",
							p, a, b, c, sum, exact(a, b, c));
					std::cerr << "error = " << err << std::endl;
				}
			}
		}
	}
	return pass;
}

double FeketeTriangleQuadrature::exact(int a, int b, int c)
{
	return fact(a) * fact(b) * fact(c) / fact(a + b + c + 2);
}

double FeketeTriangleQuadrature::fact(int x)
{
	if (x == 0)
		return 1.0;
	return x * fact(x - 1);
}

