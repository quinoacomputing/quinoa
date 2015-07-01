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


#include "SundanceBrickQuadrature.hpp"
#include "SundanceOut.hpp"
#include "SundanceGauss1D.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;
using namespace Teuchos;

void BrickQuadrature::getPoints(int order, Array<double>& wgt,
		Array<double>& x, Array<double>& y, Array<double>& z)
{
	int p = order + 1;
	p = p + (p % 2);
	int nNodes = p / 2;
	Gauss1D rule(nNodes, 0.0, 1.0);
	Array<double> d1 = rule.nodes();
	Array<double> d2 = d1;
	Array<double> d3 = d1;
	Array<double> w = rule.weights();
	int n = rule.nPoints();

	wgt.resize(n * n * n);
	x.resize(n * n * n);
	y.resize(n * n * n);
	z.resize(n * n * n);

	int k = 0;
	for (int i = 0; i < n; i++)
	{
		double p = d1[i];
		for (int j = 0; j < n; j++)
		{
			double q = d2[j]; //similar to the p value, caz we have quad
			for (int l = 0; l < n; l++, k++)
			{
				double r = d3[l];
				x[k] = p;
				y[k] = q;
				z[k] = r;
				wgt[k] = w[i] * w[j] * w[l];
			}
		}
	}
}

bool BrickQuadrature::test(int p)
{
	Array<double> w;
	Array<double> x;
	Array<double> y;
	Array<double> z;

	getPoints(p, w, x, y, z);
	bool pass = true;
	for (int a = 0; a <= p; a++)
	{
		int bMax = p - a;
		for (int b = 0; b <= bMax; b++)
		{
			int cMax = bMax - b;
			for (int c = 0; c <= cMax; c++)
			{
				double sum = 0.0;
				for (int q = 0; q < w.length(); q++)
				{
					sum += w[q] * pow(x[q], (double) a) * pow(y[q], (double) b)
							* pow(z[q], (double) c);
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

double BrickQuadrature::exact(int a, int b, int c)
{
	return 1.0 / (a + 1) / (b + 1) / (c + 1);
}

