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


#include "SundanceQuadQuadrature.hpp"
#include "SundanceOut.hpp"
#include "SundanceGauss1D.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;
using namespace Teuchos;

void QuadQuadrature::getPoints(int order, Array<double>& wgt, Array<double>& x,
		Array<double>& y)
{
	int p = order + 1;
	p = p + (p % 2);
	int nNodes = p / 2;
	Gauss1D rule(nNodes, 0.0, 1.0);
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

bool QuadQuadrature::test(int p)
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

double QuadQuadrature::exact(int a, int b)
{
	return 1.0 / (a + 1) / (b + 1);
}

