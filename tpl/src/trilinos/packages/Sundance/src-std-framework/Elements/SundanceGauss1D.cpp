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


#include "SundanceGauss1D.hpp"
#ifdef _MSC_VER
# include "winmath.h"
#endif

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

Gauss1D::Gauss1D(int n)
	: nodes_(n), weights_(n)
{
	computeWeights(n, -1.0, 1.0);
}

Gauss1D::Gauss1D(int n, double a, double b)
	: nodes_(n), weights_(n)
{
	computeWeights(n, a, b);
}



void Gauss1D::computeWeights(int n, double a, double b)
{
	int m = (n+1)/2;

	double xMid = (b+a)/2.0;
	double halfWidth = (b-a)/2.0;
	
	for (int i=0; i<m; i++)
		{
			// initial guess
			double z = cos(M_PI*(i+0.75)/(n+0.5));
			double dP;
			double zOld;
			double tol = 1.0e-14;
			// newton's method
			do
				{
					double p1 = 1.0;
					double p2 = 0.0;
					for (int j=1; j<=n; j++)
						{
							double p3 = p2;
							p2 = p1;
							p1 = ((2.0*j-1.0)*z*p2 - (j-1.0)*p3)/j;
						}
					dP = n*(z*p1-p2)/(z*z-1.0);
					zOld =z;
					z = zOld - p1/dP;
				}
			while ( fabs(z-zOld) > tol );
			nodes_[i] = xMid - halfWidth*z;
			nodes_[n-i-1] = xMid + halfWidth*z;
			weights_[i] = 2.0*halfWidth/((1.0-z*z)*dP*dP);
			weights_[n-i-1] = weights_[i];
		}
}
					
	
bool Gauss1D::unitTest()
{
	std::cerr << "------------------ Gauss1D unit test ----------------------" 
			 << std::endl;

	Gauss1D q(20, 0.0, M_PI);

	double sum = 0.0;
	for (int i=0; i<q.nPoints(); i++)
		{
			sum += q.weights()[i]*sin(q.nodes()[i]);
		}
	std::cerr << "integral of sin(x) over [0, pi] = " << sum << std::endl;
	double sumErr = fabs(sum - 2.0);
	bool sumPass = sumErr < 1.0e-10;
	std::cerr << "error = " << sumErr << std::endl;
	if (sumPass) std::cerr << "Gauss1D sine test PASSED" << std::endl;
	else std::cerr << "Gauss1D sine test FAILED" << std::endl;
	std::cerr << "------------------ End Gauss1D unit test ----------------------" 
			 << std::endl;
	return sumPass;
}


