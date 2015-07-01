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


#include "SundanceTriangleQuadrature.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceGauss1D.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

void TriangleQuadrature::getPoints(int order, Array<double>& wgt,
                                   Array<double>& x,
                                   Array<double>& y)
{
  if (!getSymmetricPoints(order, wgt, x, y)) 
    {
      getNonsymmetricPoints(order, wgt, x, y);
    }
}


bool TriangleQuadrature::getSymmetricPoints(int order, Array<double>& wgt,
                                            Array<double>& x,
                                            Array<double>& y)
{

	int np;
	Array<double> w;
	Array<int> multiplicity;
	Array<Array<double> > q;

	if (order==1)
		{
			multiplicity = tuple(1);
			np = 1;
			w = tuple(1.0);
			q.resize(1);
			q[0] = tuple(1.0/3.0, 1.0/3.0, 1.0/3.0);
		}
	else if (order==2)
		{
			multiplicity = tuple(3);
			np = 3;
			w = tuple(1.0/3.0);
			q.resize(1);
			q[0] = tuple(2.0/3.0, 1.0/6.0, 1.0/6.0);
		}
	else if (order==3)
		{
			multiplicity = tuple(6);
			np = 6;
			w = tuple(1.0/6.0);
			q.resize(1);
			q[0] = tuple(0.659027622374092, 0.231933368553031, 0.109039009072877);
		}
	else if (order==4)
		{
			multiplicity = tuple(3, 3);
			np = 6;
			w = tuple(0.109951743655322, 0.223381589678011);
			q.resize(2);
			q[0] = tuple(0.816847572980459, 0.091576213509771, 0.091576213509771);
			q[1] = tuple(0.108103018168070, 0.445948490915965, 0.445948490915965);
		}
	else if (order==5)
		{
			multiplicity = tuple(1, 3, 3);
			np = 7;
			q.resize(3);
			w = tuple(0.22500000000000, 0.125939180544827, 0.132394152788506);
			q[0] = tuple(1.0/3.0, 1.0/3.0, 1.0/3.0);
			q[1] = tuple(0.797426985353087, 0.101286507323456, 0.101286507323456);
			q[2] = tuple(0.059715871789770, 0.470142064105115, 0.470142064105115);
		}
	else if (order==6)
		{
			multiplicity = tuple(3, 3, 6);
			np = 12;
			q.resize(3);
			w = tuple(0.050844906370207, 0.116786275726379, 0.082851075618374);
			q[0] = tuple(0.873821971016996, 0.063089014491502, 0.063089014491502);
			q[1] = tuple(0.501426509658179, 0.249286745170910, 0.249286745170910);
			q[2] = tuple(0.636502499121399, 0.310352451033784, 0.053145049844817);
		}
	else
		{
      return false;
		}

	for (int i=0; i<q.length(); i++)
		{
			Array<Array<double> > qPerm;
			permute(multiplicity[i], q[i], qPerm);
			for (int j=0; j<multiplicity[i]; j++)
				{
					x.append(qPerm[j][0]);
					y.append(qPerm[j][1]);
					wgt.append(w[i]);
				}
		}
	
	return true;
}




void TriangleQuadrature::getNonsymmetricPoints(int order, Array<double>& wgt,
                                               Array<double>& x,
                                               Array<double>& y)
{
  int nNodes = (order+3)/2;
  Gauss1D rule(nNodes, -1.0, 1.0);
  Array<double> s = rule.nodes();
  Array<double> t = s;
  Array<double> w = rule.weights();
  int n = rule.nPoints();

  wgt.resize(n*n);
  x.resize(n*n);
  y.resize(n*n);

  int k=0;
  for (int i=0; i<n; i++)
    {
      double p = (1.0+s[i])/2.0;
      double J = 1.0-p;
      for (int j=0; j<n; j++, k++)
        {
          double q = (1.0 - p)*(1.0+t[j])/2.0;
          x[k] = p;
          y[k] = q;
          wgt[k] = 0.5*w[i]*w[j]*J;
        }
    }
}


void TriangleQuadrature::permute(int m, const Array<double>& q,
																 Array<Array<double> >& qPerm)
{
	qPerm.resize(m);
	if (m==1)
		{
			qPerm[0] = q;
		}
	else if (m==3)
		{
			qPerm[0] = tuple(q[0], q[1], q[2]);
			qPerm[1] = tuple(q[1], q[0], q[2]);
			qPerm[2] = tuple(q[2], q[1], q[0]);
		}
	else if (m==6)
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
                     " in TriangleQuadrature::permute()");
#else
			SUNDANCE_ERROR7("invalid multiplicity " 
                     << m <<
                     " in TriangleQuadrature::permute()");
#endif
		}
}

bool TriangleQuadrature::test(int p)
{
	Array<double> w;
	Array<double> x;
	Array<double> y;

	getPoints(p, w, x, y);
	bool pass = true;
	
	for (int a=0; a<=p; a++)
		{
			for (int b=0; b<p-a; b++)
				{
					int cMax = p - a - b;
					for (int c=0; c<=cMax; c++)
						{
							double sum = 0.0;
							for (int q=0; q<w.length(); q++)
								{
									sum += 0.5*w[q] * pow(x[q], (double) a) * pow(y[q], (double) b)
										* pow(1.0 - x[q] - y[q], (double) c);
								}
              double err = fabs(sum - exact(a,b,c));
              bool localPass = err < 1.0e-14;
							pass = pass && localPass;
              if (!localPass)
                {
                  fprintf(stderr, "order=%d m (%d, %d, %d) q=%22.15g exact=%22.15g\n", p, a, b, c, sum, exact(a, b, c));
                  std::cerr << "error = " << err << std::endl;
                }
						}
				}
		}
	return pass;
}

double TriangleQuadrature::exact(int a, int b, int c)
{
	return fact(a)*fact(b)*fact(c)/fact(a+b+c+2);
}

double TriangleQuadrature::fact(int x)
{
	if (x==0) return 1.0;
	return x*fact(x-1);
}


