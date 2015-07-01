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


#include "SundanceTetQuadrature.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


void TetQuadrature::getPoints(int order, Array<double>& wgt,
                              Array<double>& x,
                              Array<double>& y,
                              Array<double>& z)
{
	Array<double> w;
	Array<int> multiplicity;
	Array<Array<double> > q;

	if (order==1)
		{
			multiplicity = tuple(1);
			w = tuple(1.0);
			q.resize(1);
			q[0] = tuple(0.25);
		}
	else if (order==2)
		{
			multiplicity = tuple(4);
			w = tuple(0.25);
			q.resize(1);
			q[0] = tuple(0.5854101966249685, 0.1381966011250105);
		}
	else if (order==4)
		{
			multiplicity = tuple(4, 12);
			w = tuple(0.05037379410012282, 0.06654206863329239);
			q.resize(2);
			q[0] = tuple(0.7716429020672371, 0.7611903264425430e-01);
			q[1] = tuple(0.1197005277978019, 0.7183164526766925e-01, 0.4042339134672644);
		}
	else if (order==6)
		{
      multiplicity = tuple(1, 4, 12, 12);
      w = tuple(0.9040129046014750e-01, 0.1911983427899124e-01,
                0.4361493840666568e-01, 0.2581167596199161e-01);
      q.resize(4);
      q[0] = tuple(0.25);
      q[1] = tuple(0.8277192480479295, 0.5742691731735683e-01);
      q[2] = tuple(0.5135188412556341e-01, 0.4860510285706072, 0.2312985436519147);
      q[3] = tuple(0.2967538129690260, 0.6081079894015281, 0.4756909881472290e-01);
    }
	else
		{
#ifndef TRILINOS_7
			SUNDANCE_ERROR("symmetric quadrature rule order " 
                     << order << 
                     " not available for triangles");
#else
			SUNDANCE_ERROR7("symmetric quadrature rule order " 
                     << order << 
                     " not available for triangles");
#endif
		}

	for (int i=0; i<q.length(); i++)
		{
			Array<Array<double> > qPerm;
			permute(multiplicity[i], q[i], qPerm);
			for (int j=0; j<multiplicity[i]; j++)
				{
					x.append(qPerm[j][0]);
					y.append(qPerm[j][1]);
					z.append(qPerm[j][2]);
					wgt.append(w[i]);
				}
		}
}

bool TetQuadrature::supportsOrder(int order)
{
  if (order==1 || order==2 || order==4 || order==6) return true;
  return false;
}

void TetQuadrature::permute(int m, const Array<double>& q,
                            Array<Array<double> >& qPerm)
{
	qPerm.resize(m);
	if (m==1)
		{
			qPerm[0] = tuple(q[0], q[0], q[0], q[0]);
		}
	else if (m==4)
		{
			qPerm[0] = tuple(q[0], q[1], q[1], q[1]);
			qPerm[1] = tuple(q[1], q[0], q[1], q[1]);
			qPerm[2] = tuple(q[1], q[1], q[0], q[1]);
			qPerm[3] = tuple(q[1], q[1], q[1], q[0]);
		}
	else if (m==12)
		{
      qPerm[0] = tuple(q[0], q[1], q[2], q[2]);
      qPerm[1] = tuple(q[0], q[2], q[1], q[2]);
      qPerm[2] = tuple(q[0], q[2], q[2], q[1]);
      qPerm[3] = tuple(q[1], q[0], q[2], q[2]);
      qPerm[4] = tuple(q[2], q[0], q[1], q[2]);
      qPerm[5] = tuple(q[2], q[0], q[2], q[1]);
      qPerm[6] = tuple(q[1], q[2], q[0], q[2]);
      qPerm[7] = tuple(q[2], q[1], q[0], q[2]);
      qPerm[8] = tuple(q[2], q[2], q[0], q[1]);
      qPerm[9] = tuple(q[1], q[2], q[2], q[0]);
      qPerm[10] = tuple(q[2], q[1], q[2], q[0]);
      qPerm[11] = tuple(q[2], q[2], q[1], q[0]);
		}
	else
		{
#ifndef TRILINOS_7
			SUNDANCE_ERROR("invalid multiplicity " 
                     << m <<
                     " in TetQuadrature::permute()");
#else
			SUNDANCE_ERROR7("invalid multiplicity " 
                     << m <<
                     " in TetQuadrature::permute()");
#endif
		}
}

bool TetQuadrature::test(int p)
{
	Array<double> w;
	Array<double> x;
	Array<double> y;
	Array<double> z;

	getPoints(p, w, x, y, z);
	bool pass = true;
	
	for (int a=0; a<=p; a++)
		{
			for (int b=0; b<p-a; b++)
				{
					int cMax = p - a - b;
					for (int c=0; c<=cMax; c++)
						{
              int dMax = p - a - b - c;
              for (int d=0; d<=dMax; d++)
                {
                  double sum = 0.0;
                  for (int q=0; q<w.length(); q++)
                    {
                      sum += (1.0/6.0)*w[q] * pow(x[q], (double) a) * pow(y[q], (double) b)
                        * pow(z[q], (double) c)
                        * pow(1.0 - x[q] - y[q] - z[q], (double) d);
                    }
                  double err = fabs(sum - exact(a,b,c,d));
                  bool localPass = err < 1.0e-14;
                  pass = pass && localPass;
                  if (!localPass)
                    {
                      fprintf(stderr, "order=%d m (%d, %d, %d %d) q=%22.15g exact=%22.15g\n", p, a, b, c, d, sum, exact(a, b, c, d));
                      std::cerr << "error = " << err << std::endl;
                    }
                }
            }
        }
    }
	return pass;
}

double TetQuadrature::exact(int a, int b, int c, int d)
{
	return fact(a)*fact(b)*fact(c)*fact(d)/fact(a+b+c+d+3);
}

double TetQuadrature::fact(int x)
{
	if (x==0) return 1.0;
	return x*fact(x-1);
}

