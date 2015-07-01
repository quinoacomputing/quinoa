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

#include "SundancePointwiseUserDefFunctor.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


PointwiseUserDefFunctor0::PointwiseUserDefFunctor0(const std::string& name, 
                                                   int domainDim, 
                                                   int rangeDim)
  :  UserDefFunctor(name, domainDim, rangeDim)
{}

void PointwiseUserDefFunctor0::evaluationCallback(int nPoints, int maxDiffOrder,
                                    const double** in,
                                    double** out) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(maxDiffOrder > 0, std::runtime_error,
                     "diff order = " << maxDiffOrder 
                     << " not supported for functor "
                     << name());

  static Array<double> x;
  x.resize(domainDim());
  
  static Array<double> f;
  f.resize(rangeDim());

  for (int i=0; i<nPoints; i++)
    {
      for (int j=0; j<domainDim(); j++) x[j] = in[j][i];
      const double* xp = &(x[0]);
      double* fp = &(f[0]);
      eval0(xp, fp);
      for (int j=0; j<rangeDim(); j++) out[j][i] = fp[j];
    }
}

PointwiseUserDefFunctor1::PointwiseUserDefFunctor1(const std::string& name, 
                                                   int domainDim, 
                                                   int rangeDim)
  : PointwiseUserDefFunctor0(name, domainDim, rangeDim)
{}

void PointwiseUserDefFunctor1::evaluationCallback(int nPoints, int maxDiffOrder,
                                    const double** in,
                                    double** out) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(maxDiffOrder > 1, std::runtime_error,
                     "diff order = " << maxDiffOrder 
                     << " not supported for functor "
                     << name());

  static Array<double> x;
  x.resize(domainDim());
  
  static Array<double> f;
  if (maxDiffOrder==1) f.resize(rangeDim() * (1 + domainDim()) );
  else f.resize(rangeDim());
  double* fp = &(f[0]);

  for (int i=0; i<nPoints; i++)
    {
      for (int j=0; j<domainDim(); j++) x[j] = in[j][i];
      const double* xp = &(x[0]);

      if (maxDiffOrder==1) 
        {
          double* dfp = &(f[rangeDim()]);
          eval1(xp, fp, dfp);
        }
      else eval0(xp, fp);
      for (int j=0; j<f.size(); j++) out[j][i] = fp[j];
    }
}


void PointwiseUserDefFunctor1::eval0(const double* in, double* out) const 
{
  static Array<double> dummy;
  dummy.resize(domainDim() * rangeDim());

  eval1(in, out, &(dummy[0]));
}



PointwiseUserDefFunctor2::PointwiseUserDefFunctor2(const std::string& name, 
                                                   int domainDim, 
                                                   int rangeDim)
  : PointwiseUserDefFunctor1(name, domainDim, rangeDim)
{}

void PointwiseUserDefFunctor2::evaluationCallback(int nPoints, int maxDiffOrder,
                                                  const double** in,
                                                  double** out) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(maxDiffOrder > 2 || maxDiffOrder < 0, std::runtime_error,
                     "diff order = " << maxDiffOrder 
                     << " not supported for functor "
                     << name());

  int nTotal = 1;
  int numFirst = domainDim();
  int numSecond = domainDim()*(domainDim()+1)/2;

  static Array<double> x;
  x.resize(domainDim());
  
  static Array<double> f;

  if (maxDiffOrder > 0) nTotal += numFirst;
  if (maxDiffOrder > 1) nTotal += numSecond;

  f.resize(rangeDim() * nTotal);

  double* fp = &(f[0]);

  for (int i=0; i<nPoints; i++)
    {
      for (int j=0; j<domainDim(); j++) x[j] = in[j][i];
      const double* xp = &(x[0]);

      if (maxDiffOrder==0)
        {
          eval0(xp, fp);
        }
      else if (maxDiffOrder==1)
        {
          double* dfp = &(f[rangeDim()]);
          eval1(xp, fp, dfp);
        }
      else if (maxDiffOrder==2)
        {
          double* dfp = &(f[rangeDim()]);
          double* d2fp = &(f[rangeDim()*(1 + domainDim())]);
          eval2(xp, fp, dfp, d2fp);
        }
      else
        {
          TEUCHOS_TEST_FOR_EXCEPT(true);
        }
      for (int j=0; j<f.size(); j++) out[j][i] = fp[j];
    }
}


void PointwiseUserDefFunctor2::eval0(const double* in, double* f) const 
{
  static Array<double> dummy1;
  static Array<double> dummy2;
  dummy1.resize(rangeDim() *  domainDim());
  dummy2.resize(rangeDim() *  domainDim()*(domainDim()+1)/2);

  eval2(in, f, &(dummy1[0]), &(dummy2[0]));
}


void PointwiseUserDefFunctor2::eval1(const double* in, double* f, double* df) const 
{
  static Array<double> dummy2;
  dummy2.resize(rangeDim() *  domainDim()*(domainDim()+1)/2);

  eval2(in, f, df, &(dummy2[0]));
}
