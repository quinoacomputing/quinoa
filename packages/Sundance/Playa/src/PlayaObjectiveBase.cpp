/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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


#include "PlayaObjectiveBase.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaGhostView.hpp"
#include "PlayaLinearCombinationImpl.hpp"



namespace Playa
{
using std::endl;
using std::ostream;
using std::setw;

bool ObjectiveBase::fdCheck(const Vector<double>& xIn,
  double tol,
  int verbosity) const 
{
  double f0, fPlus, fMinus;
  /* get a FD step appropriate to this problem */
  double h = fdStep();

  PLAYA_MSG1(verbosity, "running fdCheck with stepsize=" << h);

  Vector<double> x = xIn.copy();
  Vector<double> x0 = x.copy();
  Vector<double> gf = x.copy();
  evalGrad(xIn, f0, gf);

  int nTot = x.space().dim();
  int n = x.space().numLocalElements();
  int lowestIndex = x.space().baseGlobalNaturalIndex();

  Array<double> df_dx(n);

  for (int globalIndex=0; globalIndex<nTot; globalIndex++)
  {
    double tmp=0.0;
    bool isLocal = globalIndex >= lowestIndex 
      && globalIndex < (lowestIndex+n);
    int localIndex = globalIndex - lowestIndex;
    /* set point to xPlus */
    if (isLocal)
    {
      tmp = x[localIndex];
      x[localIndex] = tmp + h;
    }

    /** eval at xPlus*/
    eval(x, fPlus);

    /* set point to xMinus */
    if (isLocal)
    {
      x[localIndex] = tmp - h;
    }

    /* eval at xMinus */
    eval(x, fMinus);
      
    /* check error */
    if (isLocal)
    {
      df_dx[localIndex] = (fPlus - fMinus)/2.0/h;
      PLAYA_MSG2(verbosity,
        "g=" << setw(5) << globalIndex << ", l=" << setw(5) << localIndex 
        << " f0=" << setw(12) << f0 
         << " fPlus=" << setw(12) << fPlus 
        << " fMinus=" << setw(12) << fMinus << " df_dx="
        << setw(12) << df_dx[localIndex]);
      PLAYA_MSG3(verbosity, "i " << globalIndex << " x_i=" << tmp 
        << " f(x)=" << f0 
        << " f(x+h)=" << fPlus 
        << " f(x-h)=" << fMinus
        );
      /* restore point */
      x[localIndex] = tmp;
    }
  }
  
  double localMaxErr = 0.0;

  VectorSpace<double> space = x.space();
  for (int i=0; i<space.numLocalElements(); i++)
  {
    double num =  fabs(df_dx[i]-gf[i]);
    double den = fabs(df_dx[i]) + fabs(gf[i]) + 1.0e-14;
    double r = 0.0;
    if (fabs(den) > 1.0e-16) r = num/den;
    else r = 1.0;
    PLAYA_MSG2(verbosity, "i=" << setw(5) << i
      << " FD=" << setw(16) << df_dx[i] 
      << " grad=" << setw(16) << gf[i]
      << " r=" << setw(16) << r);
    if (localMaxErr < r) localMaxErr = r;
  }
  PLAYA_MSG2(verbosity, "local max error = " << localMaxErr);
  double maxErr = localMaxErr;
  PLAYA_MSG1(verbosity, "fd check: max error = " << maxErr);
  return maxErr <= tol;
}

}
