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


#include "SundanceEquationSet.hpp"
#include "SundanceExpr.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceUnknownParameter.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceQuadratureFamilyStub.hpp"
#include "SundanceCellFilterStub.hpp"
#include "SundanceIntegral.hpp"



using namespace Sundance;
using namespace Teuchos;

using Sundance::List;

int main()
{
  Expr dx = new Derivative(0);
  Expr dy = new Derivative(1);

  Expr u = new UnknownFunctionStub("u");
  Expr alpha = new UnknownParameter("alpha");
  Expr alpha0 = new Parameter(3.14, "alpha0");
  Expr beta = new UnknownParameter("beta");
  Expr beta0 = new Parameter(2.72, "beta0");
  Expr v = new TestFunctionStub("v");

  Out::os() << "u=" << u << std::endl;
  Out::os() << "v=" << v << std::endl;
  Out::os() << "alpha=" << alpha << std::endl;

  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);

  Expr u0 = new DiscreteFunctionStub("u0");
  Expr zero = new ZeroExpr();

  RCP<CellFilterStub> cells = rcp(new CellFilterStub());
  RCP<QuadratureFamilyStub> quad = rcp(new QuadratureFamilyStub(1));

  WatchFlag watchMe("watch eqn");
  watchMe.setParam("symbolic preprocessing", 6);
  

  Expr w = Integral(cells, v*u*alpha, quad, watchMe);
  Expr dum;
  Array<Expr> dum2;

  EquationSet eqn(w, dum, tuple(v), tuple(u), tuple(u0), 
    dum, dum, alpha, alpha0,  dum2, dum2);

  Out::os() << "num unk params=" << eqn.numUnkParams() << std::endl;
  Out::os() << "num fixed params=" << eqn.numFixedParams() << std::endl;

  for (int r=0; r<eqn.numRegions(); r++)
  {
    const RegionQuadCombo& rqc = eqn.regionQuadCombos()[r];
    const DerivSet& derivs = eqn.nonzeroFunctionalDerivs(Sensitivities, rqc);

    for (DerivSet::const_iterator i=derivs.begin(); i!=derivs.end(); i++)
    {
      const MultipleDeriv& d = *i;
      Out::os() << "d=" << d << std::endl;
      for (MultipleDeriv::const_iterator j=d.begin(); j!=d.end(); j++)
      {
        Out::os() << "j=" << *j;
        if (j->isParameter()) Out::os() << " (parameter) ";
        Out::os() << std::endl;
      }
    }
  }
}


