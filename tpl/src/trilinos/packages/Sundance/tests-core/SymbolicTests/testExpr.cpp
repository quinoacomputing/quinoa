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


#include "PlayaMPISession.hpp"
#include "SundanceExpr.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceHermiteSpectralBasis.hpp"
#include "SundanceSpectralExpr.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceCellFilterStub.hpp"
#include "SundanceQuadratureFamilyStub.hpp"

using std::cout;
using std::cerr;
using std::exception;
using namespace Teuchos;
using namespace Playa;
using namespace Sundance;

int main(int argc, char** argv)
{
  try
    {
      MPISession::init(&argc, (void***)&argv);
      

      Expr::showAllParens() = true;

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);

      Expr dx = new Derivative(0);

      int ndim = 2;
      int order = 2;

      SpectralBasis SB = new HermiteSpectralBasis(ndim, order);

      Expr u = new UnknownFunctionStub("u", SB);
      Expr v = new TestFunctionStub("v", SB);
      Expr w = new UnknownFunctionStub("w", SB);
     
      Array<Expr> Ex1(6);
      Ex1[0] = 1.0;
      Ex1[1] = x;
      Ex1[2] = 0.0;
      Ex1[3] = x*y;
      Ex1[4] = x+y;
      Ex1[5] = y;

      Expr SE1 = new SpectralExpr(SB, Ex1);


      Array<Expr> Ex2(6);
      Ex2[0] = -3.0*x;
      Ex2[1] = 0.0;
      Ex2[2] = -y;
      Ex2[3] = x-y;
      Ex2[4] = -4.0*y + 2.0*x;
      Ex2[5] = 0.0;

      Expr SE2 = new SpectralExpr(SB, Ex2);

      Expr G = x*x;

      Expr Sum  = (dx*v) * (dx*u) + v*x;

      Handle<CellFilterStub> domain = rcp(new CellFilterStub());
      Handle<QuadratureFamilyStub> quad = rcp(new QuadratureFamilyStub(1));

      Expr eqn = Integral(domain, Sum, quad);

      const SpectralExpr* se = dynamic_cast<const SpectralExpr*>(Sum.ptr().get());
      if (se != 0)
        {
          SpectralBasis basis = se->getSpectralBasis();

          for(int i=0; i< basis.nterms(); i++)
            cout << se->getCoeff(i) << std::endl;
        }

      cout << Sum << std::endl << std::endl; 

      cout << eqn << std::endl << std::endl; 


    }
  
  catch(std::exception& e)
    {
      std::cerr << e.what() << std::endl;
    }
  MPISession::finalize();
}

