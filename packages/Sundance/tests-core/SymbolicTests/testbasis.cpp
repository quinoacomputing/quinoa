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


#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "Sundance.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceSpectralExpr.hpp"

int main(int argc, void** argv)
{
  try
    {
      Sundance::init(&argc, &argv);
      
      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      int ndim = 2;
      int order = 2;

      SpectralBasis SB(ndim, order);
     
      Array<Expr> Ex(6);
      Ex[0] = x;
      Ex[1] = y;
      Ex[2] = x*y;
      Ex[3] = 0.0;
      Ex[4] = x*x;
      Ex[5] = y*y;

      Expr SE = new SpectralExpr(SB, Ex);

      const SpectralExpr* se = dynamic_cast<const SpectralExpr*>(SE.ptr().get());
      SpectralBasis basis = se->getSpectralBasis();

      double e;

      for (int i=0; i< basis.nterms(); i++)
	for (int j=0; j< basis.nterms(); j++)
	  for(int k=0; k< basis.nterms(); k++)
	    {
	      e = basis.expectation(basis.getElement(i), basis.getElement(j), basis.getElement(k));
	      cout << i << " " << j  << " " << k << " " << e << std::endl;
	    }

    }
  
  catch(std::exception& e)
    {
      Sundance::handleException(e);
    }
  Sundance::finalize();
  return Sundance::testStatus();
}

