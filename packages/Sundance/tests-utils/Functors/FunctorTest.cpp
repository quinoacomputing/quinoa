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


#include "SundanceStdMathFunctors.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace Sundance;
using namespace Teuchos;

template <class F> bool functorTest(int nx, double tol)
{
  RCP<UnaryFunctor> f = rcp(new F());

  return f->test(nx, tol);
}

bool powTest(double a, int nx, double tol)
{
  RCP<UnaryFunctor> f = rcp(new PowerFunctor(a));

  return f->test(nx, tol);
}



int main(int argc, char** argv)
{
  int stat = 0;
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      int nx = 5;
      double tol = 1.0e-6;

      UnaryFunctor::checkResults() = true;

      UnaryFunctor::fdStep() = 1.0e-3;

      bool isOK = functorTest<StdReciprocal>(nx, tol);

      isOK = functorTest<StdExp>(nx, tol) && isOK ;

      isOK = functorTest<StdLog>(nx, tol) && isOK ;

      isOK = functorTest<StdSqrt>(nx, tol) && isOK ;

      /* power functor */
      isOK = powTest(-2.0, nx, tol) && isOK;
      isOK = powTest(-2.5, nx, tol) && isOK;
      isOK = powTest(0.0, nx, tol) && isOK;
      isOK = powTest(0.5, nx, tol) && isOK;
      isOK = powTest(1.0, nx, tol) && isOK;
      isOK = powTest(2.0, nx, tol) && isOK;
      isOK = powTest(4.0, nx, tol) && isOK;

      /* trig functions */

      isOK = functorTest<StdSin>(nx, tol) && isOK ;

      isOK = functorTest<StdCos>(nx, tol) && isOK ;

      isOK = functorTest<StdTan>(nx, tol) && isOK ;

      isOK = functorTest<StdASin>(nx, tol) && isOK ;

      isOK = functorTest<StdACos>(nx, tol) && isOK ;

      isOK = functorTest<StdATan>(nx, tol) && isOK ;

      /* hyperbolic functions */

      isOK = functorTest<StdSinh>(nx, tol) && isOK ;

      isOK = functorTest<StdCosh>(nx, tol) && isOK ;

      isOK = functorTest<StdTanh>(nx, tol) && isOK ;

      isOK = functorTest<StdASinh>(nx, tol) && isOK ;

      isOK = functorTest<StdACosh>(nx, tol) && isOK ;

      isOK = functorTest<StdATanh>(nx, tol) && isOK ;

      std::cerr << "done with all tests!" << std::endl;

      if (isOK) 
        {
          std::cerr << "all tests PASSED" << std::endl;
        }
      else
      {
        stat = -1;
          std::cerr << "a test has FAILED" << std::endl;
        }

    }
	catch(std::exception& e)
		{
      stat = -1;
      std::cerr << "detected exception " << e.what() << std::endl;
		}

  return stat;
}
