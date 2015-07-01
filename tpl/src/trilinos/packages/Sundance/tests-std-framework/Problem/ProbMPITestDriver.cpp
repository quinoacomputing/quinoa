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


#include "Sundance.hpp"

#define DO_TEST(testName)\
  { \
    bool testName(); \
    bool pass = testName(); \
    if (pass) \
    {\
      numPass++; \
      Out::root() << "test " #testName " PASSED!" << endl; \
    } \
    else \
    { \
      numFail++; \
      Out::root() << "test " #testName " FAILED!" << endl; \
    }\
  }


int main(int argc, char** argv)
{
  bool allPass = true;

  try
  {
    Sundance::init(&argc, &argv);

    int numPass = 0;
    int numFail = 0;

    DO_TEST(BlockStochPoissonTest1D);
    DO_TEST(HighOrderPoisson2D);
    DO_TEST(HighOrderPoissonBernstein2D);
    DO_TEST(HighOrderProjection2D);
    DO_TEST(NonlinReducedIntegration);
    DO_TEST(SpectralPoisson1D);
    DO_TEST(SpectralSqrt);

    Out::root() 
      << "==================================================================="
      << endl
      << "==================================================================="
      << endl;

    if (numFail != 0)
    {
      Out::root() << "FAILURES detected!" << endl;
      Out::root() << numFail << " out of " << numFail + numPass 
                  << " tests failed" << endl;
    }
    else
    {
      Out::root() << "All tests PASSED!" << endl;
    }
    Out::root()
      << "==================================================================="
      << endl
      << "==================================================================="
      << endl;

  }
	catch(std::exception& e)
  {
    allPass = false;
    std::cerr << e.what() << std::endl;
  }
  Sundance::finalize();
  if (allPass) return 0;
  return -1;
}
