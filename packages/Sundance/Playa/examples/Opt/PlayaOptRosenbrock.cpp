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

#include "PlayaOptBuilder.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaOut.hpp"
#include "PlayaLinearCombinationImpl.hpp"

#include "PlayaRosenbrock.hpp"

#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaMPIComm.hpp"


/* ------------------------------------------------------------------------
 *
 * This example shows the unconstrained optimization of a simple objective
 * function. The function is defined in the files PlayaRosenbock.[hpp/cpp].
 * 
 * The optimizer object is created by the OptBuilder::createOptimizer()
 * function, which reads parameters from an XML file. In this example,
 * a limited-memory BFGS optmizer is used. Its parameters are read from the
 * file basicLMBFGS.xml.
 * 
 * ------------------------------------------------------------------------ */

using namespace Playa;

int main(int argc, char *argv[])
{
  int rtn = 0;

  try
  {
    /* Initialize MPI */
    GlobalMPISession session(&argc, &argv);


    /* The VectorType object will be used when we create vector spaces, 
     * specifying what type of low-level linear algebra implementation
     * will be used. */
    VectorType<double> vecType = new EpetraVectorType();

    /* Construct the objective function */
    int M = 6;
    double alpha = 40.0;
    RCP<ObjectiveBase> obj = rcp(new Rosenbrock(M, alpha, vecType));

    Out::root() << "Objective function is " << obj->description()
                << endl;

    /* Get the starting point for the optimization run */
    Vector<double> xInit = obj->getInit();
    
    /* Run a finite-difference calculation of the function's gradient. This
     * will be expensive, but is a valuable test when developing new objective
     * functions */
    Out::root() << "Doing FD check of gradient..." << endl;
    bool fdOK = obj->fdCheck(xInit, 1.0e-6, 0);
    TEUCHOS_TEST_FOR_EXCEPTION(!fdOK, std::runtime_error,
      "finite difference test of Rosenbrock function gradient FAILED");
    Out::root() << "FD check OK!" << endl << endl;


    /* Create an optimizer object. */
    RCP<UnconstrainedOptimizerBase> opt 
      = OptBuilder::createOptimizer("basicLMBFGS.xml");

    /* Set the optimizer's verbosity to a desired volume of output */
    opt->setVerb(1);

    /* Run the optimizer with xInit s the initial guess. 
     * The results (location and value of min) and 
     * diagnostic information about the run are stored in the OptState
     * object returned by the run() function. */
    OptState state = opt->run(obj, xInit);

    /* Check the output */
    if (state.status() != Opt_Converged)
    {
      Out::root() << "optimization failed: " << state.status() << endl;
      rtn = -1;
    }
    else /* We converged, so let's make sure we got the right solution */
    {
      Out::root() << "optimization converged!" << endl; 
      Out::root() << "Iterations taken: " << state.iter() << endl;
      Out::root() << "Approximate minimum value: " << state.fCur() << endl;

      /* The exact solution is [1,1,\cdots, 1, 1]. */
      Vector<double> exactAns = xInit.space().createMember();
      exactAns.setToConstant(1.0);
      
      /* Compute the norm of the error in the location of the minimum */
      double locErr = (exactAns - state.xCur()).norm2();
      Out::root() << "||x-x^*||=" << locErr << endl;
      
      /* Compare the error to a desired tolerance */
      double testTol = 1.0e-4;
      if (locErr > testTol) 
      {
        rtn = -1;
      }
    }
     
    if (rtn == 0)
    {
      Out::root() << "test PASSED" << endl;
    }
    else
    {
      Out::root() << "test FAILED" << endl;
    }
  }
  catch(std::exception& e)
  {
    std::cerr << "Caught exception: " << e.what() << endl;
    rtn = -1;
  }
  return rtn;
}
