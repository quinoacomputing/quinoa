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

#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaPoissonBoltzmannOp.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaNOXSolver.hpp"
#include "PlayaLinearCombinationImpl.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"

#endif

using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;




int main(int argc, char *argv[]) 
{
  try
  {
    GlobalMPISession session(&argc, &argv);


    MPIComm::world().synchronize();

    /* create the nonlinear operator */
    VectorType<double> type = new EpetraVectorType();
    int nProc = MPIComm::world().getNProc();
    int nLocalRows = 128/nProc;
    PoissonBoltzmannOp* prob = new PoissonBoltzmannOp(nLocalRows, type);
    NonlinearOperator<double> F = prob;

    /* create the nox solver */
    ParameterXMLFileReader reader("nox.xml");
    ParameterList noxParams = reader.getParameters();

    Out::root() << "solver params = " << noxParams << std::endl;

    NOXSolver solver(noxParams);

    Vector<double> soln;
    SolverState<double> stat = solver.solve(F, soln);
    TEUCHOS_TEST_FOR_EXCEPTION(stat.finalState() != SolveConverged,
      runtime_error, "solve failed");

    Out::root() << "numerical solution = " << std::endl;
    Out::os() << soln << std::endl;

    Vector<double> exact = prob->exactSoln();

    Out::root() << "exact solution = " << std::endl;
    Out::os() << exact << std::endl;

//bvbw reddish port hack
    double temp_val = nLocalRows*nProc;
    double err = (exact-soln).norm2()/sqrt(temp_val);
    Out::root() << "error norm = " << err << std::endl;
      

    double tol = 1.0e-6;
    if (err > tol)
    {
      Out::root() << "NOX Poisson-Boltzmann test FAILED" << std::endl;
      return 1;
    }
    else
    {
      Out::root() << "NOX Poisson-Boltzmann test PASSED" << std::endl;
      return 0;
    }
  }
  catch(std::exception& e)
  {
    Out::root() << "Caught exception: " << e.what() << std::endl;
    return -1;
  }
}

