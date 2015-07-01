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


#include <iostream>
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaMPISession.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaVectorSpaceImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaPCGSolver.hpp"
#include "PlayaICCPreconditionerFactory.hpp"
#include "PlayaMatrixMarketIO.hpp"

using std::cout;
using std::endl;
using std::setw;
using std::exception;
using namespace Playa;

int main(int argc, char** argv)
{
  int rtn = 0;
  try
    {
      MPISession::init(&argc, &argv);     
      Tabs::showDepth() = false;
      

      LinearOperator<double> A = readMatrixMarketFile("nut-stiffness-n-2-p-1.mtx");
      VectorSpace<double> spc = A.domain();

      Vector<double> ans = spc.createMember();
      ans.randomize();

      Vector<double> b = A*ans;

      ParameterList cgParams("Linear Solver");
      cgParams.set("Verbosity", 2);
      cgParams.set("Max Iterations", 100);
      cgParams.set("Tolerance", 1.0e-6);

      PreconditionerFactory<double> prec
      	= new ICCPreconditionerFactory<double>();
      LinearSolver<double> solver = new PCGSolver(cgParams, prec);

      Vector<double> x = b.copy();
      SolverState<double> state = solver.solve(A, b, x);

      if (state.finalState()==SolveConverged)
	{
	  Out::root() << "solve succeeded!" << endl;
	  Out::root() << "|error|=" << norm2(x-ans) << endl;
	}
      else
	{
	  Out::root() << "solve FAILED! Message=" << state.finalMsg() << endl;
	  rtn = -1;
	}
    }
  catch(std::exception& e)
    {
      cerr << "exception detected: " << e.what() << endl;
      rtn = -1;
    }

  MPISession::finalize();
  return rtn;
}
