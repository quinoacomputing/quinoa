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

#include "PlayaExceptions.hpp"
#include "PlayaNonlinearSolverBuilder.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "PlayaNOXSolver.hpp"
#include "PlayaNewtonArmijoSolverImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#endif

using namespace Playa;
using namespace PlayaExprTemplates;
using namespace Teuchos;


NonlinearSolver<double> NonlinearSolverBuilder::createSolver(const std::string& filename)
{
  ParameterXMLFileReader reader(filename);
  ParameterList solverParams = reader.getParameters();
  return createSolver(solverParams);
}



NonlinearSolver<double> NonlinearSolverBuilder::createSolver(const ParameterList& params)
{
  if (params.isSublist("NOX Solver"))
  {
    return new NOXSolver(params);
  }
  else if (params.isSublist("Nonlinear Solver"))
  {
    ParameterList sub = params.sublist("Nonlinear Solver");
    Array<string> names = tuple<string>("Newton Armijo Solver", "Newton-Armijo Solver", "NewtonArmijoSolver");
    for (int i=0; i<names.size(); i++)
    {
      if (sub.isSublist(names[i]))
      {
        ParameterList subsub = sub.sublist(names[i]);
        LinearSolver<double> linSolver;
        if (subsub.isParameter("Linear Solver"))
        {
          string solverFile = subsub.get<string>("Linear Solver");
          linSolver = LinearSolverBuilder::createSolver(solverFile);
        }
        else if (subsub.isSublist("Linear Solver"))
        {
          linSolver = LinearSolverBuilder::createSolver(subsub);
        }
        else
        {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
            "Nonlinear solver parameter list " << sub
            << " does not appear to specify a solver for the linear subproblems");
        }
        return new NewtonArmijoSolver<double>(subsub, linSolver);
      }
    }
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
      "Nonlinear solver parameter list " << params
      << " can't be parsed to find a nonlinear solver");
  }

  return NonlinearSolver<double>();
    
}

