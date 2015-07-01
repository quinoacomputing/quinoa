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
#include "PlayaLinearSolverBuilder.hpp"
#include "PlayaAmesosSolver.hpp"
#include "PlayaAztecSolver.hpp"
#include "PlayaBelosSolver.hpp"
#include "PlayaBICGSTABSolverDecl.hpp"
#include "PlayaBlockTriangularSolverDecl.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaBICGSTABSolverImpl.hpp"
#include "PlayaBlockTriangularSolverImpl.hpp"
#endif

using namespace Playa;
using namespace PlayaExprTemplates;
using namespace Teuchos;


LinearSolver<double> LinearSolverBuilder::createSolver(const std::string& filename)
{
  ParameterXMLFileReader reader(filename);
  ParameterList solverParams = reader.getParameters();
  return createSolver(solverParams);
}



LinearSolver<double> LinearSolverBuilder::createSolver(const ParameterList& params, int verb)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!params.isSublist("Linear Solver"), std::runtime_error,
    "did not find Linear Solver sublist in " << params);

  ParameterList solverSublist = params.sublist("Linear Solver");

  const std::string& solverType = getParameter<string>(solverSublist, "Type");

  Tabs tab;
  PLAYA_MSG1(verb, tab << "Solver builder creating a solver of type="
    << solverType);
  Tabs tab2;
  PLAYA_MSG2(verb, tab2 << "params = " << solverSublist);

  if (solverType=="Aztec")
  {
    return new AztecSolver(solverSublist);
  }
  else if (solverType=="Playa" || solverType=="TSF")
  {
    const std::string& solverMethod = getParameter<string>(solverSublist, "Method");
    if (solverMethod=="BICGSTAB") 
    {
      return new BICGSTABSolver<double>(solverSublist);
    }
    else if (solverMethod=="GMRES")
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, RuntimeError, "Playa GMRES solver not implemented");
    }
  }
  else if (solverType=="Amesos")
  {
    return new AmesosSolver(solverSublist);
  }
  else if (solverType=="Belos")
  {
    return new BelosSolver(solverSublist);
  }
  else if (solverType=="Block Triangular")
  {
    ParameterList subSolverParams = solverSublist.sublist("Sub Solver");
    LinearSolver<double> subSolver = createSolver(subSolverParams);
    return new BlockTriangularSolver<double>(subSolver);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, 
    "Could not create a solver from parameter list " 
    << params);
  return LinearSolver<double>();
    
}

