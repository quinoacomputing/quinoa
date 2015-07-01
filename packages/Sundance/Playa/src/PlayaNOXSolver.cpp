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


// $Id$ 
// $Source$ 


//   


#include "PlayaNOXSolver.hpp"         
#include "NOX_StatusTest_SafeCombo.hpp"         
#include "NOX.H"         
//#include "NOX_Parameter_Teuchos2NOX.H"         
#include "PlayaLinearSolverBuilder.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaExceptions.hpp"
#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif


using namespace NOX;
using namespace NOX::NOXPlaya;
using namespace Teuchos;
using namespace Playa;
using std::runtime_error;
using std::cout;
using std::ostream;


static Time& noxSolverTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("NOX solve"); 
  return *rtn;
}

NOXSolver::NOXSolver(const ParameterList& params)
  : linSolver_(),
    statusTest_(),
    params_(),
    printParams_()
{
  TEUCHOS_TEST_FOR_EXCEPTION(!params.isSublist("NOX Solver"), runtime_error,
                     "did not find NOX Solver sublist in " << params);
  
  params_ = params.sublist("NOX Solver");
  /* NOX wants to have the process ID in a parameter list???? */
  params_.sublist("Printing").set("MyPID", MPIComm::world().getRank());

  if (params_.isSublist("Status Test"))
    {
      statusTest_ = StatusTestBuilder::makeStatusTest(params_);
    }
  else
    {
      RCP<StatusTest::Generic> A = rcp(new StatusTest::NormF(1.0e-12));
      RCP<StatusTest::Generic> B = rcp(new StatusTest::MaxIters(20));
      statusTest_ = 
        rcp(new StatusTest::SafeCombo(StatusTest::SafeCombo::OR, A, B));
    }
  
  if (params_.isSublist("Linear Solver"))
    {
      linSolver_ = LinearSolverBuilder::createSolver(params_);
    }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(!params_.isSublist("Linear Solver"),
      RuntimeError, "no linear solver specified in NOX parameters");
  }
  
  if (params_.isSublist("Printing"))
    {
      printParams_ = params_.sublist("Printing");
    }
  
  TEUCHOS_TEST_FOR_EXCEPTION(linSolver_.ptr().get()==0, runtime_error,
                     "null linear solver object in NOXSolver ctor");

  TEUCHOS_TEST_FOR_EXCEPTION(statusTest_.get()==0, runtime_error,
                     "null status test object in NOXSolver ctor");

}

NOXSolver::NOXSolver(const ParameterList& nonlinParams,
      const LinearSolver<double>& linSolver)
  : linSolver_(linSolver),
    statusTest_(),
    params_(),
    printParams_()
{
  Tabs tab(0);
  TEUCHOS_TEST_FOR_EXCEPTION(!nonlinParams.isSublist("NOX Solver"), runtime_error,
                     "did not find NOX Solver sublist in " << nonlinParams);
  
  params_ = nonlinParams.sublist("NOX Solver");
  /* NOX wants to have the process ID in a parameter list???? */
  params_.sublist("Printing").set("MyPID", MPIComm::world().getRank());

  if (params_.isSublist("Status Test"))
    {
      statusTest_ = StatusTestBuilder::makeStatusTest(params_);
    }
  else
    {
      RCP<StatusTest::Generic> A = rcp(new StatusTest::NormF(1.0e-12));
      RCP<StatusTest::Generic> B = rcp(new StatusTest::MaxIters(20));
      statusTest_ = 
        rcp(new StatusTest::SafeCombo(StatusTest::SafeCombo::OR, A, B));
    }
  
  if (params_.isSublist("Linear Solver"))
    {
      Out::root() << tab << "WARNING: linear solver in NOX parameter list "
        "will be overridden by alternate solver" << std::endl;
    }
  
  if (params_.isSublist("Printing"))
    {
      printParams_ = params_.sublist("Printing");
    }
  
  TEUCHOS_TEST_FOR_EXCEPTION(linSolver_.ptr().get()==0, runtime_error,
                     "null linear solver object in NOXSolver ctor");

  TEUCHOS_TEST_FOR_EXCEPTION(statusTest_.get()==0, runtime_error,
                     "null status test object in NOXSolver ctor");

}




SolverState<double>
NOXSolver::solve(const NonlinearOperator<double>& F, 
                 Playa::Vector<double>& solnVec) const 
{
  TimeMonitor timer(noxSolverTimer());

  Vector<double> x0 = F.getInitialGuess();
  RCP<NOX::NOXPlaya::Group> grp = rcp(new NOX::NOXPlaya::Group(x0, F, linSolver_));
  RCP<Teuchos::ParameterList> noxParams 
    = Teuchos::rcp(&params_, false);
  RCP<NOX::Solver::Generic> solver 
    = NOX::Solver::buildSolver(grp, statusTest_, noxParams);

  NOX::StatusTest::StatusType rtn = solver->solve();


  const NOX::NOXPlaya::Group* solnGrp 
    = dynamic_cast<const NOX::NOXPlaya::Group*>(&(solver->getSolutionGroup()));

  TEUCHOS_TEST_FOR_EXCEPTION(solnGrp==0, runtime_error,
                     "Solution group could not be cast to NOX::NOXPlaya::Group");

  double resid = solnGrp->getNormF();
  int itersUsed = solver->getNumIterations();



  const NOX::NOXPlaya::Vector* x 
    = dynamic_cast<const NOX::NOXPlaya::Vector*>(&(solnGrp->getX()));

  TEUCHOS_TEST_FOR_EXCEPTION(x==0, runtime_error,
    "Solution vector could not be cast to NOX::NOXPlaya::Vector");
  
  solnVec = x->getPlayaVector();

  if (rtn==NOX::StatusTest::Converged)
  {
    return SolverState<double>(SolveConverged, "Solve converged", itersUsed, resid);
  }
  else if (rtn==NOX::StatusTest::Unconverged)
  {
    return SolverState<double>(SolveFailedToConverge, "Solve failed to converge", itersUsed, resid);
  }
  else
  {
    return SolverState<double>(SolveCrashed, "Solve crashed", itersUsed, resid);
  }

}
