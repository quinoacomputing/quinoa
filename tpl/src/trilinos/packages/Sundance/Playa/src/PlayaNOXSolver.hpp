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

#ifndef PLAYA_NOXSOLVER_HPP
#define PLAYA_NOXSOLVER_HPP

#include "PlayaDefs.hpp"
#include "PlayaNonlinearSolverBase.hpp"
#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_Playa_Group.hpp"
#include "NOX_Playa_StatusTestBuilder.hpp"
#include "NOX_Multiphysics_Solver_Manager.H"
#include "Teuchos_Assert.hpp"   
#include "Teuchos_ParameterList.hpp"

namespace Playa
{
using namespace Teuchos;
using std::runtime_error;


/**
 * Playa wrapper for NOX solver
 */
class NOXSolver : public NonlinearSolverBase<double>
{
public:
  /** */
  NOXSolver(){;}
  /** */
  NOXSolver(const ParameterList& params);
  /** */
  NOXSolver(const ParameterList& nonlinParams,
    const LinearSolver<double>& linSolver);

  /** */
  SolverState<double>  solve(const NonlinearOperator<double>& F, 
    Vector<double>& soln) const ;

  /** */
  const LinearSolver<double>& linSolver() const 
    {return linSolver_;}

  /* */
  GET_RCP(NonlinearSolverBase<double>);


private:

  LinearSolver<double> linSolver_;
  mutable RCP<NOX::StatusTest::Generic> statusTest_;
  mutable ParameterList params_;
  mutable ParameterList printParams_;
};
}

#endif
