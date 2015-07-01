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

#include "SundanceL2Projector.hpp"
#include "PlayaAztecSolver.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"

#include "SundanceDerivative.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceGaussianQuadrature.hpp"

#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace Playa;


L2Projector::L2Projector(const DiscreteSpace& space, 
                         const Expr& expr, 
                         const LinearSolver<double>& solver)
  : prob_(), solver_()
{
  CoordinateSystem cs = new CartesianCoordinateSystem();
  init(space, cs, expr, solver, new GaussianQuadrature(4));
}


L2Projector::L2Projector(const DiscreteSpace& space, 
  const CoordinateSystem& cs,
  const Expr& expr, 
  const LinearSolver<double>& solver)
  : prob_(), solver_()
{
  init(space, cs, expr, solver, new GaussianQuadrature(4));
}


L2Projector::L2Projector(const DiscreteSpace& space, 
  const Expr& expr, 
  const LinearSolver<double>& solver,
  const QuadratureFamily& quad)
  : prob_(), solver_()
{
  CoordinateSystem cs = new CartesianCoordinateSystem();
  init(space, cs, expr, solver, quad);
}

L2Projector::L2Projector(const DiscreteSpace& space, 
                         const Expr& expr)
  : prob_(), solver_()
{
  CoordinateSystem cs = new CartesianCoordinateSystem();

  /* Create an Aztec solver for solving the linear subproblems */
  std::map<int,int> azOptions;
  std::map<int,double> azParams;
  
  azOptions[AZ_solver] = AZ_cg;
  azOptions[AZ_precond] = AZ_dom_decomp;
  azOptions[AZ_subdomain_solve] = AZ_icc;
  azOptions[AZ_graph_fill] = 1;
  azOptions[AZ_max_iter] = 1000;
  azOptions[AZ_output] = AZ_none;
  azParams[AZ_tol] = 1.0e-13;
  
  LinearSolver<double> solver = new AztecSolver(azOptions,azParams);

  init(space, cs, expr, solver, new GaussianQuadrature(4));
}

L2Projector::L2Projector(const DiscreteSpace& space, 
  const Expr& expr,
  const QuadratureFamily& quad)
  : prob_(), solver_()
{
  CoordinateSystem cs = new CartesianCoordinateSystem();

  /* Create an Aztec solver for solving the linear subproblems */
  std::map<int,int> azOptions;
  std::map<int,double> azParams;
  
  azOptions[AZ_solver] = AZ_cg;
  azOptions[AZ_precond] = AZ_dom_decomp;
  azOptions[AZ_subdomain_solve] = AZ_icc;
  azOptions[AZ_graph_fill] = 1;
  azOptions[AZ_max_iter] = 1000;
  azOptions[AZ_output] = AZ_none;
  azParams[AZ_tol] = 1.0e-13;
  
  LinearSolver<double> solver = new AztecSolver(azOptions,azParams);

  init(space, cs, expr, solver, quad);
}


L2Projector::L2Projector(const DiscreteSpace& space, 
  const CoordinateSystem& cs,
  const Expr& expr)
  : prob_(), solver_()
{
  /* Create an Aztec solver for solving the linear subproblems */
  std::map<int,int> azOptions;
  std::map<int,double> azParams;
  
  azOptions[AZ_solver] = AZ_cg;
  azOptions[AZ_precond] = AZ_dom_decomp;
  azOptions[AZ_subdomain_solve] = AZ_icc;
  azOptions[AZ_graph_fill] = 1;
  azOptions[AZ_max_iter] = 1000;
  azOptions[AZ_output] = AZ_none;
  azParams[AZ_tol] = 1.0e-13;
  
  LinearSolver<double> solver = new AztecSolver(azOptions,azParams);

  init(space, cs, expr, solver, new GaussianQuadrature(4));
}




void L2Projector::init(const DiscreteSpace& space,        
  const CoordinateSystem& coordSys,
  const Expr& expr, 
  const LinearSolver<double>& solver,
  const QuadratureFamily& quad)
{
  TEUCHOS_TEST_FOR_EXCEPTION(space.basis().size() != expr.size(),
                     std::runtime_error,
                     "mismatched vector structure between basis and expr");
  
  TEUCHOS_TEST_FOR_EXCEPTION(space.basis().size() == 0,
                     std::runtime_error,
                     "Empty basis?");
  
  Expr v = new TestFunction(space.basis()[0], "dummy_v[0]");
  Expr u = new UnknownFunction(space.basis()[0], "dummy_u[0]");
  
  for (int i=1; i<space.basis().size(); i++)
    {
      v.append(new TestFunction(space.basis()[i], "dummy_v[" 
                                + Teuchos::toString(i)+"]"));
      u.append(new UnknownFunction(space.basis()[i], "dummy_u[" 
                                + Teuchos::toString(i)+"]"));
    }

  CellFilter interior = new MaximalCellFilter();

  Expr eqn = 0.0;
  Expr J = coordSys.jacobian();

  for (int i=0; i<space.basis().size(); i++)
    {
      eqn = eqn + Integral(space.cellFilters(i), 
        J*v[i]*(u[i]-expr[i]), 
        quad);
    }
  Expr bc;

  prob_ = LinearProblem(space.mesh(), eqn, bc, v, u, space.vecType());
  solver_ = solver;
}



