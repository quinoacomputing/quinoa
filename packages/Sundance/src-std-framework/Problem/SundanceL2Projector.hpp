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

#ifndef SUNDANCE_L2PROJECTOR_H
#define SUNDANCE_L2PROJECTOR_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceLinearProblem.hpp"
#include "SundanceCoordinateSystem.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaVectorType.hpp"

namespace Sundance
{
using namespace Playa;
using namespace Teuchos;

/**
 * L2Projector projects an expression onto a DiscreteSpace. 
 */
class L2Projector
{
public:
  /** */
  L2Projector(){;}
  /** */
  L2Projector(const DiscreteSpace& space, 
    const Expr& expr);
  /** */
  L2Projector(const DiscreteSpace& space, 
    const Expr& expr,
    const QuadratureFamily& quad);
  /** */
  L2Projector(const DiscreteSpace& space, 
    const Expr& expr,
    const LinearSolver<double>& solver);
  /** */
  L2Projector(const DiscreteSpace& space, 
    const Expr& expr,
    const LinearSolver<double>& solver,
    const QuadratureFamily& quad);
  /** */
  L2Projector(const DiscreteSpace& space, 
    const CoordinateSystem& coordSys,
    const Expr& expr);
  /** */
  L2Projector(const DiscreteSpace& space, 
    const CoordinateSystem& coordSys,
    const Expr& expr,
    const LinearSolver<double>& solver);

  /** */
  Expr project() const {return prob_.solve(solver_);}

  /** */
  const LinearProblem& prob() const {return prob_;}

private:

  void init(const DiscreteSpace& space,       
    const CoordinateSystem& coordSys,
    const Expr& expr,
    const LinearSolver<double>& solver,
    const QuadratureFamily& quad);
    
  LinearProblem prob_;

  LinearSolver<double> solver_;
    
};
}


#endif
