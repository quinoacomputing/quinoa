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

#ifndef SUNDANCE_FUNCTIONAL_H
#define SUNDANCE_FUNCTIONAL_H

#include "SundanceDefs.hpp"
#include "SundanceLinearProblem.hpp"
#include "SundanceVectorCalculus.hpp"
#include "SundanceNonlinearProblem.hpp"
#include "SundanceFunctionalEvaluator.hpp"
#include "PlayaNonlinearOperator.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaVectorType.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;

/**
 *
 */
class Functional
{
public:
  /** */
  Functional(){;}

  /** */
  Functional(
    const Mesh& mesh, 
    const Expr& integral, 
    const Playa::VectorType<double>& vecType);

  /** */
  Functional(
    const Mesh& mesh, 
    const Expr& integral, 
    const Expr& essentialBC,
    const Playa::VectorType<double>& vecType);

  /** */
  LinearProblem linearVariationalProb(const Expr& var,
    const Expr& varEvalPts,
    const Expr& unk,
    const Expr& fixed,
    const Expr& fixedEvalPts) const ;

    
  /** */
  NonlinearProblem
  nonlinearVariationalProb(const Expr& var,
    const Expr& varEvalPts,
    const Expr& unk,
    const Expr& unkEvalPts,
    const Expr& fixed,
    const Expr& fixedEvalPts) const ;


  /** */
  FunctionalEvaluator evaluator(const Expr& var,
    const Expr& varEvalPts,
    const Expr& fixed,
    const Expr& fixedEvalPts) const ;


  /** */
  FunctionalEvaluator evaluator(const Expr& var,
    const Expr& varEvalPts) const ;

  /** */
  const Mesh& mesh() const {return mesh_;}
    

private:
  Mesh mesh_;

  Expr integral_;

  Expr bc_;

  Playa::VectorType<double> vecType_;
    
};

/** \relates Functional */
double L2Norm(const Mesh& mesh, const CellFilter& domain,
  const Expr& expr, const QuadratureFamily& quad,
  const WatchFlag& watch=WatchFlag());

/** \relates Functional */
double H1Seminorm(
  const Mesh& mesh,
  const CellFilter& filter,
  const Expr& f,
  const QuadratureFamily& quad,
  const WatchFlag& watch=WatchFlag());

/** \relates Functional */
double H1Norm(
  const Mesh& mesh,
  const CellFilter& filter,
  const Expr& f,
  const QuadratureFamily& quad,
  const WatchFlag& watch=WatchFlag());
}


#endif
