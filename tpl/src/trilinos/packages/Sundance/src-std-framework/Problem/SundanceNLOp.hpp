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

#ifndef SUNDANCE_NLOP_H
#define SUNDANCE_NLOP_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceBlock.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "PlayaNonlinearOperatorBase.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaVectorType.hpp"

namespace Sundance
{
using namespace Teuchos;

class Assembler;


/** 
 * NLOp encapsulates a discrete nonlinear problem, and can
 * be passed to a nonlinear solver such as NOX.
 */
class NLOp 
  : public ObjectWithClassVerbosity<NLOp>,
    public Playa::NonlinearOperatorBase<double>
{
public:
  /** Empty ctor */
  NLOp();

  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * and a vector type */
  NLOp(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const Expr& test, const Expr& unk, const Expr& u0, 
    const Playa::VectorType<double>& vecType,
    bool partitionBCs = false);

  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * and a vector type */
  NLOp(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const BlockArray& test, const BlockArray& unk, const Expr& u0);


  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * parameters, and a vector type */
  NLOp(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const Expr& test, const Expr& unk, const Expr& u0, 
    const Expr& params, const Expr& paramVals,  
    const Playa::VectorType<double>& vecType,
    bool partitionBCs = false);


  /** */
  NLOp(const RCP<Assembler>& assembler, 
    const Expr& u0);

  /** Compute the residual and Jacobian at the current evaluation point */
  LinearOperator<double> computeJacobianAndFunction(Vector<double>& functionValue) const ;

  /** Write the Jacobian and residual into the objects provided */
  void computeJacobianAndFunction(LinearOperator<double>& J,
    Vector<double>& resid) const ;

  /** Compute direct sensitivities to parameters */
  Expr computeSensitivities(const LinearSolver<double>& solver) const ;
      

  /** Return the current evaluation point as a Sundance expression */
  Expr getU0() const {return u0_;}

  /** Compute the residual at the current eval point */
  Playa::Vector<double> computeFunctionValue() const ;
      
  /** Write the residual into the object provided */
  void computeFunctionValue(Vector<double>& resid) const ;

  /** Get an initial guess */
  Playa::Vector<double> getInitialGuess() const ;

  /** Create the Jacobian object, but don't fill it in. */
  LinearOperator<double> allocateJacobian() const ;

  /** Set an initial guess */
  void setInitialGuess(const Expr& u0New);

  /** This function forces the assembler to reassemble the matrix */
  void reAssembleProblem() const;

  /* Handle boilerplate */
  GET_RCP(Playa::NonlinearOperatorBase<double>);

protected:
  /** */
  void updateDiscreteFunctionValue(const Vector<double>& vec) const ;
private:
      
  /** */
  RCP<Assembler> assembler_;

  /** */
  mutable Playa::LinearOperator<double> J_;

  /** */
  Expr u0_;

  /** */
  Expr params_;

  /** */
  Expr paramVals_;
};
}




#endif
