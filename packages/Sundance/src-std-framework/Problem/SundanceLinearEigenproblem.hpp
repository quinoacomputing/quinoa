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

#ifndef SUNDANCE_LINEAREIGENPROBLEM_H
#define SUNDANCE_LINEAREIGENPROBLEM_H

#include "SundanceDefs.hpp"
#include "SundanceLinearProblem.hpp"
#include "SundanceEigensolution.hpp"
#include "SundanceFunctionalEvaluator.hpp"
#include "PlayaEigensolver.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace Playa;


/**
 *
 */
class LinearEigenproblem
{
public:
  /** */
  LinearEigenproblem(){}
  /** */
  LinearEigenproblem(
    const Mesh& mesh, const Expr& eqn,
    const Expr& v, const Expr& u,
    const VectorType<double>& vecType) ;
  /** */
  LinearEigenproblem(
    const Mesh& mesh, const Expr& eqn,
    const Expr& v, const Expr& u,
    const VectorType<double>& vecType,
    bool lumpMass) ;
  /** */
  LinearEigenproblem(
    const Mesh& mesh, const Expr& eqn,
    const Expr& massExpr,
    const Expr& v, const Expr& u,
    const VectorType<double>& vecType,
    bool lumpMass) ;

  /** */
  Eigensolution solve(const Eigensolver<double>& solver) const ;

  /** */
  LinearOperator<double> getK() const {return kProb_.getOperator();}

  /** */
  LinearOperator<double> getM() const {return mProb_.getOperator();}
    
  
private:
  /** */
  Array<Expr> makeEigenfunctions(Array<Vector<double> >& ev) const ;

  /** */
  LinearProblem makeMassProb(
    const Mesh& mesh,
    const Expr& massExpr,
    const Expr& v, const Expr& u,
    const VectorType<double>& vecType) const ;

  /** */
  LinearOperator<double> 
  lumpedOperator(const LinearOperator<double>& M) const ;

  bool lumpMass_;
  LinearProblem kProb_;
  LinearProblem mProb_;
  LinearOperator<double> M_;
  LinearOperator<double> MUnlumped_;
  DiscreteSpace discSpace_;
};




}


#endif
