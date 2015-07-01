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


#include "SundancePCDPreconditioner.hpp"
#include "Sundance.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaSimpleIdentityOpDecl.hpp"
#include "PlayaSimpleIdentityOpImpl.hpp"
#include "PlayaSimpleComposedOpDecl.hpp"
#include "PlayaSimpleComposedOpImpl.hpp"
#include "PlayaSimpleBlockOpDecl.hpp"
#include "PlayaSimpleBlockOpImpl.hpp"
#include "PlayaSimpleScaledOpDecl.hpp"
#include "PlayaSimpleScaledOpImpl.hpp"
#include "PlayaGenericRightPreconditioner.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaInverseOperatorDecl.hpp"
#include "PlayaInverseOperatorImpl.hpp"

using namespace Playa;
using namespace Sundance;

PCDPreconditionerFactory::PCDPreconditionerFactory(
  const ParameterList& params,
  const LinearProblem& MpProb,
  const LinearProblem& ApProb,
  const LinearProblem& FpProb
  )
  : MpProb_(MpProb),
    ApProb_(ApProb),
    FpProb_(FpProb),
    MpSolver_(),
    ApSolver_(),
    FSolver_()
{
  ParameterList msParams = params.sublist("MpSolver");
  MpSolver_ = LinearSolverBuilder::createSolver(msParams);
  ParameterList asParams = params.sublist("ApSolver");
  ApSolver_ = LinearSolverBuilder::createSolver(asParams);
  ParameterList fsParams = params.sublist("FSolver");
  FSolver_ = LinearSolverBuilder::createSolver(fsParams);
}

Preconditioner<double> 
PCDPreconditionerFactory::
createPreconditioner(const LinearOperator<double>& K) const
{
  Tabs tab;

  LinearOperator<double> F = K.getBlock(0,0);
//  F.setName("F");
  LinearOperator<double> FInv = inverse(F, FSolver_);
//  FInv.setName("FInv");
  LinearOperator<double> Bt = K.getBlock(0,1);
//  Bt.setName("Bt");


  LinearOperator<double> Fp = FpProb_.getOperator();

  LinearOperator<double> Mp = MpProb_.getOperator();
//  Mp.setName("Mp");

  LinearOperator<double> MpInv = inverse(Mp, MpSolver_);
//  MpInv.setName("MpInv");

  LinearOperator<double> Ap = ApProb_.getOperator();
//  Ap.setName("Ap");

  LinearOperator<double> ApInv = inverse(Ap, ApSolver_);
//  ApInv.setName("ApInv");


  VectorSpace<double> pDomain = Bt.domain();
  VectorSpace<double> uDomain = F.domain();

  LinearOperator<double> Iu = identityOperator(uDomain);
//  Iu.setName("Iu");
  LinearOperator<double> Ip = identityOperator(pDomain);
//  Ip.setName("Ip");

  LinearOperator<double> XInv = MpInv * Fp * ApInv;

  VectorSpace<double> rowSpace = K.range();
  VectorSpace<double> colSpace = K.domain();
   
  LinearOperator<double> Q1 = makeBlockOperator(colSpace, rowSpace);
//  Q1.setName("Q1");
  LinearOperator<double> Q2 = makeBlockOperator(colSpace, rowSpace);
  // Q2.setName("Q2");
  LinearOperator<double> Q3 = makeBlockOperator(colSpace, rowSpace);
  //Q3.setName("Q3");
   
  Q1.setBlock(0, 0, FInv);
  Q1.setBlock(1, 1, Ip);
  Q1.endBlockFill();
   
  Q2.setBlock(0, 0, Iu);
  Q2.setBlock(0, 1, -1.0*Bt);
  Q2.setBlock(1, 1, Ip);
  Q2.endBlockFill();
   
  Q3.setBlock(0, 0, Iu);
  Q3.setBlock(1, 1, -1.0*XInv);
  Q3.endBlockFill();
   
  LinearOperator<double> P1 = Q2 * Q3;
  LinearOperator<double> PInv = Q1 * P1;

  return new GenericRightPreconditioner<double>(PInv);
}


