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

#ifndef SUNDANCE_INTERNALSOLVEMANAGER_HPP
#define SUNDANCE_INTERNALSOLVEMANAGER_HPP

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceBlock.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaSolverState.hpp"

namespace Sundance
{

class LinearSolveDriver
{
public:
  /** */
  LinearSolveDriver() {}
    

  /** */
  Expr formSolutionExpr(const Array<Vector<double> >& solnVector,
    const Array<RCP<DiscreteSpace> >& solutionSpace,
    const Array<Array<string> >& names,
    int verb) const ;

  /** */
  SolverState<double> solve(const LinearSolver<double>& solver,
    const LinearOperator<double>& A,
    const Array<Vector<double> >& rhs,
    const Array<RCP<DiscreteSpace> >& solutionSpace,
    const Array<Array<string> >& names,
    int verb,
    Expr& soln) const ;


  /** */
  void writeIntoSolutionExpr(
    const Array<Vector<double> >& solnVector,
    Expr soln, int verb) const ;

  /** Filename for dump of bad matrix */
  static std::string& badMatrixFilename() 
    {static std::string rtn = "badMatrix.dat"; return rtn;}

  /** Filename for dump of bad vector */
  static std::string& badVectorFilename() 
    {static std::string rtn = "badVector.dat"; return rtn;}

  /** Whether a solve failure throws an exception */
  static bool& solveFailureIsFatal()
    {static bool rtn=true; return rtn;}

  /** Whether to dump a matrix upon solve failure */
  static bool& dumpBadMatrix()
    {static bool rtn=false; return rtn;}

  
  

private:
};

}

#endif
