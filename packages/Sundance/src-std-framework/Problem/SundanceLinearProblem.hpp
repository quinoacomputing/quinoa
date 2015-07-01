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

#ifndef SUNDANCE_LINEARPROBLEM_H
#define SUNDANCE_LINEARPROBLEM_H

#include "SundanceDefs.hpp"
#include "SundanceLinearSolveDriver.hpp"
#include "SundanceObjectWithVerbosity.hpp"

namespace Sundance
{
using namespace Teuchos;

class Assembler;

/** 
 * LinearProblem encapsulates all information needed to form
 * a discrete linear problem. 
 *
 */
class LinearProblem 
{
public:
  /** Empty ctor */
  LinearProblem();
    
  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * and a vector type. */
  LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const Expr& test, const Expr& unk, 
    const Playa::VectorType<double>& vecType
    );
    
  /** Construct with a mesh, equation set, bcs, and blocks of variables */
  LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const BlockArray& test, const BlockArray& unk);
    
  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * parameters, and a vector type. */
  LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const Expr& test, const Expr& unk, 
    const Expr& unkParams, const Expr& unkParamVals, 
    const Playa::VectorType<double>& vecType);
    
  /** Construct with a mesh, equation set, bcs, parameters, and blocks of
      variables */
  LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
    const BlockArray& test, const BlockArray& unk, 
    const Expr& unkParams, const Expr& unkParamVals);

  /** */
  LinearProblem(const RCP<Assembler>& assembler);

  /** Solve the problem using the specified solver algorithm */
  Expr solve(const LinearSolver<double>& solver) const ;

  /** Solve the problem, writing the solution into the given function */
  SolverState<double> solve(const LinearSolver<double>& solver,
    Expr& soln) const ;


  /** Return the multivector on the right-hand side of the linear equation */
  Array<Vector<double> > getRHS() const ;

  /** Return the vector on the right-hand side of the linear equation */
  Vector<double> getSingleRHS() const {return getRHS()[0];}

  /** Return the operator on the left-hand side of the equation */
  LinearOperator<double> getOperator() const ;

  /** Return the map from cells and functions to row indices */
  const RCP<DOFMapBase>& rowMap(int blockRow) const ;
    
  /** Return the map from cells and functions to column indices */
  const RCP<DOFMapBase>& colMap(int blockCol) const ;

  /** Return the discrete space in which solutions live */
  const Array<RCP<DiscreteSpace> >& solnSpace() const ;

    
  /** Return the set of row indices marked as 
   * essential boundary conditions */
  const RCP<Set<int> >& bcRows(int blockRow) const ;

  /** Return the number of block rows in the problem  */
  int numBlockRows() const ;

  /** Return the number of block cols in the problem  */
  int numBlockCols() const ;

  /** with this function we can force the assembler to reassemble the matrix */
  void reAssembleProblem() const ;

  /** */
  Expr formSolutionExpr(const Array<Vector<double> >& vec) const ;

  /** Flag indicating whether to stop on a solve failure */
  static bool& solveFailureIsFatal()
    {return LinearSolveDriver::solveFailureIsFatal();}
    

  /** Flag indicating whether to write out the matrix and vector
   * after a solve failure */
  static bool& dumpBadMatrix() 
    {return LinearSolveDriver::dumpBadMatrix();}

  /** Filename for dump of bad matrix */
  static std::string& badMatrixFilename() 
    {return LinearSolveDriver::badMatrixFilename();}

  /** Filename for dump of bad vector */
  static std::string& badVectorFilename() 
    {return LinearSolveDriver::badVectorFilename();}

    

private:

      
  /** */
  RCP<Assembler> assembler_;

  /** */
  mutable LinearOperator<double> A_;

  /** */
  mutable Array<Vector<double> > rhs_;

  /** */
  Array<Array<string> > names_;

  /** */
  LinearSolveDriver solveDriver_;

  /** */
  Expr params_;

};

}


#endif
