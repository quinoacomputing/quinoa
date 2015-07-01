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

#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearCombinationDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaInverseOperatorDecl.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "PlayaLoadableVector.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "PlayaAztecSolver.hpp"
#include "PlayaMatrixLaplacian1D.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"


#include "PlayaLinearCombinationImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"

#include "PlayaInverseOperatorImpl.hpp"
#endif

using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;


int main(int argc, char *argv[]) 
{
  typedef Teuchos::ScalarTraits<double> ST;

  int status=-1;

  try
  {
    GlobalMPISession session(&argc, &argv);

    MPIComm::world().synchronize();

    VectorType<double> type = new EpetraVectorType();

    ParameterXMLFileReader reader("userPrecParams.xml");

    ParameterList solverParams = reader.getParameters();
    ParameterList innerSolverParams = solverParams.sublist("Inner Solve");
    ParameterList outerSolverParams = solverParams.sublist("Outer Solve");

    /* create the range space  */
    int nLocalRows = solverParams.get<int>("nLocal");

    MatrixLaplacian1D builder(nLocalRows, type);

    LinearOperator<double> A = builder.getOp();

    Vector<double> x = A.domain().createMember();
    int myRank = MPIComm::world().getRank();
    int nProcs = MPIComm::world().getNProc();
      
    x.randomize();
    if (myRank==0) loadable(x)->setElement(0, 0);
    if (myRank==nProcs-1) loadable(x)->setElement(nProcs * nLocalRows - 1, 0.0);

    cout << "x=" << std::endl;
    x.print(cout);
      
    Vector<double> y = A*x;
    cout << "y=" << std::endl;
    y.print(cout);

    Vector<double> ans = A.range().createMember();

    LinearSolver<double> innerSolver 
      = LinearSolverBuilder::createSolver(innerSolverParams);

    LinearSolver<double> outerSolver 
      = LinearSolverBuilder::createSolver(outerSolverParams);

    /* call the setUserPrec() function to set the operator and solver 
     * to be used for preconditioning */
    outerSolver.setUserPrec(A, innerSolver);

    LinearOperator<double> AInv = inverse(A, outerSolver);
      

    ans = AInv * y;

    //      SolverState<double> state = solver.solve(A, y, ans);
      

      
    //      cout << state << std::endl;

    cout << "answer is " << std::endl;
    ans.print(cout);
      
    double err = (x-ans).norm2();
    cout << "error norm = " << err << std::endl;

    double tol = 1.0e-6;
    if (err > tol)
    {
      cout << "User-defined preconditioner test FAILED" << std::endl;
      status = 1;
    }
    else
    {
      cout << "User-defined preconditioner test PASSED" << std::endl;
      status = 0;
    }
  }
  catch(std::exception& e)
  {
    cout << "Caught exception: " << e.what() << std::endl;
    status = -1;
  }

  return status;
}

