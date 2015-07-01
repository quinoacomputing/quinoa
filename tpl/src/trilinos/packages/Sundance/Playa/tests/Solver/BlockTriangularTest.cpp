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
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "PlayaBlockTriangularSolverDecl.hpp"
#include "PlayaDefaultBlockVectorSpaceDecl.hpp"
#include "PlayaLinearCombinationDecl.hpp"
#include "PlayaBlockVectorSpaceDecl.hpp"
#include "PlayaAztecSolver.hpp"
#include "PlayaMatrixLaplacian1D.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaSimpleBlockOpDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaSimpleBlockOpImpl.hpp"
#include "PlayaDefaultBlockVectorImpl.hpp"
#include "PlayaBlockTriangularSolverImpl.hpp"

#endif


using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;


int main(int argc, char *argv[]) 
{
  typedef Teuchos::ScalarTraits<double> ST;

  try
    {
      GlobalMPISession session(&argc, &argv);
      

      MPIComm::world().synchronize();

      VectorType<double> type = new EpetraVectorType();

      /* create the range space  */
      int nLocalRows = 10;

      MatrixLaplacian1D builder(nLocalRows, type);

      LinearOperator<double> A = builder.getOp();

      int nBlocks = 3;
      Array<Vector<double> > x(nBlocks);
      Array<VectorSpace<double> > space(nBlocks);
      for (int i=0; i<nBlocks; i++)
        {
          space[i] = A.domain();
          x[i] = A.domain().createMember();
          x[i].randomize();
        }

      VectorSpace<double> blkSpace = blockSpace(space);

      LinearOperator<double> bigA = makeBlockOperator(blkSpace, blkSpace);
      Vector<double> bigRHS = blkSpace.createMember();
      Vector<double> bigX = blkSpace.createMember();
      
      for (int i=0; i<nBlocks; i++)
        {
          bigX.setBlock(i, x[i]);
          for (int j=i; j<nBlocks; j++)
            {
              MatrixLaplacian1D builder(nLocalRows, type);
              LinearOperator<double> Aij = builder.getOp();
              bigA.setBlock(i,j,Aij);
            }
        }
      bigA.endBlockFill();
      
      bigRHS = bigA * bigX;
      Vector<double> bigSoln = blkSpace.createMember();

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(Playa::searchForFile("SolverParameters/poissonParams.xml"));
#else
      ParameterXMLFileReader reader("poissonParams.xml");
#endif

      ParameterList solverParams = reader.getParameters();
      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);
      LinearSolver<double> blockSolver 
        = new BlockTriangularSolver<double>(solver);
      
      SolverState<double> state = blockSolver.solve(bigA, bigRHS, bigSoln);
      
      std::cerr << state << std::endl;

      double err = (bigSoln - bigX).norm2();
      std::cerr << "error norm = " << err << std::endl;

      double tol = 1.0e-8;
      if (err > tol)
        {
          std::cerr << "Poisson solve test FAILED" << std::endl;
          return 1;
        }
      else
        {
          std::cerr << "Poisson solve test PASSED" << std::endl;
          return 0;
        }
    }
  catch(std::exception& e)
    {
      std::cerr << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
  return 0;
}

