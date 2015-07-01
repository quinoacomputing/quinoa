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



#include <cstdlib>
#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "PlayaEpetraVectorSpace.hpp"
#include "Teuchos_Time.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaEpetraMatrix.hpp"
#include "PlayaMatrixLaplacian1D.hpp"
#include "PlayaRandomSparseMatrixBuilderDecl.hpp"
#include "PlayaRandomBlockMatrixBuilderDecl.hpp"
#include "PlayaCompoundTester.hpp"
#include "PlayaOut.hpp"
#include "PlayaDefaultBlockVectorSpaceDecl.hpp"
#include "PlayaLinearCombinationImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaDefaultBlockVectorImpl.hpp"
#include "PlayaRandomBlockMatrixBuilderImpl.hpp"
#endif



using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;
using std::endl;

int main(int argc, char *argv[]) 
{
  int stat = 0;
  try
  {
    GlobalMPISession session(&argc, &argv);
    MPIComm::world().synchronize();

    Out::os() << "go!" << std::endl;
    VectorType<double> type = new EpetraVectorType();

    Array<int> domainBlockSizes = tuple(2,3,4);
    Array<int> rangeBlockSizes = tuple(2,2);

    Array<VectorSpace<double> > domainBlocks(domainBlockSizes.size());
    Array<VectorSpace<double> > rangeBlocks(rangeBlockSizes.size());

    for (int i=0; i<domainBlocks.size(); i++)
    {
      domainBlocks[i] = type.createEvenlyPartitionedSpace(MPIComm::world(),
        domainBlockSizes[i]);
    }

    for (int i=0; i<rangeBlocks.size(); i++)
    {
      rangeBlocks[i] = type.createEvenlyPartitionedSpace(MPIComm::world(),
        rangeBlockSizes[i]);
    }
      
    VectorSpace<double> domain = blockSpace(domainBlocks);
    VectorSpace<double> range = blockSpace(rangeBlocks);

    double blockDensity = 0.75;
    double onProcDensity = 0.5;
    double offProcDensity = 0.1;
      
    RandomBlockMatrixBuilder<double> builder(domain, range, 
      blockDensity,
      onProcDensity,
      offProcDensity,
      type);

    LinearOperator<double> A = builder.getOp();

    Out::os() << "A num block rows = " << A.numBlockRows() << std::endl;
    Out::os() << "A num block cols = " << A.numBlockCols() << std::endl;

    Vector<double> x = domain.createMember();
    Out::os() << "randomizing trial vector" << std::endl;
    x.randomize();

    Array<Vector<double> > xBlock(domain.numBlocks());
    for (int i=0; i<xBlock.size(); i++)
    {
      xBlock[i] = x.getBlock(i);
    }

    Vector<double> xx = x.copy();

      

    Out::os() << "------------------------------------------------------------" << std::endl;
    Out::os() << "computing A*x..." << std::endl;
    Vector<double> y0 = A * x;
    for (int i=0; i<y0.space().numBlocks(); i++)
    {
      Out::os() << "y0[" << i << "] = " << std::endl << y0.getBlock(i) << std::endl;
    }
      

    Vector<double> y1 = range.createMember();
    Out::os() << "------------------------------------------------------------" << std::endl;
    Out::os() << "computing A*x block-by-block..." << std::endl;
    Array<Vector<double> > yBlock(range.numBlocks());
    for (int i=0; i<yBlock.size(); i++)
    {
      yBlock[i] = range.getBlock(i).createMember();
      yBlock[i].zero();
      for (int j=0; j<xBlock.size(); j++)
      {
        LinearOperator<double> Aij = A.getBlock(i,j);
        if (Aij.ptr().get() != 0)
        {
          Out::os() << "A(" << i << ", " << j << ") = " << std::endl 
                    << Aij << std::endl;
        }
        else
        {
          Out::os() << "A(" << i << ", " << j << ") = 0 " << std::endl;
        }
        Out::os() << "x[" << j << "] = " << std::endl << xBlock[j] << std::endl;
        if (Aij.ptr().get()==0) continue;
        yBlock[i] = yBlock[i] + Aij * xBlock[j];
      }
      y1.setBlock(i, yBlock[i]);
    }

    for (int i=0; i<y1.space().numBlocks(); i++)
    {
      Out::os() << "y1[" << i << "] = " << std::endl << y1.getBlock(i) << std::endl;
    }
    double err = (y1 - y0).norm2();
    Out::os() << "error = " << err << std::endl;

    double tol = 1.0e-13;
    if (err < tol)
    {
      Out::os() << "block op test PASSED" << std::endl;
    }
    else
    {
      stat = -1;
      Out::os() << "block op test FAILED" << std::endl;
    }
  }
  catch(std::exception& e)
  {
    stat = -1;
    Out::os() << "Caught exception: " << e.what() << std::endl;
  }
  return stat;
}



