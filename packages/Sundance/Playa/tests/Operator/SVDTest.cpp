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
#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaSimpleDiagonalOpDecl.hpp"
#include "PlayaSimpleComposedOpDecl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaOut.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaVectorOpsImpl.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaSimpleDiagonalOpImpl.hpp"
#include "PlayaSimpleComposedOpImpl.hpp"
#include "PlayaSimpleTransposedOpImpl.hpp"
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

      VectorType<double> vecType = new SerialVectorType();

      VectorSpace<double> domain 
        = vecType.createEvenlyPartitionedSpace(MPIComm::world(), 3);
      VectorSpace<double> range
        = vecType.createEvenlyPartitionedSpace(MPIComm::world(), 5);

      RCP<MatrixFactory<double> > mf 
        = vecType.createMatrixFactory(domain, range);

      LinearOperator<double> A = mf->createMatrix();
      RCP<DenseSerialMatrix> APtr 
        = rcp_dynamic_cast<DenseSerialMatrix>(A.ptr());

      APtr->setRow(0, tuple(1.0,   2.0,  3.0));
      APtr->setRow(1, tuple(4.0,   5.0,  6.0));
      APtr->setRow(2, tuple(7.0,   8.0,  9.0));
      APtr->setRow(3, tuple(10.0, 11.0,  12.0));
      APtr->setRow(4, tuple(13.0, 14.0,  15.0));

      Out::os() << "A = " << std::endl;
      A.setVerb(10);
      Out::os() << A << std::endl;

      LinearOperator<double> U;
      LinearOperator<double> Vt;
      Vector<double> sigma;

      denseSVD(A, U, sigma, Vt);

      Out::os() << "U = " << std::endl;
      U.setVerb(10);
      Out::os() << U << std::endl;

      Out::os() << "sigma = " << std::endl;
      Out::os() << sigma << std::endl; 

      Out::os() << "Vt = " << std::endl;
      Vt.setVerb(10);
      Out::os() << Vt << std::endl;

      int nSamples = 10;
      bool allOK = true;
      double tol = 1.0e-13;
      for (int i=0; i<nSamples; i++)
      {
        Out::os() << "Sample #" << i << " of " << nSamples << std::endl;
        Vector<double> x = domain.createMember();
        x.randomize();
        
        U.setVerb(0);
        Vt.setVerb(0);
        A.setVerb(0);
        
        LinearOperator<double> Sigma = diagonalOperator(sigma);
        
        Vector<double> z = (U * Sigma * Vt)*x - A*x;
        double ez = z.norm2();
        Out::os() << "|| (U Sigma Vt - A)*x || = " << ez << std::endl;
        
        Vector<double> y = (U.transpose() * U)*x - x;
        double ey = y.norm2();
        Out::os() << "|| (U^T U - I)*x || = " << ey << std::endl;
        
        Vector<double> w = (Vt * Vt.transpose())*x - x;
        double ew = w.norm2();
        Out::os() << "|| (V^T*V - I)*x || = " << ew << std::endl;
        if (ew > tol || ez > tol || ey > tol) allOK = false;
      }

      if (allOK)
      {
        Out::os() << "SVD test PASSED" << std::endl;
      }
      else
      {
        Out::os() << "SVD test FAILED" << std::endl;
        stat = -1;
      }

      
    }
  catch(std::exception& e)
    {
      stat = -1;
      Out::os() << "Caught exception: " << e.what() << std::endl;
    }
  return stat;
}



