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
#include "PlayaGlobalAnd.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "PlayaSerialVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaMultiVectorOperatorDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaVectorOpsImpl.hpp"
#include "PlayaMultiVectorOperatorImpl.hpp"
#endif


using namespace Playa;
using namespace PlayaExprTemplates;

int main(int argc, char *argv[]) 
{
  int stat = 0;
  try
    {
      GlobalMPISession session(&argc, &argv);
 
      /* create a distributed vector space for the multivector's vectors */
      VectorType<double> rowType = new EpetraVectorType();
      int nLocalRows = 2;
      VectorSpace<double> space 
        = rowType.createEvenlyPartitionedSpace(MPIComm::world(), nLocalRows);
      
      /* create a replicated vector space for the small space of columns */
      int nVecs = 3;
      VectorType<double> colType = new SerialVectorType();
      VectorSpace<double> replSpace 
        = colType.createEvenlyPartitionedSpace(MPIComm::world(), nVecs);

      /* create some random vectors */
      Teuchos::Array<Vector<double> > vecs(nVecs);
      for (int i=0; i<nVecs; i++)
      {
        vecs[i] = space.createMember();
        vecs[i].randomize();
      }

      /* Test multiplication by a multivector operator. We will compute
       * y1 by directly summing columns, and y2 by applying the operator */
      LinearOperator<double> A = multiVectorOperator<double>(vecs, replSpace);

      
      Vector<double> y1 = space.createMember();
      Vector<double> y2 = space.createMember();
      y1.zero(); 
      y2.zero(); 
      

      /* Sum columns, putting the weights into x */
      Vector<double> x = replSpace.createMember();

      Out::os() << "A=" << A << std::endl;
      
      for (int j=0; j<replSpace.numLocalElements(); j++)
        {
          
          x[j] = 2.0*(drand48()-0.5);
          y1 = y1 + x[j] * vecs[j];
          Out::os() << "x[" << j << "]=" << x[j] << std::endl;
          Out::os() << "vecs[j]=" << vecs[j] << std::endl;
        }  
      Out::os() << "y1=" << std::endl << y1 << std::endl;
      
      /* Apply the operator to the vector of weights */
      y2 = A * x;
      Out::os() << "y2=A*x=" << std::endl << y2 << std::endl;

      Vector<double> y12 = y1-y2;
      Vector<double> y21 = y2-y1;
      Out::os() << "y1-y2=" << std::endl << y12 << std::endl;
      Out::os() << "y2-y1=" << std::endl << y21 << std::endl;
      
      double errA = (y1-y2).norm2();

      Out::root() << "error in A*x = " << errA << std::endl;


      /* Now test z = A^T * y */
      LinearOperator<double> At = A.transpose();
      
      Vector<double> z1 = replSpace.createMember();
      z1.zero();
      Vector<double> z2 = replSpace.createMember();
      z2.zero();

      Vector<double> y = y1.copy();

      /* compute by vectorwise multiplication */
      for (int j=0; j<replSpace.numLocalElements(); j++)
      {
        z1[j] = vecs[j].dot(y);
      }
      /* compute with operator */
      z2 = At * y;
      

      double errAt = (z1-z2).normInf();
      Out::root() << "error in At*y = " << errA << std::endl;

      double tol = 1.0e-13;
      bool pass = errA + errAt < tol;
      pass = globalAnd(pass);
      if (pass)
        {
          Out::root() << "multivector op test PASSED" << std::endl;
        }
      else
        {
          stat = -1;
          Out::root() << "multivector op test FAILED" << std::endl;
        }
    }
  catch(std::exception& e)
    {
      stat = -1;
      std::cerr << "Caught exception: " << e.what() << std::endl;
    }
  return stat;
}



