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
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaRandomSparseMatrixBuilderDecl.hpp"
#include "PlayaCompoundTester.hpp"
#include "PlayaMatrixMatrixTester.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#endif

STREAM_OUT(Vector<double>)

using namespace Playa;
using namespace PlayaExprTemplates;

int main(int argc, char *argv[]) 
{
  int stat = 0;
  try
    {
      GlobalMPISession session(&argc, &argv);
 

      VectorType<double> type = new EpetraVectorType();

      int nLocalRows = 4;
      
      double onProcDensity = 0.5;
      double offProcDensity = 0.1;

      RandomSparseMatrixBuilder<double> ABuilder(nLocalRows, nLocalRows, 
        onProcDensity, offProcDensity, type);
      RandomSparseMatrixBuilder<double> BBuilder(nLocalRows, nLocalRows, 
        onProcDensity, offProcDensity, type);

      /* Build some rectangular matrices to test products */
      RandomSparseMatrixBuilder<double> CBuilder(2*nLocalRows, nLocalRows, 
        onProcDensity, offProcDensity, type);
      RandomSparseMatrixBuilder<double> DBuilder(3*nLocalRows, 2*nLocalRows, 
        onProcDensity, offProcDensity, type);

      LinearOperator<double> A = ABuilder.getOp();
      LinearOperator<double> B = BBuilder.getOp();

      LinearOperator<double> C = CBuilder.getOp();
      LinearOperator<double> D = DBuilder.getOp();

      Out::root() << "A = " << std::endl;
      Out::os() << A << std::endl;
      Out::root() << "B = " << std::endl;
      Out::os() << B << std::endl;

      Out::root() << "C = " << std::endl;
      Out::os() << C << std::endl;
      Out::root() << "D = " << std::endl;
      Out::os() << D << std::endl;
      
      CompoundTester<double> tester(A, B, 
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10));

      bool allPass =  tester.runAllTests();

      Out::root() << std::endl << std::endl 
                  << "testing multiplication of square matrices " 
                  << std::endl << std::endl;

      MatrixMatrixTester<double> mmTester(A, B, 
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10));

      allPass = mmTester.runAllTests() && allPass;

      Out::root() << std::endl << std::endl 
                  << "testing multiplication of rectangular matrices " 
                  << std::endl << std::endl;

      MatrixMatrixTester<double> rectMMTester(C, D, 
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10),
        TestSpecifier<double>(true, 1.0e-13, 1.0e-10));

      allPass = rectMMTester.runAllTests() && allPass;

     if (!allPass) stat = -1;
    }
  catch(std::exception& e)
    {
      stat = 0;
      std::cerr << "Caught exception: " << e.what() << std::endl;
    }
  return stat;
}



