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
#include "PlayaInverseOperatorDecl.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "PlayaGlobalAnd.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "PlayaSerialVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "PlayaDenseLUSolver.hpp"
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


bool runit(const VectorType<double>& vecType,
  const LinearSolver<double>& solver)
{
  typedef Teuchos::ScalarTraits<double> ST;

  /* create the range space  */
  int nLocalRows = 10;
  
  MatrixLaplacian1D builder(nLocalRows, vecType);

  LinearOperator<double> A = builder.getOp();

  Out::root() << "matrix is " << std::endl;
  Out::os() << A << std::endl;

  Vector<double> x = A.domain().createMember();

  x.randomize();

  Out::root() << "input is " << std::endl;
  Out::os() << x << std::endl;
  Vector<double> y = A*x;

  Out::root() << "rhs is " << std::endl;
  Out::os() << y << std::endl;

  Vector<double> ans = A.range().createMember();
  
  Out::root() << "slot for solution is " << std::endl;
  Out::os() << ans << std::endl;

  LinearOperator<double> AInv = inverse(A, solver);

  ans = AInv * y;

  Out::root() << "answer is " << std::endl;
  Out::os() << ans << std::endl;
      
  double err = (x-ans).norm2();
  Out::root() << "error norm = " << err << std::endl;

  double tol = 1.0e-7;
    
  if (err <= tol)
  {
    Out::root() << "Poisson solve test PASSED" << std::endl;
    return true;
  }
  else
  {
    Out::root() << "Poisson solve test FAILED" << std::endl;
    return false;
  }
}


int main(int argc, char *argv[]) 
{
  bool status = 0;

  try
  {
    GlobalMPISession session(&argc, &argv);

    int nProc = session.getNProc();
    int rank = session.getRank();

    VectorType<double> epetra = new EpetraVectorType();
    VectorType<double> serial = new SerialVectorType();

    LinearSolver<double> denseLU = new DenseLUSolver();
    LinearSolver<double> amesos = LinearSolverBuilder::createSolver("amesos.xml");
    LinearSolver<double> belos_ml = LinearSolverBuilder::createSolver("belos-ml.xml");
    LinearSolver<double> belos_ifpack = LinearSolverBuilder::createSolver("belos-ifpack.xml");
    LinearSolver<double> aztec_ml = LinearSolverBuilder::createSolver("aztec-ml.xml");
    LinearSolver<double> aztec_ifpack = LinearSolverBuilder::createSolver("aztec-ifpack.xml");
    LinearSolver<double> bicgstab = LinearSolverBuilder::createSolver("bicgstab.xml");

    bool allOK = true;

    Out::root() << "Running Belos/ML" << std::endl;
    allOK = runit(epetra, belos_ml) && allOK;

    Out::root() << "Running Belos/Ifpack" << std::endl;
    allOK = runit(epetra, belos_ifpack) && allOK;

    Out::root() << "Running Aztec/ML" << std::endl;
    allOK = runit(epetra, aztec_ml) && allOK;

    Out::root() << "Running Aztec/Ifpack" << std::endl;
    allOK = runit(epetra, aztec_ifpack) && allOK;


    Out::root() << "Running BICGSTAB" << std::endl;
    allOK = runit(epetra, bicgstab) && allOK;

#ifdef BLAH
    if (nProc == 1)
    {
      Out::root() << "Running Amesos (serial)" << std::endl;
      allOK = runit(epetra, amesos) && allOK;
    }

    if (rank==0)
    {
      Out::root() << "Running dense LU (serial)" << std::endl;
      allOK = runit(serial, denseLU) && allOK;
    }
#endif
    allOK = globalAnd(allOK);

    if (allOK) 
    {
      Out::root() << "all Poisson solve tests PASSED!" << std::endl;
    }
    else
    {
      status = -1;
      Out::root() << "some Poisson solve tests FAILED!" << std::endl;
    }
  }
  catch(std::exception& e)
  {
    cout << "Caught exception: " << e.what() << std::endl;
    return -1;
  }
  return status;
}

