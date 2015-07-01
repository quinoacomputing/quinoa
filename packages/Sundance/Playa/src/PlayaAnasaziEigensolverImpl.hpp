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

#ifndef PLAYA_ANASAZIEIGENSOLVER_IMPL_HPP
#define PLAYA_ANASAZIEIGENSOLVER_IMPL_HPP

#include "PlayaDefs.hpp"
#include "PlayaAnasaziEigensolverDecl.hpp" 
#include "PlayaParameterListPreconditionerFactory.hpp" 
#include "PlayaPreconditionerFactory.hpp" 
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "PlayaAnasaziAdapter.hpp"


namespace Playa
{
using Teuchos::ParameterList;
using Anasazi::SimpleMV;
/** */
template <class MV, class OP> 
class InitTraits
{
public:
  /** */
  static RCP<OP> opPtr(const LinearOperator<double>& A);

  /** */
  static RCP<MV> makeMV(int numVecs, const VectorSpace<double>& space);

  /** */
  static Vector<double> vec(const RCP<MV>& mv, int i);
};


/** */
template <> class InitTraits<SimpleMV, LinearOperator<double> >
{
public:
  typedef SimpleMV            MV;
  typedef LinearOperator<double>            OP;

  /** */
  static RCP<OP> opPtr(const LinearOperator<double>& A)
    {
      if (A.ptr().get() != 0)
        return rcp(new LinearOperator<double>(A));
      else
      {
        RCP<LinearOperator<double> > rtn;
        return rtn;
      }
    }

  /** */
  static RCP<MV> makeMV(int blockSize, const VectorSpace<double>& space)
    {
      RCP<MV> mv = rcp(new MV(blockSize));
      for (int i=0; i<blockSize; i++) (*mv)[i] = space.createMember();
      return mv;
    }

  /** */
  static Vector<double> vec(const RCP<MV>& mv, int i)
    {
      return (*mv)[i];
    }

  
};




template <class Scalar>  
inline void AnasaziEigensolver<Scalar>::solve(
  const LinearOperator<Scalar>& K,
  const LinearOperator<Scalar>& M,
  Array<Vector<Scalar> >& evecs,
  Array<std::complex<Scalar> >& ew) const 
{
  typedef SimpleMV            MV;
  typedef LinearOperator<Scalar>            OP;

  typedef Anasazi::MultiVecTraits<Scalar,MV>     MVT;
  typedef Anasazi::OperatorTraits<Scalar,MV,OP>  OPT;

  TimeMonitor timer(solveTimer());
  VectorSpace<Scalar> KDomain = K.domain();

  RCP<OP> KPtr = InitTraits<MV, OP>::opPtr(K);
  RCP<OP> MPtr = InitTraits<MV, OP>::opPtr(M);



  
  // Eigensolver parameters
  std::string method = this->params().template get<string>("Method");
  int numEigs = this->params().template get<int>("Number of Eigenvalues");
  int blockSize = this->params().template get<int>("Block Size");
  bool usePrec = this->params().template get<bool>("Use Preconditioner");
  bool hermitian = this->params().template get<bool>("Is Hermitian");


  
  /* Make a multivector with row space = domain of K, column 
   * space = multiVec Space*/
  RCP<MV> mv = InitTraits<MV, OP>::makeMV(blockSize, KDomain);

  /* Fill the multivector with random values */
  MVT::MvRandom( *mv );

  /* Create a preconditioner */
  ParameterList precParams = this->params().sublist("Preconditioner");
  PreconditionerFactory<double> precFactory 
    = new ParameterListPreconditionerFactory(precParams);

  LinearOperator<Scalar> P;
  if (usePrec) 
  {
    TimeMonitor pTimer(precondBuildTimer());
    P = precFactory.createPreconditioner(K).right();
  }

  /* Create eigenproblem */
  RCP<Anasazi::Eigenproblem<Scalar,MV,OP> > problem;

  if (MPtr.get() != 0)
  {
    problem =
      rcp( new Anasazi::BasicEigenproblem<Scalar,MV,OP>(KPtr, MPtr, mv) );
  }
  else
  {
    problem =
      rcp( new Anasazi::BasicEigenproblem<Scalar,MV,OP>(KPtr, mv) );
  }

  ParameterList eigParams = this->params();
  problem->setHermitian(hermitian);
  problem->setNEV(numEigs);
  if (usePrec) problem->setPrec(InitTraits<MV, OP>::opPtr(P));

  bool ret = problem->setProblem();
  TEUCHOS_TEST_FOR_EXCEPTION(!ret, std::runtime_error,
    "Eigenproblem not setup correctly");
  

  // Create the solver manager
  RCP<Anasazi::SolverManager<Scalar,MV,OP> > MySolverMan;
  if (method=="Block Davidson")
  {
    MySolverMan = rcp(new Anasazi::BlockDavidsonSolMgr<Scalar,MV,OP>(problem, eigParams));
  }
  else if (method=="Block Krylov Schur")
  {
    MySolverMan = rcp(new Anasazi::BlockKrylovSchurSolMgr<Scalar,MV,OP>(problem, eigParams));
  }
  else if (method=="Simple LOBPCG")
  {
    MySolverMan = rcp(new Anasazi::SimpleLOBPCGSolMgr<Scalar,MV,OP>(problem, eigParams));
  }
  else if (method=="LOBPCG")
  {
    MySolverMan = rcp(new Anasazi::LOBPCGSolMgr<Scalar,MV,OP>(problem, eigParams));
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
      "solver method [" << method << "] not recognized");
  }

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMan->solve();
  Out::os() << "return code = " << returnCode << endl;
  TEUCHOS_TEST_FOR_EXCEPTION(returnCode != Anasazi::Converged, 
    std::runtime_error, "Anasazi did not converge!");
  
  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<Scalar,MV> sol = problem->getSolution();
  RCP<MV> evecs_mv = sol.Evecs;
  int numev = sol.numVecs;
  
  /* Copy the columns of the eigenvector MV into an array of Playa vectors */
  ew.resize(numev);
  evecs.resize(numev);

  for (int i=0; i<numev; i++)
  {
    Vector<Scalar> tmp = InitTraits<MV, OP>::vec(evecs_mv, i);

    evecs[i] = KDomain.createMember();
    evecs[i].acceptCopyOf(tmp);
    /* record the associated eigenvalue. The matrix is Hermitian so
     * we know the eigenvalue is real. */
    //evals[i] = sol.Evals[i].realpart;
    // if matrix might not be hermitian
    ew[i].real() = sol.Evals[i].realpart;
    ew[i].imag() = sol.Evals[i].imagpart;
  }
}


}


#endif
