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


#include "SundanceStochBlockJacobiSolver.hpp"
#include "Sundance.hpp"



#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#endif



namespace Sundance
{

  using std::ends;
  using std::setprecision;

void
StochBlockJacobiSolver::solve(const Array<LinearOperator<double> >& KBlock,
  const Array<Vector<double> >& fBlock,
  Array<Vector<double> >& xBlock) const
{
  Array<int> hasNonzeroMatrix(KBlock.size());
  for (int i=0; i<KBlock.size(); i++) hasNonzeroMatrix[i] = true;
  
  solve(KBlock, hasNonzeroMatrix, fBlock, xBlock);
}


void
StochBlockJacobiSolver::solve(const Array<LinearOperator<double> >& KBlock,
  const Array<int>& hasNonzeroMatrix,
  const Array<Vector<double> >& fBlock,
  Array<Vector<double> >& xBlock) const
{
  int L = KBlock.size();
  int P = pcBasis_.nterms();
  int Q = fBlock.size();

  /*
   * Solve the equations using block Gauss-Jacobi iteration
   */
  Array<Vector<double> > uPrev(P);
  Array<Vector<double> > uCur(P);

  for (int i=0; i<P; i++)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(fBlock[0].ptr().get()==0, 
      std::runtime_error, "empty RHS vector block i=[" << i << "]");
    uPrev[i] = fBlock[0].copy();
    uCur[i] = fBlock[0].copy();
    uPrev[i].zero();
    uCur[i].zero();
  }

  if (verbosity_) Out::root() << "starting Jacobi loop" << std::endl;
  bool converged = false;
  for (int iter=0; iter<maxIters_; iter++)
  {
    if (verbosity_) Out::root() << "Jacobi iter=" << iter << std::endl;
    bool haveNonConvergedBlock = false;
    double maxErr = -1.0;
    int numNonzeroBlocks = 0;
    for (int i=0; i<P; i++)
    {
      if (verbosity_) Out::root() << "Iter " << iter << ": block row i=" << i << " of " << P << " ..." << ends;
      Vector<double> b = fBlock[0].copy();
      b.zero();
      int nVecAdds = 0;
      for (int j=0; j<Q; j++)
      {
        double c_ij0 = pcBasis_.expectation(i,j,0);
        if (std::fabs(c_ij0) > 0.0) 
        {
          b = b + c_ij0 * fBlock[j];
          nVecAdds++;
        }
        if (j>=L) continue; 
        if (!hasNonzeroMatrix[j]) continue;
        Vector<double> tmp = fBlock[0].copy();
        tmp.zero();
        bool blockIsNeeded = false;
        for (int k=0; k<P; k++)
        {
          if (j==0 && k==i) continue;
          double c_ijk = pcBasis_.expectation(i,j,k);
          if (std::fabs(c_ijk) > 0.0)
          {
            tmp = tmp + c_ijk * uPrev[k];
            nVecAdds++;
            blockIsNeeded = true;
          }
        }
        numNonzeroBlocks += blockIsNeeded;
        b = (b - KBlock[j]*tmp);
        nVecAdds++;
      }
      b = b * (1.0/pcBasis_.expectation(i,i,0));
      if (verbosity_) Out::root() << "num vec adds = " << nVecAdds << std::endl;
      diagonalSolver_.solve(KBlock[0], b, uCur[i]);
      double err = (uCur[i]-uPrev[i]).norm2();
      if (err > convTol_) haveNonConvergedBlock=true;
      if (err > maxErr) maxErr = err;
    }

    /* update solution blocks */
    for (int i=0; i<P; i++) uPrev[i] = uCur[i].copy();
      
    /* done all block rows -- check convergence */
    if (!haveNonConvergedBlock)
    {
      if (verbosity_) Out::root() << "=======> max err=" << maxErr << std::endl;
      if (verbosity_) Out::root() << "=======> converged! woo-hoo!" << std::endl;
      if (verbosity_) Out::root() << "estimated storage cost: " 
                  << setprecision(3) << 100*((double) L)/((double) numNonzeroBlocks) 
                  << " percent of monolithic storage" << std::endl;
      converged = true;
      break;
    }
    else
    {
      if (verbosity_) Out::root() << "maxErr=" << maxErr << ", trying again" << std::endl;
    }
  }

  TEUCHOS_TEST_FOR_EXCEPT(!converged);
  xBlock = uCur;
}
}
