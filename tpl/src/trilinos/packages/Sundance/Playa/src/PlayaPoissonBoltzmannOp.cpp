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

#include "PlayaPoissonBoltzmannOp.hpp"
#include "PlayaTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using namespace Teuchos;

namespace Playa
{

PoissonBoltzmannOp::PoissonBoltzmannOp(int nLocal, const VectorType<double>& vecType)
  : NonlinearOperatorBase<double>(), J_(nLocal, vecType), importer_(),
    uLeftBC_(0.0), uRightBC_(2.0*log(cosh(1.0/sqrt(2.0))))
{
  setDomainAndRange(J_.domain(), J_.range());

  int rank = MPIComm::world().getRank();
  int nProc = MPIComm::world().getNProc();
  if (nProc > 1)
  {
    Array<int> ghosts;
    int low = J_.domain().baseGlobalNaturalIndex();
    int high = low + J_.domain().numLocalElements();
    if (rank != nProc - 1)
    {
      ghosts.append(high);
    }
    if (rank != 0) 
    {
      ghosts.append(low-1);
    }

    importer_ = vecType.createGhostImporter(J_.domain(), ghosts.size(), &(ghosts[0]));
  }
  else
  {
    importer_ = vecType.createGhostImporter(J_.domain(), 0, 0);
  }
}

Vector<double> PoissonBoltzmannOp::getInitialGuess() const
{
  Vector<double> rtn = J_.domain().createMember();

  rtn.setToConstant(0.5);

  return rtn;
}


LinearOperator<double> 
PoissonBoltzmannOp::computeJacobianAndFunction(Vector<double>& functionValue) const 
{
  Tabs tab;
  Out::root() << tab << "in PBOp::computeJacAndVec" << std::endl;
  Vector<double> x = currentEvalPt() ;
  Out::root() << tab << "eval pt = " << std::endl;
  Out::os() << x << std::endl;
  J_.setEvalPoint(x);
  

  RCP<GhostView<double> > u;
  Out::root() << tab << "importing view" << std::endl;
  importer_->importView(currentEvalPt(), u);
  Out::root() << tab << "done importing view" << std::endl;
  int low = J_.domain().baseGlobalNaturalIndex();
  int high = low + J_.domain().numLocalElements();
  Out::os() << tab << "my indices are: " << low << ", " << high-1 << std::endl;

  functionValue = J_.range().createMember();
  double h= J_.h();

  for (int r=low; r<high; r++)
  {
    Tabs tab1;
    double u_i = u->getElement(r);
    double f = 0.0;
    if (r==0) 
    {
      f = u_i - uLeftBC_;
    }
    else if (r==J_.domain().dim()-1)
    {
      f = u_i - uRightBC_;
    }
    else
    {
      double u_plus = u->getElement(r+1);
      double u_minus = u->getElement(r-1);
      f = (u_plus + u_minus - 2.0*u_i)/h/h - exp(-u_i);
    }
    functionValue[r-low] = f;
  }

  Out::root() << tab << "done PBOp::computeJacAndVec" << std::endl;
  return J_.getOp();
}

Vector<double> PoissonBoltzmannOp::exactSoln() const
{
  Tabs tab;
  Out::root() << tab << "in PBOp::exactSoln" << std::endl;
  Vector<double> rtn = J_.domain().createMember();

  int low = J_.domain().baseGlobalNaturalIndex();
  int high = low + J_.domain().numLocalElements();

  double h= J_.h();
  
  for (int r=low; r<high; r++)
  {
    double x = r*h;
    double u = 2.0*log(cosh(x/sqrt(2.0)));
    rtn[r-low] = u;
  }

  Out::os() << tab << "done PBOp::exactSoln" << std::endl;
  return rtn;
}

}
