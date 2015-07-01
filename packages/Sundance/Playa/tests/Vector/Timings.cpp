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

#include "PlayaGlobalAnd.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "PlayaEpetraVector.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaOut.hpp"
#include "PlayaRand.hpp"
#include <fstream>
#include "Teuchos_GlobalMPISession.hpp"

#include "Teuchos_Time.hpp"
#include "PlayaMPIComm.hpp"
using std::cout;
using std::endl;

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaBlockIteratorImpl.hpp"
#endif

using std::setw;
using std::ofstream;
using namespace Playa;
using namespace PlayaExprTemplates;

Array<Vector<double> > vecMaker(int nVecs, int n,
  int nProc, int rank, const VectorType<double>& vecType)
{
  /* This VS will go out of scope when the function is exited, but
   * its vectors will remember it */
  VectorSpace<double> space 
    = vecType.createEvenlyPartitionedSpace(MPIComm::world(), n);

  Rand::setLocalSeed(space.comm(), 314159);

  Array<Vector<double> > rtn(nVecs);
  for (int i=0; i<rtn.size(); i++)
  {
    rtn[i] = space.createMember();
    rtn[i].randomize();
  }
  return rtn;

}

bool runTest(int nProc, int n, int rank, const VectorType<double>& vecType,
  Array<Array<double> >& tList)
{
  bool pass = true;

  Array<Vector<double> > vecs = vecMaker(2, n, nProc, rank, vecType);
  Vector<double> x = vecs[0];
  Vector<double> y = vecs[1];

  
  Vector<double> z1 = x.copy();
  Vector<double> z2 = x.copy();
  
  Epetra_Vector& epZ = EpetraVector::getConcrete(z2);
  const Epetra_Vector& epY = EpetraVector::getConcrete(y);


  int nReps = 20;
  int nOutReps = 5;

  Time tEp("epetra time");
  Time tEx("expr time");
  Time tEpDot("epetra dot time");
  Time tExDot("expr dot time");

  double res1 = 0.0;  
  double res2 = 0.0;  

  for (int k=0; k<nOutReps; k++)
  {
    tEp.start();
    for (int i=0; i<nReps; i++)
    {
      epZ.Update(2.0, epY, 1.0);
    }
    tEp.stop();

    tEx.start();
    for (int i=0; i<nReps; i++)
    {
      z1 += 2.0*y;
    }
    tEx.stop();


    tEpDot.start();
    for (int i=0; i<nReps; i++)
    {
      epY.Dot(epY, &res1);
    }
    tEpDot.stop();

    tExDot.start();
    for (int i=0; i<nReps; i++)
    {
      res2 = y*y;
    }
    tExDot.stop();
  }

  tList.append(tuple(
                 (double) n,
                 tEx.totalElapsedTime()/nReps/nOutReps,
                 tEp.totalElapsedTime()/nReps/nOutReps,
                 tExDot.totalElapsedTime()/nReps/nOutReps,
                 tEpDot.totalElapsedTime()/nReps/nOutReps
      ));
      

  double updateErr = (z2 - z1).norm2();
  double dotErr = ::fabs(res1-res2);
  Out::root() << "update error (ex, ep) = " << updateErr << endl;
  Out::root() << "dot error (ex, ep) = " << dotErr << endl;
  
  pass = dotErr <= 1.0e-8 && dotErr <= 1.0e-8;
  return pass;
}


int main(int argc, char *argv[])
{
  int stat = 0;
  try
  {
    GlobalMPISession session(&argc, &argv);
    int nProc = session.getNProc();
    int rank = session.getRank();

    VectorType<double> type1 = new EpetraVectorType();

    bool allPass = true;

    Array<Array<double> > tList;
    

    int n = 1;
    int i=0;
    while (n*nProc < 40000)
    {
      Out::root() << "running n=" << n << endl;
      allPass = runTest(nProc, n, rank, type1, tList) && allPass;
      Out::root() << tList[i] << endl;
      n = (int) ceil(1.8 * n);
      i++;
    }

    ofstream of("playa-timings.dat");
    for (int i=0; i<tList.size(); i++)
    {
      for (int j=0; j<tList[i].size(); j++)
      {
        if (j!=0) cout << " ";
        of << setw(14) << tList[i][j];
      }
      of << endl;
    }

    allPass = globalAnd(allPass);

    if (!allPass)
    {
      Out::root() << "detected a test that FAILED" << std::endl;
      stat = -1;
    }
    else
    {
      Out::root() << "all tests PASSED" << std::endl;
    }

    TimeMonitor::summarize();
  }
  catch(std::exception& e)
  {
    std::cerr << "Caught exception: " << e.what() << std::endl;
    stat = -1;
  }
  return stat;
}
