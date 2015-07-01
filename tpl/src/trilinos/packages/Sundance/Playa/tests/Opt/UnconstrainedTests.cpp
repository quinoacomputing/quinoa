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

#include "PlayaSteepestDescent.hpp"
#include "PlayaOptBuilder.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaOut.hpp"
#include "PlayaLinearCombinationImpl.hpp"

#include "Teuchos_GlobalMPISession.hpp"

#include "Teuchos_Time.hpp"
#include "PlayaMPIComm.hpp"

#include <fstream>

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaVectorImpl.hpp"
#include "PlayaBlockIteratorImpl.hpp"
#endif

using namespace Playa;
using std::ofstream;
using std::ostringstream;

class UncTestProb : public ObjectiveBase
{
public:
  /** */
  virtual ~UncTestProb() {}

  /** */
  virtual Vector<double> exactSoln() const = 0 ;
};

class Rosenbrock : public UncTestProb
{
public:
  Rosenbrock(int N, double alpha, const VectorType<double>& vecType)
    : N_(N), vs_(vecType.createEvenlyPartitionedSpace(MPIComm::self(), 2*N)),
      alpha_(alpha)
    {}

  void  evalGrad(const Vector<double>& x, double& f, 
    Vector<double>& grad) const ;

  void eval(const Vector<double>& x, double& f) const ;

  Vector<double> getInit() const ;

  Vector<double> exactSoln() const ;

  string description() const 
    {
      ostringstream oss;
      oss << "Rosenbrock[n=" << N_ << ", alpha=" << alpha_ << "]";
      return oss.str();
    }

private:
  int N_;
  VectorSpace<double> vs_;
  double alpha_;
};

void Rosenbrock::eval(const Vector<double>& x, double& f) const
{
  f = 0.0;
  for (int i=0; i<N_; i++)
  {
    double p = x[2*i+1] - x[2*i]*x[2*i];
    double q = 1.0-x[2*i];
    f += alpha_ * p*p + q*q;
  }
}

void Rosenbrock::evalGrad(const Vector<double>& x, double& f,
  Vector<double>& df) const
{
  f = 0.0;
  df.zero();
  for (int i=0; i<N_; i++)
  {
    double p = x[2*i+1] - x[2*i]*x[2*i];
    double q = 1.0-x[2*i];
    f += alpha_ * p*p + q*q;
    df[2*i] += -4.0*alpha_*p*x[2*i] - 2.0*q;
    df[2*i+1] += 2.0*alpha_*p; 
  }  
}

Vector<double> Rosenbrock::getInit() const
{
  Vector<double> rtn = vs_.createMember();
  rtn.setToConstant(-1.0);
  return rtn;
}

Vector<double> Rosenbrock::exactSoln() const
{
  Vector<double> rtn = vs_.createMember();
  rtn.setToConstant(1.0);
  return rtn;
}


class Ellipsoid : public UncTestProb
{
public:
  Ellipsoid(int n, const VectorType<double>& vecType)
    : n_(n), vs_(vecType.createEvenlyPartitionedSpace(MPIComm::self(), n))
    {}

  void  evalGrad(const Vector<double>& x, double& f, 
    Vector<double>& grad) const ;

  void eval(const Vector<double>& x, double& f) const ;

  Vector<double> getInit() const ;

  virtual Vector<double> exactSoln() const ;

  string description() const 
    {
      ostringstream oss;
      oss << "Ellipsoid[n=" << n_ << "]";
      return oss.str();
    }

private:
  int n_;
  VectorSpace<double> vs_;
};

void Ellipsoid::eval(const Vector<double>& x, double& f) const
{
  f = 0.0;
  for (int i=0; i<n_; i++)
  {
    f += ::sqrt(2+i)*x[i]*x[i];
  }
}

void Ellipsoid::evalGrad(const Vector<double>& x, double& f, 
  Vector<double>& grad) const
{
  f = 0.0;
  for (int i=0; i<n_; i++)
  {
    f += ::sqrt(i+2)*x[i]*x[i];
    grad[i] = 2.0*::sqrt(i+2)*x[i];
  }
}

Vector<double> Ellipsoid::getInit() const
{
  Vector<double> rtn = vs_.createMember();
  rtn.setToConstant(1.0);
  return rtn;
}

Vector<double> Ellipsoid::exactSoln() const
{
  Vector<double> rtn = vs_.createMember();
  rtn.setToConstant(0.0);
  return rtn;
}


int main(int argc, char *argv[])
{
  int stat = 0;
  try
  {
    GlobalMPISession session(&argc, &argv);
    Tabs::showDepth() = false;
    VectorType<double> vecType = new EpetraVectorType();
    double testTol = 1.0e-4;

    
    int n = 200;
    RCP<UncTestProb> ell = rcp(new Ellipsoid(n, vecType));
    RCP<UncTestProb> rosen = rcp(new Rosenbrock(2, 1.0, vecType));

    Array<RCP<UncTestProb> > probs = tuple(ell, rosen);
    Array<string> algs 
      = tuple<string>("basicLMBFGS", "steepestDescent");

    int numFail = 0;
    int testCount = 0;

    for (int i=0; i<probs.size(); i++)
    {
      RCP<UncTestProb> obj = probs[i];
      Out::root() << "====================================================="
                  << endl
                  << "  running test case: " << obj->description() << endl
                  << "====================================================="
                  << endl << endl;
      Tabs tab1;
      
      
      Vector<double> xInit = obj->getInit();
      Out::root() << tab1 << "Doing FD check of gradient..." << endl;
      bool fdOK = obj->fdCheck(xInit, 1.0e-6, 0);
      if (!fdOK)
      {
        Tabs tab2;
        Out::root() << tab2 << "FD check FAILED!" << endl;
        obj->fdCheck(xInit, 1.0e-6, 2);
        numFail++;
        continue;
      }
      else
      {
        Tabs tab2;
        Out::root() << tab2 << "FD check OK" << endl;
      }
      Out::root() << endl;

      for (int j=0; j<algs.size(); j++)
      {
        Tabs tab2;
        RCP<UnconstrainedOptimizerBase> opt 
          = OptBuilder::createOptimizer(algs[j] + ".xml");

        Out::root() << tab1 << "Running opt algorithm [" 
                    << algs[j] << "]" << endl;
        RCP<ConvergenceMonitor> mon = rcp(new ConvergenceMonitor());

        opt->setVerb(0);

        OptState state = opt->run(obj, xInit, mon);

        bool ok = true;

        if (state.status() != Opt_Converged)
        {
          Out::root() << tab2 << "[" << algs[j] 
                      << "] optimization failed: " << state.status() << endl;
          ok = false;
        }
        else
        {
          Out::root() << tab2 
                      << "[" << algs[j] 
                      << "] optimization converged!" << endl ;
          
          Vector<double> exactAns = obj->exactSoln();
          
          double locErr = (exactAns - state.xCur()).norm2();
          Out::root() << tab2 << "||x-x^*||=" << locErr << endl;
          if (locErr > testTol) 
          {
            ok = false;
          }
          
          string monName = algs[j] + ".dat";
          ofstream os(monName.c_str());
          mon->write(os);
        }
        
        if (ok) 
        {
          Out::root() << tab2 << " test  PASSED " << endl;
        }
        else
        {
          numFail++;
          Out::root() << tab2 << " ****** test FAILED ****** " << endl;
        }
        Out::root() << endl;
        testCount++;
      }
    }
    
    if (numFail > 0)
    {
      Out::root() << "detected " << numFail 
                  << " FAILURES out of " << testCount << " tests" << endl;
      stat = -1;
    }
    else
    {
      Out::root() << "all tests PASSED" << endl;
    }
  }
  catch(std::exception& e)
  {
    std::cerr << "Caught exception: " << e.what() << endl;
    stat = -1;
  }
  return stat;
}
