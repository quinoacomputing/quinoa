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


#include "SundanceChainRuleEvaluator.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceCombinatorialUtils.hpp"

using namespace Teuchos;
using std::cout;
using std::exception;
using Sundance::List;


void write(const MultipleDeriv& md,
           const Array<MultiSet<int> >& K,
           const Array<MultipleDeriv>& L);

int main(int argc, char** argv)
{
  typedef Array<OrderedPair<Array<MultiSet<int> >, Array<MultipleDeriv> > > CR;
  try
		{
      GlobalMPISession session(&argc, &argv);

			Expr u = new UnknownFunctionStub("u");
			Expr v = new UnknownFunctionStub("v");
			Expr w = new UnknownFunctionStub("w");

      int nArgs = 2;
      MultipleDeriv md = makeDeriv(u,u);
      int order = md.order();
      

      for (int l=1; l<=order; l++)
        {
          Array<int> s(l, nArgs);
          Array<Array<int> > distinctQ = indexCombinations(s);
          Set<MultiSet<int> > q;
          for (int p=0; p<distinctQ.size(); p++)
            {
              q.put(makeMultiSet(distinctQ[p]));
            }
          if (l > 1) cout << " + " << std::endl;
          for (Set<MultiSet<int> >::const_iterator 
                 i=q.begin(); i!=q.end(); i++)
            {
              const MultiSet<int>& lambda = *i;
              if (lambda != *(q.begin())) cout << " + " << std::endl;
              cout << "f_" << lambda << " * [";
              for (int s=1; s<=md.order(); s++)
                {
                  CR p = chainRuleTerms(s, lambda, md);
                  bool firstTerm = true;
                  for (CR::const_iterator j=p.begin(); j!=p.end(); j++)
                    {
                      if (!firstTerm) cout << "+";
                      firstTerm = false;
                      Array<MultiSet<int> > K = j->first();
                      Array<MultipleDeriv> L = j->second();
                      write(md, K, L);
                    }
                }
              cout << "]" << std::endl;
            }
        }
    }
	catch(std::exception& e)
		{
			Out::println(e.what());
		}
}


void write(const MultipleDeriv& md,
           const Array<MultiSet<int> >& K,
           const Array<MultipleDeriv>& L)
{
  int factor = chainRuleMultiplicity(md, K, L);
  if (factor != 1) cout << factor << "*";
  bool firstTerm = true;
  
  for (int j=0; j<K.size(); j++)
    {
      if (!firstTerm) cout << "*";
      firstTerm = false;
      bool firstFactor = true;
      for (MultiSet<int>::const_iterator k=K[j].begin(); k!=K[j].end(); k++)
        {
          if (!firstFactor) cout << "*";
          firstFactor = false;
          int q = *k;
          cout << "D[q_" << q << ", " << L[j] << "]";
        }
    }
}
