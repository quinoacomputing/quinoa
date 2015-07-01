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


#include "SundanceCombinatorialUtils.hpp"
#include "SundanceIntVec.hpp"
#include "Teuchos_GlobalMPISession.hpp"


using namespace Sundance;
using namespace Teuchos;



#define TEST_MS(x) \
  {\
    Set<MultiSet<int> > subs = multisetSubsets(x);\
    write(x, subs);                                                \
    Set<MultiSet<MultiSet<int> > > parts = multisetPartitions(x);\
    write(x, parts);                                                 \
    Array<Array<MultiSet<int> > > comps = multisetCompositions(x);\
    write(x, comps);                                               \
  }



void write(const MultiSet<int>& x, 
           const Set<MultiSet<MultiSet<int> > >& y)
{
  std::cout << "---- Partitions of " << x << " ------------------"
       << std::endl;
  for (Set<MultiSet<MultiSet<int> > >::const_iterator 
         i=y.begin(); i!=y.end(); i++)
    {
      std::cout << *i << std::endl;
    }
}

void write(const MultiSet<int>& x, 
           const Array<Array<MultiSet<int> > >& y)
{
  std::cout << "---- Compositions of " << x << " ------------------"
       << std::endl;
  for (int i=0; i<y.size(); i++)
    {
      std::cout << y[i] << std::endl;
    }
}

void write(const MultiSet<int>& x, 
           const Set<MultiSet<int> >& y)
{
  std::cout << "---- Subsets of " << x << " ------------------"
       << std::endl;
  for (Set<MultiSet<int> >::const_iterator 
         i=y.begin(); i!=y.end(); i++)
    {
      std::cout << *i << std::endl;
    }
}



int main(int argc, char** argv)
{
  int stat = 0;
  try
		{
      GlobalMPISession session(&argc, &argv);


      bool bad = false;

      for (int n=1; n<=4; n++)
        {
          Array<Array<Array<int> > > c = compositions(n);
          std::cout << "N=" << n << " compositions=" << c << std::endl;

          MultiSet<int> mu;
          for (int m=1; m<=n; m++)
            {
              mu.put(m);

            }
          for (int m=1; m<=n; m++)
            {
              Array<Array<Array<int> > > b = binnings(mu, m);
              std::cout << "binnings = " << b << std::endl;
            }

          std::cout << "--------- non-neg compositions" << std::endl;
          for (int m=1; m<=n; m++)
            {
              for (int k=1; k<=n; k++)
                {
                  Array<Array<int> > a = nonNegCompositions(m, k);
                  std::cout << m << " " << k << " " << std::endl;
                  for (int l=0; l<a.size(); l++)
                    {
                      std::cout << "         " << a[l] << std::endl;
                    }
                }
            }
          
          std::cout << "-------- index combs ---- " << std::endl;
          Array<int> s = tuple(2,3,2);
          Array<Array<int> > C = indexCombinations(s);
          for (int m=0; m<C.size(); m++)
            {
              std::cout << C[m] << std::endl;
            }
        }

      std::cout << "--------- index tuples ----------------" << std::endl;

      Array<Array<int> > x = distinctIndexTuples(2, 6);

      std::cout << "num choices = " << x.size() << std::endl;

      for (int i=0; i<x.size(); i++) 
        {
          if ((i % 5)==0) std::cout << std::endl;
          std::cout << x[i] << std::endl;
        }

      std::cout << "--------- int vec parts ----------------" << std::endl;
      IntVec iv = intVec(1,3,2);
      Array<Array<IntVec> > parts;
      iv.getPartitions(3, parts);
      std::cout << "----------- partitions of " << iv << " ---------" 
                << std::endl;
      for (int i=0; i<parts.size(); i++)
      {
        std::cout << i << " -- " << std::endl;
        for (int j=0; j<parts[i].size(); j++) 
        {
          std::cout << "\t\t" << parts[i][j] << std::endl;
        }
      }

      std::cout << "--------- weighted parts ----------------" << std::endl;
      Array<int> wgts = tuple(3,1);
      Array<Array<int> > wParts;
      int M = 3;
      weightedPartitions(M, wgts, wParts);
      std::cout << "----------- partitions of " << M << " ---------" 
                << std::endl;
      for (int i=0; i<wParts.size(); i++)
      {
        std::cout << "\t" << wParts[i] << std::endl;
      }

      std::cout << "--------- weighted ordered parts ---------" << std::endl;

      IntVec iv2 = intVec(1,3,2,1);
      Array<Array<IntVec> > vParts;

      weightedOrderedPartitions(iv2, wgts, vParts);
      std::cout << "----------- partitions of " << iv2 << " ---------" 
                << std::endl;
      for (int i=0; i<vParts.size(); i++)
      {
        std::cout << "\t" << vParts[i] << std::endl;
      }
      
      IntVec lam = intVec(1,2,1);
      IntVec nu = intVec(2,3);
      Array<Array<IntVec> > K;
      Array<Array<IntVec> > L;
      int S = 2;
      pSet(lam, nu, S, K, L);

      std::cout << "pSet(" << lam << ", " << nu << ", " << S 
                << ")" << std::endl;
      for (int i=0; i<K.size(); i++)
      {
        std::cout << "i=" << i << std::endl;
        IntVec kSum(K[i][0].size());
        IntVec lSum(L[i][0].size());
        
        for (int j=0; j<K[i].size(); j++)
        {
          kSum = kSum + K[i][j];
          lSum = lSum + K[i][j].abs() * L[i][j];
          std::cout << "\t\t\t K=" << K[i][j] << ", L="
                    << L[i][j] << std::endl;
        }
        std::cout << "\t\tKSum=" << kSum << std::endl;
        std::cout << "\t\tLSum=" << lSum << std::endl;
      }
      

#ifdef BLAH
      TEST_MS(makeMultiSet(1));
      TEST_MS(makeMultiSet(1, 1));
      TEST_MS(makeMultiSet(1, 2));
      TEST_MS(makeMultiSet(1, 1, 2));
      TEST_MS(makeMultiSet(1, 1, 2, 2));
      TEST_MS(makeMultiSet(1, 2, 3));
      TEST_MS(makeMultiSet(1, 2, 2, 3));
      TEST_MS(makeMultiSet(1, 2, 2, 3, 3));
      TEST_MS(makeMultiSet(1, 2, 2, 3, 3, 3));
#endif // BLAH

      if (!bad) 
        {
          std::cerr << "all tests PASSED" << std::endl;
        }
      else
        {
          stat = -1;
          std::cerr << "a test has FAILED" << std::endl;
        }
    }
	catch(std::exception& e)
		{
      stat = -1;
      std::cerr << "detected exception " << e.what() << std::endl;
		}

  return stat;
}
