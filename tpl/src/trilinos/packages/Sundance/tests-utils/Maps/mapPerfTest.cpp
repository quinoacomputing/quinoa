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


#include "SundanceIntHashSet.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Time.hpp"

using namespace Sundance;
using namespace Teuchos;

int main(int argc, char** argv)
{
  try
    {
      GlobalMPISession session(&argc, &argv);
      
      int nRow = 40000;
      int nReps = 4;
      
      Time tSet("STL set time");
      Time tHash("hash set time");
      
      for (int nData=20; nData<360; nData+=20)
      {
        {
        tSet.start();
        Array<Set<int> > sets(nRow);


        for (int r=0; r<nRow; r++)
          {
            Set<int>& s = sets[r];
            for (int rep=0; rep<nReps; rep++)
              {
                for (int x=0; x<nData; x++)
                  {
                    int y = rand() % nData;
                    s.put(y);
                  }
              }
          }
        tSet.stop();
        }
        {
          tHash.start();
          Array<IntHashSet> sets(nRow+1);
          for (int r=0; r<nRow; r++)
            {
              sets[r].setCapacity(nData+1);
            }


          for (int r=0; r<nRow; r++)
            {
            IntHashSet& s = sets[r];
            for (int rep=0; rep<nReps; rep++)
              {
                for (int x=0; x<nData; x++)
                  {
                    int y = rand() % nData;
                    s.put(y);
                  }
              }
          }
        tHash.stop();
        }

        
        std::cerr << nData << "\t set=" << tSet.totalElapsedTime()
             << "\t hash=" << tHash.totalElapsedTime() 
             << "\t ratio=" 
             << tHash.totalElapsedTime()/tSet.totalElapsedTime()
             << std::endl;
      }

      
    }
  catch(std::exception& e)
    {
      std::cerr << "caught exception " << e.what() << std::endl;
    }
  
}
